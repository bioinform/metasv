import logging

logger = logging.getLogger(__name__)

import sys
import os
import argparse
import subprocess
import pysam
import bisect
import vcf
from fasta_utils import *
from sv_interval import *


def print_header(header, file_fd):
    for line in header:
        file_fd.write("%s\n" % (line))

def get_template():
    return vcf.Reader(open(os.path.join(os.path.dirname(os.path.realpath(__file__)), "resources/template.vcf")))

def merge_vcfs(in_vcfs_dir, contigs, out_vcf):
    logger.info("Mergings per-chromosome VCFs from %s" % (in_vcfs_dir))
    header_done = False
    out_vcf_file = open(out_vcf, "w")
    for contig in contigs:
        chr_vcf = os.path.join(in_vcfs_dir, "%s.vcf.gz" % (contig.name))
        if os.path.isfile(chr_vcf):
            chr_tabix_file = pysam.Tabixfile(chr_vcf)
            if not header_done:
                print_header(chr_tabix_file.header, out_vcf_file)
            for entry in chr_tabix_file.fetch():
                out_vcf_file.write("%s\n" % (entry))
            chr_tabix_file.close()
    out_vcf_file.close()
    pysam.tabix_index(out_vcf, force=True, preset="vcf")


def parse_info(info):
    info_fields = info.split(";")
    info = {}
    for field in info_fields:
        if field.find("=") >= 0:
            info_key, info_value = field.split("=")
            info[info_key] = info_value.split(",")
        else:
            info[field] = True
    return info


def load_gap_intervals(gap_file):
    if gap_file is None: return []
    logger.info("Loading the gaps in the genome from %s" % (gap_file))
    with open(gap_file) as gap_file_fd:
        gap_intervals = [SVInterval(it.contig, it.start, it.end, it.name, "gap") for it in
                         pysam.tabix_file_iterator(gap_file_fd, parser=pysam.asBed())]
    return merge_intervals(gap_intervals)


def get_gt(gt, fmt):
    fmt_index = [i for i, field in enumerate(fmt.split(":")) if field == "GT"][0]
    return gt.split(":")[fmt_index]


def load_intervals(in_vcf, intervals={}, gap_intervals=[], include_intervals=[], source=None, contig_whitelist=[],
                   minsvlen=50, wiggle=100, inswiggle=100):
    if not os.path.isfile(in_vcf): return intervals
    logger.info("Loading SV intervals from %s" % (in_vcf))

    vcf_reader = vcf.Reader(open(in_vcf))
    # Assume single sample for now
    sample = vcf_reader.samples[0]
    for vcf_record in vcf_reader:
        if vcf_record.CHROM not in contig_whitelist: continue

        gt = vcf_record.genotype(sample)
        print vcf_record.FILTER, vcf_record.ALT, vcf_record.REF, vcf_record.INFO, vcf_reader.samples

        if source in ["HaplotypeCaller"]:
            if vcf_record.FILTER and "PASS" not in vcf_record.FILTER: continue

            # Ignore tri-allelic stuff
            if len(vcf_record.ALT) > 1: continue

            # Ignore MNPs
            if len(vcf_record.REF) != 1 and len(vcf_record.ALT[0]) != 1: continue

            # Check for SV length
            if len(vcf_record.REF) < minsvlen and len(vcf_record.ALT[0]) < minsvlen: continue

            if not vcf_record.is_indel: continue

            if vcf_record.is_deletion:
                interval = SVInterval(vcf_record.contig,
                                      vcf_record.pos,
                                      vcf_record.pos + len(vcf_record.ref) - 1,
                                      source,
                                      "DEL",
                                      len(vcf_record.ref) - 1,
                                      sources=set([source]),
                                      wiggle=wiggle,
                                      gt=gt)
            else:
                interval = SVInterval(vcf_record.contig,
                                      vcf_record.pos + 1,
                                      vcf_record.pos + 1,
                                      source,
                                      "INS",
                                      len(vcf_record.alt) - 1,
                                      sources=set([source]),
                                      wiggle=max(inswiggle,wiggle),
                                      gt=gt)

        else:
            if source == "BreakSeq" and "PASS" not in vcf_record.FILTER: continue
            if len(vcf_record.ALT) > 1: continue
            sv_type = vcf_record.INFO["SVTYPE"]
            if sv_type == "DUP:TANDEM": sv_type = "DUP"
            if "SVLEN" not in vcf_record.INFO:
                if source == "BreakSeq" and sv_type == "INS":
                    vcf_record.INFO["SVLEN"] = 0
                else:
                    continue

            svlen = abs(vcf_record.INFO["SVLEN"])

            if svlen < minsvlen: continue
            wiggle = max(inswiggle,wiggle) if (source in ["Pindel", "BreakSeq", "HaplotypeCaller"] and sv_type == "INS") else wiggle
            if source == "Pindel" and sv_type == "INS": vcf_record.pos += 1
            interval = SVInterval(vcf_record.contig, vcf_record.pos, int(vcf_record.INFO["END"]), source, sv_type, svlen,
                                  sources=set([source]), wiggle=wiggle, gt=gt)
        if interval_overlaps_interval_list(interval, gap_intervals):
            logger.warn("Skipping " + str(interval) + " due to overlap with gaps")
            continue
        if not interval_overlaps_interval_list(interval, include_intervals, min_fraction_self=1.0):
            logger.warn("Skipping " + str(interval) + " due to being outside the include regions")
            continue
        interval.info = copy.deepcopy(vcf_record.INFO)

        if interval.sv_type not in intervals:
            intervals[interval.sv_type] = [interval]
        else:
            intervals[interval.sv_type].append(interval)
    return intervals

"""
    tabix_file = pysam.Tabixfile(in_vcf, parser=pysam.asVCF())

    for vcf_record in (record for record in tabix_file.fetch() if record.contig in contig_whitelist):

        fmt = "GT"
        gt = "./1"
        try:
            gt = get_gt(vcf_record[0], vcf_record.format)
        except IndexError:
            pass

        info = parse_info(vcf_record.info)

        if source in ["HaplotypeCaller"]:
            if vcf_record.filter != "PASS" and vcf_record.filter != ".": continue

            # ignore tri-allelic stuff
            if vcf_record.alt.find(",") >= 0: continue

            # ignore MNPs
            if len(vcf_record.ref) != 1 and len(vcf_record.alt) != 1: continue

            # Check for SV length
            if len(vcf_record.ref) < minsvlen and len(vcf_record.alt) < minsvlen: continue

            if len(vcf_record.ref) == 1:
                # This is the case for insertions
                interval = SVInterval(vcf_record.contig,
                                      vcf_record.pos + 1,
                                      vcf_record.pos + 1,
                                      source,
                                      "INS",
                                      len(vcf_record.alt) - 1,
                                      sources=set([source]),
                                      wiggle=max(inswiggle,wiggle),
                                      gt=gt)
            else:
                interval = SVInterval(vcf_record.contig,
                                      vcf_record.pos,
                                      vcf_record.pos + len(vcf_record.ref) - 1,
                                      source,
                                      "DEL",
                                      len(vcf_record.ref) - 1,
                                      sources=set([source]),
                                      wiggle=wiggle,
                                      gt=gt)
        else:
            if source == "BreakSeq" and vcf_record.filter != "PASS": continue
            if vcf_record.alt.find(",") != -1: continue
            # logger.info(str(vcf_record.alt))
            sv_type = info["SVTYPE"][0]
            if sv_type == "DUP:TANDEM": sv_type = "DUP"
            if "SVLEN" not in info:
                if source == "BreakSeq" and sv_type == "INS":
                    info["SVLEN"] = [0]
                else:
                    continue
            svlen = abs(int(info["SVLEN"][0]))
            if svlen < minsvlen: continue
            wiggle = max(inswiggle,wiggle) if (source in ["Pindel", "BreakSeq", "HaplotypeCaller"] and sv_type == "INS") else wiggle
            if source == "Pindel" and sv_type == "INS": vcf_record.pos += 1
            interval = SVInterval(vcf_record.contig, vcf_record.pos, int(info["END"][0]), source, sv_type, svlen,
                                  sources=set([source]), wiggle=wiggle, gt=gt)
        if interval_overlaps_interval_list(interval, gap_intervals):
            logger.warn("Skipping " + str(interval) + " due to overlap with gaps")
            continue
        if not interval_overlaps_interval_list(interval, include_intervals, min_fraction_self=1.0):
            logger.warn("Skipping " + str(interval) + " due to being outside the include regions")
            continue
        interval.info = info

        if interval.sv_type not in intervals:
            intervals[interval.sv_type] = [interval]
        else:
            intervals[interval.sv_type].append(interval)
    return intervals
    """
