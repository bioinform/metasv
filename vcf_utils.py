import logging

logger = logging.getLogger(__name__)

import sys
import os
import argparse
import subprocess
import pysam
import bisect
from fasta_utils import *
from sv_interval import *

def print_header(header, file_fd):
  for line in header:
    file_fd.write("%s\n" % (line))

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
  info_dict = {}
  for field in info_fields:
    if field.find("=") >= 0:
      info_key, info_value = field.split("=")
      info_dict[info_key] = info_value.split(",")
    else:
      info_dict[field] = True
  return info_dict

def load_gap_intervals(gap_file):
  if gap_file is None: return []
  logger.info("Loading the gaps in the genome from %s" % (gap_file))
  with open(gap_file) as gap_file_fd:
    gap_intervals = [SVInterval(it.contig, it.start, it.end, it.name, "gap") for it in pysam.tabix_file_iterator(gap_file_fd, parser=pysam.asBed())]
  return merge_intervals(gap_intervals)

def get_gt(gt, fmt):
  fmt_index = [i for i, field in enumerate(fmt.split(":")) if field == "GT"][0]
  return gt.split(":")[fmt_index]

def load_intervals(in_vcf, intervals={}, gap_intervals=[], include_intervals=[], source=None, contig_whitelist=[], is_gatk=False):
  if not os.path.isfile(in_vcf): return intervals
  logger.info("Loading SV intervals from %s" % (in_vcf))
  tabix_file = pysam.Tabixfile(in_vcf, parser=pysam.asVCF())

  for vcf_record in (record for record in tabix_file.fetch() if record.contig in contig_whitelist):

    fmt = "GT"
    gt = "./1"
    try:
      gt = get_gt(vcf_record[0], vcf_record.format)
    except IndexError:
      pass

    info_dict = parse_info(vcf_record.info)

    if is_gatk:
      if vcf_record.filter != "PASS" and vcf_record.filter != ".": continue
      if vcf_record.alt.find(",") >= 0: continue
      if len(vcf_record.ref) != 1 and len(vcf_record.alt) != 1: continue
      if len(vcf_record.ref) < 50 and len(vcf_record.alt) < 50: continue
      if len(vcf_record.ref) == 1:
        wiggle = 100 if source in ["Pindel", "BreakSeq", "HaplotypeCaller"] else 0
        interval = SVInterval(vcf_record.contig, vcf_record.pos + 1, vcf_record.pos + 1, source, "INS", len(vcf_record.alt)-1, sources=set([source]), wiggle=wiggle, gt=gt)
      else:
        interval = SVInterval(vcf_record.contig, vcf_record.pos, vcf_record.pos + len(vcf_record.ref) - 1, source, "DEL", len(vcf_record.ref)-1, sources=set([source]), gt=gt)
    else:
      if source == "BreakSeq" and vcf_record.filter != "PASS": continue
      if vcf_record.alt.find(",") != -1: continue
      #logger.info(str(vcf_record.alt))
      sv_type = info_dict["SVTYPE"][0]
      if sv_type == "DUP:TANDEM": sv_type = "DUP"
      if "SVLEN" not in info_dict:
        if source == "BreakSeq" and sv_type == "INS": 
          info_dict["SVLEN"] = [0]
        else:
          continue
      svlen = abs(int(info_dict["SVLEN"][0]))
      if svlen < 50: continue
      wiggle = 100 if (source in ["Pindel", "BreakSeq", "HaplotypeCaller"] and sv_type == "INS") else 0
      if source == "Pindel" and sv_type == "INS": vcf_record.pos += 1
      interval = SVInterval(vcf_record.contig, vcf_record.pos, int(info_dict["END"][0]), source, sv_type, svlen, sources=set([source]), wiggle=wiggle, gt=gt)
    if interval_overlaps_interval_list(interval, gap_intervals):
      logger.warn("Skipping " + str(interval) + " due to overlap with gaps")
      continue
    if not interval_overlaps_interval_list(interval, include_intervals, min_fraction_self = 1.0):
      logger.warn("Skipping " + str(interval) + " due to being outside the include regions")
      continue
    interval.info_dict = info_dict

    if interval.sv_type not in intervals:
      intervals[interval.sv_type] = [interval]
    else:
      intervals[interval.sv_type].append(interval)
  return intervals

def print_vcf_header(outfd, reference, contigs, sample):
  contig_str = ""
  for contig in contigs:
    contig_str += "##contig=<ID=%s,length=%d>\n" % (contig.name, contig.length)

  outfd.write("""##fileformat=VCFv4.1
##reference=%s
##source=%s
##INFO=<ID=BKPTID,Number=.,Type=String,Description=\"ID of the assembled alternate allele in the assembly file\">
##INFO=<ID=CIEND,Number=2,Type=Integer,Description=\"Confidence interval around END for imprecise variants\">
##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">
##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Imprecise structural variation\">
##INFO=<ID=SVCNV,Number=0,Type=Flag,Description=\"Structural variation or copy number variation\">
##INFO=<ID=MEINFO,Number=4,Type=String,Description=\"Mobile element info of the form NAME,START,END,POLARITY\">
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">
##INFO=<ID=NORMAL_COUNT,Number=1,Type=Integer,Description=\"Number of normal reads supporting reference\">
##INFO=<ID=NUM_READS,Number=1,Type=Integer,Description=\"Number of reads supporting event\">
##INFO=<ID=SCORE,Number=1,Type=Integer,Description=\"Score reported by tool\">
##INFO=<ID=CIPOS,Number=2,Type=Integer,Description=\"Confidence interval around POS for imprecise variants\">
##INFO=<ID=CIEND,Number=2,Type=Integer,Description=\"Confidence interval around END for imprecise variants\">
##INFO=<ID=HOMLEN,Number=.,Type=Integer,Description=\"Length of base pair identical micro-homology at event breakpoints\">
##INFO=<ID=HOMSEQ,Number=.,Type=String,Description=\"Sequence of base pair identical micro-homology at event breakpoints\">
##INFO=<ID=METRANS,Number=4,Type=String,Description=\"Mobile element transduction info of the form CHR,START,END,POLARITY\">
##INFO=<ID=DGVID,Number=1,Type=String,Description=\"ID of this element in Database of Genomic Variation\">
##INFO=<ID=DBVARID,Number=1,Type=String,Description=\"ID of this element in DBVAR\">
##INFO=<ID=DBRIPID,Number=1,Type=String,Description=\"ID of this element in DBRIP\">
##INFO=<ID=MATEID,Number=.,Type=String,Description=\"ID of mate breakends\">
##INFO=<ID=PARID,Number=1,Type=String,Description=\"ID of partner breakend\">
##INFO=<ID=EVENT,Number=1,Type=String,Description=\"ID of event associated to breakend\">
##INFO=<ID=CILEN,Number=2,Type=Integer,Description=\"Confidence interval around the length of the inserted material between breakends\">
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth of segment containing breakend\">
##INFO=<ID=DPADJ,Number=.,Type=Integer,Description=\"Read Depth of adjacency\">
##INFO=<ID=CN,Number=1,Type=Integer,Description=\"Copy number of segment containing breakend\">
##INFO=<ID=CNADJ,Number=.,Type=Integer,Description=\"Copy number of adjacency\">
##INFO=<ID=CICN,Number=2,Type=Integer,Description=\"Confidence interval around copy number for the segment\">
##INFO=<ID=CICNADJ,Number=.,Type=Integer,Description=\"Confidence interval around copy number for the adjacency\">
##INFO=<ID=SVTOOL,Number=1,Type=String,Description=\"Tool used to generate the SV\">
##INFO=<ID=SOURCES,Number=.,Type=String,Description=\"List of original raw SV calls as Toolname:Start:End:Size\">
##INFO=<ID=NUM_SVMETHODS,Number=1,Type=Integer,Description=\"Number of methods supporting the event\">
##INFO=<ID=VT,Number=1,Type=String,Description=\"indicates what type of variant the line represents\">
##INFO=<ID=SVMETHOD,Number=.,Type=String,Description=\"Type of approach used to detect SV: RP (read pair), RD (read depth), SR (split read), JM (junction) or AS (assembly)\">
##INFO=<ID=natorRD,Number=1,Type=Float,Description=\"CNVnator: Normalized RD\">
##INFO=<ID=natorP1,Number=1,Type=Float,Description=\"CNVnator: e-val by t-test\">
##INFO=<ID=natorP2,Number=1,Type=Float,Description=\"CNVnator: e-val by Gaussian tail\">
##INFO=<ID=natorP3,Number=1,Type=Float,Description=\"CNVnator: e-val by t-test (middle)\">
##INFO=<ID=natorP4,Number=1,Type=Float,Description=\"CNVnator: e-val by Gaussian tail (middle)\">
##INFO=<ID=natorQ0,Number=1,Type=Float,Description=\"CNVnator: Fraction of reads with 0 mapping quality\">
##FILTER=<ID=LowQual,Description=\"Low Quality\">
%s##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
##FORMAT=<ID=CN,Number=1,Type=Integer,Description=\"Copy number genotype for imprecise events\">
##FORMAT=<ID=CNQ,Number=1,Type=Float,Description=\"Copy number genotype quality for imprecise events\">
##FORMAT=<ID=CNL,Number=.,Type=Float,Description=\"Copy number genotype likelihood for imprecise events\">
##FORMAT=<ID=NQ,Number=1,Type=Integer,Description=\"Phred style probability score that the variant is novel with respect to the genome's ancestor\">
##FORMAT=<ID=HAP,Number=1,Type=Integer,Description=\"Unique haplotype identifier\">
##FORMAT=<ID=AHAP,Number=1,Type=Integer,Description=\"Unique identifier of ancestral haplotype\">
##ALT=<ID=DEL,Description=\"Deletion\">
##ALT=<ID=DEL:ME:ALU,Description=\"Deletion of ALU element\">
##ALT=<ID=DEL:ME:L1,Description=\"Deletion of L1 element\">
##ALT=<ID=DUP,Description=\"Duplication\">
##ALT=<ID=DUP:TANDEM,Description=\"Tandem Duplication\">
##ALT=<ID=INS,Description=\"Insertion of novel sequence\">
##ALT=<ID=INS:ME:ALU,Description=\"Insertion of ALU element\">
##ALT=<ID=INS:ME:L1,Description=\"Insertion of L1 element\">
##ALT=<ID=INV,Description=\"Inversion\">
##ALT=<ID=CNV,Description=\"Copy number variable region\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n""" % (reference, " ".join(sys.argv), contig_str, sample))
