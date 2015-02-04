#!/usr/bin/env python

import logging
import sys
import argparse
import vcf
import pysam

from metasv.breakdancer_reader import BreakDancerReader
from metasv.pindel_reader import PindelReader
from metasv.cnvnator_reader import CNVnatorReader
from metasv.breakseq_reader import BreakSeqReader
from metasv.vcf_utils import get_template
from metasv.fasta_utils import get_contigs

logging.basicConfig()
logger = logging.getLogger(__name__)

tool_to_reader = {"BreakDancer": BreakDancerReader, "Pindel": PindelReader, "CNVnator": CNVnatorReader, "BreakSeq": BreakSeqReader}


def convert_svtool_to_vcf(file_name, sample, out_vcf, toolname, reference, sort=False, index=False):
    vcf_template_reader = get_template()
    vcf_template_reader.samples = [sample]

    vcf_fd = open(out_vcf, "w") if out_vcf is not None else sys.stdout
    vcf_writer = vcf.Writer(vcf_fd, vcf_template_reader)

    reference_handle = pysam.Fastafile(reference) if reference else None
    reference_contigs = get_contigs(reference)
    if sort:
        if not reference_contigs:
            logger.warn("Chromosomes will be sorted in lexicographic order since reference is missing")
        else:
            logger.info("Chromosomes will be sorted by the reference order")

    vcf_records = []
    for tool_record in tool_to_reader[toolname](file_name, reference_handle = reference_handle):
        vcf_record = tool_record.to_vcf_record(sample)
        if vcf_record is None:
            continue
        if sort:
            vcf_records.append(vcf_record)
        else:
            vcf_writer.write_record(vcf_record)

    if sort:
        if reference_contigs:
            contigs_order_dict = {contig.name: index for (index, contig) in enumerate(reference_contigs)}
            vcf_records.sort(cmp=lambda x, y: cmp((contigs_order_dict[x.CHROM], x.POS), (contigs_order_dict[y.CHROM], y.POS)))
        else:
            vcf_records.sort(cmp=lambda x, y: cmp((x.CHROM, x.POS), (y.CHROM, y.POS)))

        for vcf_record in vcf_records:
            vcf_writer.write_record(vcf_record)
    vcf_writer.close()
    if out_vcf and index:
        pysam.tabix_index(out_vcf, force=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert SV tool output file to VCF",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--input", help="SV tool output file", required=False)
    parser.add_argument("--output", help="Output VCF to create", required=False)
    parser.add_argument("--tool", help="Tool name", required=False, default="BreakDancer",
                        choices=sorted(tool_to_reader.keys()))
    parser.add_argument("--sample", help="Sample name", required=True)
    parser.add_argument("--reference", help = "Reference FASTA")
    parser.add_argument("--sort", action = "store_true", help = "Sort the VCF records before writing")
    parser.add_argument("--index", action="store_true", help="Tabix compress and index the output VCF")

    args = parser.parse_args()
    convert_svtool_to_vcf(args.input, args.sample, args.output, args.tool, args.reference, sort=args.sort, index=args.index)
