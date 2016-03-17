#!/usr/bin/env python

import argparse
import logging
from metasv.generate_final_vcf import convert_metasv_bed_to_vcf

if __name__ == "__main__":
    FORMAT = '%(levelname)s %(asctime)-15s %(name)-20s %(message)s'
    logging.basicConfig(level=logging.INFO, format=FORMAT)
    logger = logging.getLogger(__name__)

    parser = argparse.ArgumentParser(
        description="Convert MetaSV final BED to VCF",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("--sample", help="Sample name", required=True)
    parser.add_argument("--bed", help="MetaSV final BED", required=True)
    parser.add_argument("--vcf", help="Final VCF to output", required=True)
    parser.add_argument("--reference", help="Reference FASTA")
    parser.add_argument("--work", help="Work directory", default="work")
    parser.add_argument("--pass_only", action="store_true",
                        help="Output only PASS calls")

    args = parser.parse_args()

    convert_metasv_bed_to_vcf(bedfile=args.bed, vcf_out=args.vcf,
                              workdir=args.work,
                              sample=args.sample,
                              reference=args.reference,
                              pass_calls=args.pass_only)