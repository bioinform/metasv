import logging
import sys
import argparse
import os

import vcf


logger = logging.getLogger(__name__)

mydir = os.path.dirname(os.path.realpath(__file__))

from breakdancer_reader import BreakDancerReader
from pindel_reader import PindelReader
from cnvnator_reader import CNVnatorReader

tool_to_reader = {"BreakDancer": BreakDancerReader, "Pindel": PindelReader, "CNVnator": CNVnatorReader}


def convert_svtool_to_vcf(file_name, sample, out_vcf, toolname):
    vcf_template_reader = vcf.Reader(open(os.path.join(mydir, "resources/template.vcf"), "r"))
    vcf_template_reader.samples = [sample]

    vcf_fd = open(out_vcf, "w") if out_vcf is not None else sys.stdout
    vcf_writer = vcf.Writer(vcf_fd, vcf_template_reader)

    for tool_record in tool_to_reader[toolname](file_name):
        vcf_record = tool_record.to_vcf_record(sample)
        if vcf_record is None:
            continue
        vcf_writer.write_record(vcf_record)
    vcf_writer.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser("Convert SV tool output file to VCF",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--input", help="SV tool output file", required=False)
    parser.add_argument("--output", help="Output VCF to create", required=False)
    parser.add_argument("--tool", help="Tool name", required=False, default="BreakDancer",
                        choices=["BreakDancer", "CNVnator", "Pindel"])
    parser.add_argument("--sample", help="Sample name", required=True)

    args = parser.parse_args()
    convert_svtool_to_vcf(args.input, args.sample, args.output, args.tool)
