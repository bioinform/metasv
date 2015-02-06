#!/usr/bin/python

import argparse
import sys
import datetime
import os
from collections import OrderedDict
import json
import base64

import pybedtools
import pysam
import vcf

import fasta_utils


mydir = os.path.dirname(os.path.realpath(__file__))
vcf_template = os.path.join(mydir, "resources/template.vcf")


def convert_metasv_bed_to_vcf(bedfile=None, vcf_out=None, vcf_template_file=vcf_template, sample=None, reference=None, pass_calls=True):
    vcf_template_reader = vcf.Reader(open(vcf_template_file, "r"))

    # The following are hacks to ensure sample name and contig names are put in the VCF header
    vcf_template_reader.samples = [sample]
    if reference:
        contigs = fasta_utils.get_contigs(reference)
        vcf_template_reader.contigs = OrderedDict([(contig.name, (contig.name, contig.length)) for contig in contigs])
        vcf_template_reader.metadata["reference"] = reference
        vcf_template_reader.metadata["fileDate"] = str(datetime.date.today())
        vcf_template_reader.metadata["source"] = [" ".join(sys.argv)]

    vcf_writer = vcf.Writer(open(vcf_out, "w"), vcf_template_reader)

    for interval in pybedtools.BedTool(bedfile):
        chrom = interval.chrom
        pos = interval.start
        end = interval.end

        sub_names = interval.name.split(":")
        sub_lengths = map(lambda x: int(x.split(",")[1]), sub_names)

        sub_types = map(lambda x: x.split(",")[0], sub_names)
        sub_methods = [name.split(",")[2] for name in sub_names]
        svmethods = (";".join([name.split(",")[2] for name in sub_names])).split(";")
        try:
            info = json.loads(base64.b64decode(name.split(",")[3]))
        except TypeError:
            info = dict()
        if len(interval.fields) > 9:
            info.update(json.loads(base64.b64decode(interval.fields[9])))

        index_to_use = 0
        should_ignore = False  # Marghoob this is not used?
        if "DEL" in sub_types:
            index_to_use = sub_types.index("DEL")
            svmethods_s = set(svmethods) - {"SC"}
            if pass_calls and len(svmethods_s) == 1:
                continue
        elif "INV" in sub_types:
            index_to_use = sub_types.index("INV")
            svmethods_s = set(svmethods) - {"SC"}
            if pass_calls and len(svmethods_s) == 1:
                continue
        elif "INS" in sub_types and "SC" in sub_methods:
            index_to_use = sub_methods.index("SC")

        svlen = sub_lengths[index_to_use]
        if sub_types[index_to_use] == "DEL":
            svlen = -svlen

        sv_type = sub_types[index_to_use]
        if sv_type == "INS":
            if pass_calls and end != pos + 1:
                continue
            end = pos
        sv_id = "."
        ref = "."
        alt = ["<%s>" % sv_type]
        qual = "."
        sv_filter = "."
        info.update({"END": end, "SVLEN": svlen, "SVTYPE": sv_type, "SVMETHOD": svmethods, "NUM_SVMETHODS": len(svmethods)})
        sv_format = "GT"
        sample_indexes = [0]
        samples = [vcf.model._Call(None, sample, ["1/1"])]
        vcf_record = vcf.model._Record(chrom, pos, sv_id, ref, alt, qual, sv_filter, info, sv_format, sample_indexes,
                                       samples)

        vcf_writer.write_record(vcf_record)

    vcf_writer.close()
    pysam.tabix_index(vcf_out, force=True, preset="vcf")


if __name__ == "__main__":
    parser = argparse.ArgumentParser("Convert MetaSV final BED to VCF",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("--sample", help="Sample name", required=True)
    parser.add_argument("--bed", help="MetaSV final BED", required=True)
    parser.add_argument("--vcf", help="Final VCF to output", required=True)
    parser.add_argument("--vcf_template", help="VCF template", required=True)
    parser.add_argument("--reference", help="Reference FASTA", required=False)

    args = parser.parse_args()

    convert_metasv_bed_to_vcf(bedfile=args.bed, vcf_out=args.vcf, vcf_template_file=args.vcf_template, sample=args.sample,
                              reference=args.reference)
