import argparse
import sys
import datetime
import os
from collections import OrderedDict
import json
import base64
import logging

import pybedtools
import pysam
import vcf

import fasta_utils

mydir = os.path.dirname(os.path.realpath(__file__))
vcf_template = os.path.join(mydir, "resources/template.vcf")


def convert_metasv_bed_to_vcf(bedfile=None, vcf_out=None, vcf_template_file=vcf_template, sample=None, reference=None,
                              pass_calls=True):
    func_logger = logging.getLogger("%s" % (convert_metasv_bed_to_vcf.__name__))

    vcf_template_reader = vcf.Reader(open(vcf_template_file, "r"))

    # The following are hacks to ensure sample name and contig names are put in the VCF header
    vcf_template_reader.samples = [sample]
    contigs = []
    if reference:
        contigs = fasta_utils.get_contigs(reference)
        contigs_order_dict = {contig.name: index for (index, contig) in enumerate(contigs)}
        vcf_template_reader.contigs = OrderedDict([(contig.name, (contig.name, contig.length)) for contig in contigs])
        vcf_template_reader.metadata["reference"] = reference

    vcf_template_reader.metadata["fileDate"] = str(datetime.date.today())
    vcf_template_reader.metadata["source"] = [" ".join(sys.argv)]

    vcf_writer = vcf.Writer(open(vcf_out, "w"), vcf_template_reader)

    vcf_records = []
    if bedfile:
        for interval in pybedtools.BedTool(bedfile):
            chrom = interval.chrom
            pos = interval.start
            end = interval.end
            genotype = "./." if len(interval.fields) < 11 else interval.fields[10]

            if genotype == "0/0":
                func_logger.info("Skipping homozygous reference %s" % str(interval))
                continue

            sub_names = interval.name.split(":")
            sub_lengths = map(lambda x: int(x.split(",")[2]), sub_names)

            sub_types = map(lambda x: x.split(",")[1], sub_names)
            sub_methods = [name.split(",")[3] for name in sub_names]
            svmethods = (";".join([name.split(",")[3] for name in sub_names])).split(";")
            try:
                info = json.loads(base64.b64decode(name.split(",")[0]))
            except TypeError:
                info = dict()
            if len(interval.fields) > 9:
                info.update(json.loads(base64.b64decode(interval.fields[9])))
            
            
            index_to_use = 0
            is_pass = False
            svlen = -1
            if "DEL" in sub_types:
                index_to_use = sub_types.index("DEL")
                svmethods_s = set(svmethods) - {"SC","AS"}
                is_pass = len(svmethods_s) > 1
            elif "INV" in sub_types:
                index_to_use = sub_types.index("INV")
                svmethods_s = set(svmethods) - {"SC"}
                is_pass = len(svmethods_s) > 1
                if "AS" in svmethods_s:
                    pos = int(interval.fields[6])
                    end = int(interval.fields[7])
                    svlen = int(interval.fields[8])
                    is_pass = (int(interval.fields[8]) != -1) and (svlen >= 100)
                else:
                	is_pass = len(svmethods_s) > 1
                
            elif "INS" in sub_types and "SC" in sub_methods:
                # TODO: I think it should be sub_types.index
                index_to_use = sub_methods.index("SC")
                pos = int(interval.fields[6])
                end = int(interval.fields[7])
                svlen = int(interval.fields[8])
            elif "ITX" in sub_types:
                index_to_use = sub_types.index("ITX")
                svmethods_s = set(svmethods) - {"SC","AS"}
                is_pass = len(svmethods_s) > 1
                end = info["END"]
            elif "CTX" in sub_types:
                index_to_use = sub_types.index("CTX")
                svmethods_s = set(svmethods) - {"SC","AS"}
                is_pass = len(svmethods_s) > 1
                end = info["END"]
                
            if pos < 1:
                func_logger.info("Variant with pos < 1 encountered. Skipping! %s" % str(interval))
                continue

            if svlen < 0: svlen = sub_lengths[index_to_use]
            if sub_types[index_to_use] == "DEL":
                svlen = -svlen

            sv_type = sub_types[index_to_use]
            if sv_type == "INS":
                if pass_calls and end != pos + 1:
                    continue
                end = pos
                is_pass = (int(interval.fields[8]) != -1) and (svlen == 0 or svlen >= 100)
            sv_id = "."
            ref = "."
            alt = [vcf.model._SV(sv_type)]
            qual = "."
            sv_filter = ["PASS" if is_pass else "LowQual"]
            info.update(
                {"END": end, "SVLEN": svlen, "SVTYPE": sv_type, "SVMETHOD": svmethods, "NUM_SVMETHODS": len(svmethods)})
            sv_format = "GT"
            sample_indexes = [0]
            vcf_record = vcf.model._Record(chrom, pos, sv_id, ref, alt, qual, sv_filter, info, sv_format, sample_indexes)
            vcf_record.samples = vcf_template_reader._parse_samples([genotype], "GT", vcf_record)
            vcf_records.append(vcf_record)

    if contigs:
        vcf_records.sort(key=lambda x: (contigs_order_dict[x.CHROM], x.POS))
    else:
        vcf_records.sort(key=lambda x: (x.CHROM, x.POS))

    for vcf_record in vcf_records:
        vcf_writer.write_record(vcf_record)
    vcf_writer.close()

    func_logger.info("Tabix compressing and indexing %s" % vcf_out)
    pysam.tabix_index(vcf_out, force=True, preset="vcf")


if __name__ == "__main__":
    FORMAT = '%(levelname)s %(asctime)-15s %(name)-20s %(message)s'
    logging.basicConfig(level=logging.INFO, format=FORMAT)
    logger = logging.getLogger(__name__)

    parser = argparse.ArgumentParser(description="Convert MetaSV final BED to VCF",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("--sample", help="Sample name", required=True)
    parser.add_argument("--bed", help="MetaSV final BED", required=True)
    parser.add_argument("--vcf", help="Final VCF to output", required=True)
    parser.add_argument("--vcf_template", help="VCF template", default=vcf_template)
    parser.add_argument("--reference", help="Reference FASTA", required=False)
    parser.add_argument("--pass_only", action="store_true", help="Output only PASS calls")

    args = parser.parse_args()

    convert_metasv_bed_to_vcf(bedfile=args.bed, vcf_out=args.vcf, vcf_template_file=args.vcf_template,
                              sample=args.sample,
                              reference=args.reference, pass_calls=args.pass_only)
