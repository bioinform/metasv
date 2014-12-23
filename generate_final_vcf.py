#!/net/kodiak/volumes/lake/shared/users/marghoob/my_env/bin/python

import argparse
import pybedtools
import pysam
import vcf
import logging

FORMAT = '%(levelname)s %(asctime)-15s %(name)-20s %(message)s'
logging.basicConfig(level=logging.INFO, format=FORMAT)
logger = logging.getLogger(__name__)

parser = argparse.ArgumentParser("Convert MetaSV final BED to VCF", formatter_class = argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("--sample", help = "Sample name", required = True)
parser.add_argument("--bed", help = "MetaSV final BED", required = True)
parser.add_argument("--vcf", help = "Final VCF to output", required = True)
parser.add_argument("--vcf_template", help = "VCF template", required = True)

args = parser.parse_args()

vcf_template = vcf.Reader(open(args.vcf_template, "r"), compressed = True)
vcf_writer = vcf.Writer(open(args.vcf, "w"), vcf_template)

for interval in pybedtools.BedTool(args.bed):
  chrom = interval.chrom
  pos = interval.start
  end = interval.end

  sub_names = interval.name.split(":")
  sub_lengths = map(int, interval.fields[5].split(":"))

  sub_types = map(lambda x: x.split(",")[0], sub_names)
  sub_methods = [name.split(",")[2] for name in sub_names]
  svmethods = (";".join([name.split(",")[2] for name in sub_names])).split(";")

  index_to_use = 0
  should_ignore = False
  if "DEL" in sub_types:
    index_to_use = sub_types.index("DEL")
    svmethods_s = set(svmethods) - set(["SC"])
    if len(svmethods_s) == 1: continue
  elif "INV" in sub_types:
    index_to_use = sub_types.index("INV")
    svmethods_s = set(svmethods) - set(["SC"])
    if len(svmethods_s) == 1: continue
  elif "INS" in sub_types and "SC" in sub_methods:
    index_to_use = sub_methods.index("SC")
  #elif "INS" in sub_types:
  #  continue

  svlen = sub_lengths[index_to_use]
  if sub_types[index_to_use] == "DEL":
    svlen = -svlen

  sv_type = sub_types[index_to_use]
  if sv_type == "INS":
    if end != pos + 1: continue
    end = pos
  sv_id = "."
  ref = "."
  alt = ["<%s>" % (sv_type)]
  qual = "."
  sv_filter = "."
  info = {"END": end, "SVLEN": svlen, "SVTYPE": sv_type, "SVMETHOD": svmethods}
  sv_format = "GT"
  sample_indexes = [0]
  samples = [vcf.model._Call(None, args.sample, ["1/1"])]
  vcf_record = vcf.model._Record(chrom, pos, sv_id, ref, alt, qual, sv_filter, info, sv_format, sample_indexes, samples)

  vcf_writer.write_record(vcf_record)

vcf_writer.close()
