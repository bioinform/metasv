#!/net/kodiak/volumes/lake/shared/users/marghoob/my_env/bin/python

from __future__ import print_function
import pysam
import argparse
import os
import sys
import logging
import subprocess
import re
import hashlib
import uuid
import vcf
from multiprocessing import Process
import pybedtools
from pybedtools.parallel import parallel_apply
import random
from functools import partial
from itertools import compress


FORMAT = '%(levelname)s %(asctime)-15s %(message)s'
logging.basicConfig(level=logging.INFO, format=FORMAT)
logger = logging.getLogger(__name__)

parser = argparse.ArgumentParser("Annotate BED file with stuff", formatter_class = argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--in_bed", help = "BED file to annotate", required = True)
parser.add_argument("--rmask_bed", help = "Repeat masker BED")
parser.add_argument("--segdups_bed", help = "Segdups BED")
parser.add_argument("--truth_bed", help = "Truth BED file")
parser.add_argument("--out_bed", help = "Output BED file")
parser.add_argument("--mappability_bed", help = "Mappability BED file")
parser.add_argument("--conservation_bed", help = "Conservation score BED")
parser.add_argument("--dgv_bed", help = "DGV BED")
parser.add_argument("--reference", help = "Reference FASTA file")
parser.add_argument("--tmpdir", help = "Tmp directory")
parser.add_argument("--bam", help = "alignment BAM")
parser.add_argument("--pre_asm", help = "Pre-assembly interval bed")

args = parser.parse_args()

def add_bool(interval, flag = False, use_bool = True):
  fields = interval.fields
  fields.append(str(flag if use_bool else int(flag)))
  return pybedtools.create_interval_from_list(fields)

def remove_extra_nuc_fields(interval):
  fields = interval.fields[:-7]
  del fields[-2]
  return pybedtools.create_interval_from_list(fields)

def fix_pre_asm_overlap(interval):
  supports = interval.fields[-3].split(",")
  fields = interval.fields[:-6] + [interval.fields[-2]] + supports[:4] + [",".join(supports[4:])]
  return pybedtools.create_interval_from_list(fields)

def fix_truth_overlap(interval):
  interval_start = interval.start
  is_true = interval.fields[-5] != "."
  offset_start = (int(interval.fields[-4]) - interval_start) if is_true else -1
  offset_end = (interval.end - int(interval.fields[-4])) if is_true else -1
  fields = interval.fields[:-5] + map(str, [interval.fields[-4], offset_start, offset_end, interval.fields[-2], is_true])
  return pybedtools.create_interval_from_list(fields)

def annotate_bed(in_bed, other_bed, use_bool = False):
  true_bed = (in_bed + other_bed).each(partial(add_bool, flag = True, use_bool = use_bool)).saveas(os.path.join(args.tmpdir, "1.bed"))
  false_bed = (in_bed - other_bed).each(partial(add_bool, flag = False, use_bool = use_bool)).saveas(os.path.join(args.tmpdir, "2.bed"))

  return true_bed.cat(false_bed, postmerge = False, force_truncate = False).saveas(os.path.join(args.tmpdir, "cat.bed")).sort()

def add_weighted_score(in_bed, score_bed):
  out_bed = in_bed.intersect(score_bed, wao = True).saveas(os.path.join(args.tmpdir, "score.bed"))

  bed_array = []
  last_interval = pybedtools.Interval("", 0, 0)
  map_value = 0.0
  for interval in out_bed:
    if interval.chrom != last_interval.chrom or interval.start != last_interval.start or interval.end != last_interval.end:
      if last_interval.chrom:
        bed_array.append(tuple(last_interval.fields[:-5]) + (str(map_value),))
      map_value = 0.0
      last_interval = interval
    if float(interval.fields[-1]) > 0:
      map_value += float(interval.fields[-1]) * float(interval.fields[-2]) / float(interval.length)

  if last_interval.chrom:
    bed_array.append(tuple(last_interval.fields[:-5]) + (str(map_value),))

  return pybedtools.BedTool(bed_array)

def annotate_coverage(interval, bam = None):
  fields = interval.fields

  coverages = []
  very_good_coverages = []
  for start in xrange(interval.start, interval.end, 50):
    end = start + 50
    this_coverage = 0
    this_very_good_coverage = 0
    for aln in bam.fetch(interval.chrom, start, end):
      if aln.is_duplicate: continue
      #if aln.mapq < 2 or aln.is_duplicate: continue
      this_coverage += 1

      if aln.mapq >= 20 and aln.cigar is not None and int(aln.opt("NM")) < 2:
        total_softclip = sum([length for (op, length) in aln.cigar if op == 4])
        if total_softclip <= 10: this_very_good_coverage += 1

    coverages.append(this_coverage)
    very_good_coverages.append(this_very_good_coverage)

  min_cov = min(coverages)
  max_cov = max(coverages)
  average_coverage = float(sum(coverages)) / len(coverages)
  average_very_good_coverage = float(sum(very_good_coverages)) / len(very_good_coverages)

  possible = 1 if (average_coverage > 0 and min_cov / average_coverage < 0.6) else 0
  fields.append(";".join([str(i) for i in coverages]))
  fields.append(str(average_coverage))
  fields.append(str(min_cov))
  fields.append(str(max_cov))
  fields.append(str(min_cov / average_coverage if average_coverage > 0 else 1))
  fields.append(str(max_cov / average_coverage if average_coverage > 0 else 1))
  fields.append(";".join([str(i) for i in very_good_coverages]))
  fields.append(str(average_very_good_coverage))
  fields.append(str(min(very_good_coverages)))
  fields.append(str(max(very_good_coverages)))
  fields.append(str(min(very_good_coverages) / average_very_good_coverage if (average_very_good_coverage > 0) else 1))
  fields.append(str(max(very_good_coverages) / average_very_good_coverage if (average_very_good_coverage > 0) else 1))

  return pybedtools.create_interval_from_list(fields)

def add_coverage_information(in_bed, bam):
  return in_bed.each(partial(annotate_coverage, bam = bam))

pybedtools.set_tempdir(args.tmpdir)

in_bed = pybedtools.BedTool(args.in_bed)

out_bed = in_bed
logger.info("Initial feature count %d" % (out_bed.count()))

if not os.path.isdir(args.tmpdir):
  os.makedirs(args.tmpdir)

bed_fields = ["#CHROM", "START", "END"]
bed_fields += ["NUM_CONTIGS_USED", "TOTAL_CONTIGS_COUNT"]
bed_fields += map(lambda x: "VERYGOOD_%s" % (x), ["NUM_ASMS", "NUM_UNIQUE_ASMS", "HAS_ASM", "IS_CONSISTENT", "INSERTION_LENGTH", "ASM_START", "ASM_END", "EXCISION_REF", "EXCISION_ASM"])
bed_fields += map(lambda x: "GOOD_%s" % (x), ["NUM_ASMS", "NUM_UNIQUE_ASMS", "HAS_ASM", "IS_CONSISTENT", "INSERTION_LENGTH", "ASM_START", "ASM_END", "EXCISION_REF", "EXCISION_ASM"])
bed_fields += ["PAIRING_COUNT", "PAIRING_START", "PAIRING_END", "PAIRING_EXISTS_AND_CONSISTENT"]

if args.rmask_bed:
  out_bed = annotate_bed(out_bed, pybedtools.BedTool(args.rmask_bed))
  logger.info("Feature count after rmask %d" % (out_bed.count()))
  bed_fields += ["OVERLAPS_RMASK"]

if args.segdups_bed:
  out_bed = annotate_bed(out_bed, pybedtools.BedTool(args.segdups_bed))
  logger.info("Feature count after segdups %d" % (out_bed.count()))
  bed_fields += ["OVERLAPS_SEGDUPS"]

if args.mappability_bed:
  out_bed = add_weighted_score(out_bed, pybedtools.BedTool(args.mappability_bed))
  logger.info("Feature count after mappability %d" % (out_bed.count()))
  bed_fields += ["MAPPABILITY_SCORE"]

if args.conservation_bed:
  out_bed = add_weighted_score(out_bed, pybedtools.BedTool(args.conservation_bed))
  logger.info("Feature count after conservation %d" % (out_bed.count()))
  bed_fields += ["CONSERVATION_SCORE"]

if args.dgv_bed:
  out_bed = annotate_bed(out_bed, pybedtools.BedTool(args.dgv_bed))
  logger.info("Feature count after DGV %d" % (out_bed.count()))

if args.reference:
  out_bed = out_bed.nucleotide_content(C = True, fi = args.reference).each(remove_extra_nuc_fields).saveas(os.path.join(args.tmpdir, "gc.bed"))
  logger.info("Feature count after GC content %d" % (out_bed.count()))
  bed_fields += ["GC_CONTENT"]

if args.bam:
  out_bed = add_coverage_information(out_bed, pysam.Samfile(args.bam, "rb")).saveas(os.path.join(args.tmpdir, "bam.bed"))
  logger.info("Feature count after BAM %d" % (out_bed.count()))
  bed_fields += ["COVS", "AVG_COV", "MIN_COV", "MAX_COV", "MIN_COV_OVER_AVG_COV", "MAX_COV_OVER_AVG_COV"]
  bed_fields += ["VERY_GOOD_COVS", "AVG_VERY_GOOD_COV", "MIN_VERY_GOOD_COV", "MAX_VERY_GOOD_COV", "MIN_VERY_GOOD_COV_OVER_AVG_VERY_GOOD_COV", "MAX_VERY_GOOD_COV_OVER_AVG_VERY_GOOD_COV"]

if args.pre_asm:
  out_bed = out_bed.intersect(pybedtools.BedTool(args.pre_asm), r = True, f = 1.0, wo = True).each(fix_pre_asm_overlap).saveas(os.path.join(args.tmpdir, "pre_asm.bed"))
  logger.info("Feature count after pre-asm %d" % (out_bed.count()))
  bed_fields += ["READS_SUPPORTING_INTERVAL", "PLUS_SUPPORT", "MINUS_SUPPORT", "LOCATION_SPAN", "UNIQ_LOCS", "LOCS"]

if args.truth_bed:
  #out_bed = annotate_bed(out_bed.saveas(os.path.join(args.tmpdir, "out.bed")), pybedtools.BedTool(args.truth_bed), use_bool = True)
  out_bed = out_bed.intersect(pybedtools.BedTool(args.truth_bed), wao = True).each(fix_truth_overlap).saveas(os.path.join(args.tmpdir, "truth.bed"))
  logger.info("Feature count after truth %d" % (out_bed.count()))
  bed_fields += ["TRUE_POSITION", "TRUE_LENGTH", "OFFSET_INTERVAL_START", "OFFSET_INTERVAL_END", "OVERLAPS_TRUTH"]

out_fd = open(args.out_bed, "w")

with open(args.out_bed, "w") as out_fd:
  print("\t".join(bed_fields), file = out_fd)
  for interval in out_bed:
    out_fd.write(str(interval))
