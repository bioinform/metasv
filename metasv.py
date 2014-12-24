#!/net/kodiak/volumes/lake/shared/users/marghoob/my_env/bin/python

import sys
import os
import argparse
import subprocess
import pysam
import shutil
import logging
import copy
from fasta_utils import *
from vcf_utils import *
import sv_interval
from pindel_reader import PindelReader
from breakdancer_reader import BreakDancerReader
from cnvnator_reader import CNVnatorReader
import pybedtools
from generate_sv_intervals import parallel_generate_sc_intervals
from run_spades import run_spades_parallel
from run_age import run_age_parallel

FORMAT = '%(levelname)s %(asctime)-15s %(name)-20s %(message)s'
logging.basicConfig(level=logging.INFO, format=FORMAT)
logger = logging.getLogger(__name__)

parser = argparse.ArgumentParser("Merge SVs from different tools", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("--sample", metavar="Sample", help="Sample name", required=True)
parser.add_argument("--pindel_vcf", nargs="+", metavar="pindel_vcf", help="VCF file or dir for Pindel VCFs", required=False, default=[])
parser.add_argument("--pindel_native", nargs="+", metavar="File list", help="Pindel native files", required=False, default=[])
parser.add_argument("--breakdancer_vcf", nargs="+", metavar="breakdancer_vcf", help="VCF file or dir for BreakDancer VCFs", required=False, default=[])
parser.add_argument("--breakdancer_native", nargs="+", metavar="File list", help="BreakDancer native files", required=False, default=[])
parser.add_argument("--breakseq_vcf", nargs="+", metavar="breakseq_vcf", help="VCF file or dir for BreakSeq VCFs", required=False, default=[])
parser.add_argument("--cnvnator_vcf", nargs="+", metavar="cnvnator_vcf", help="VCF file or dir for CNVnator VCFs", required=False, default=[])
parser.add_argument("--cnvnator_native", nargs="+", metavar="File list", help="CNVnator native files", required=False, default=[])
parser.add_argument("--gatk_vcf", nargs="+", metavar="file", help="VCF file or dir for gatk VCFs", required=False, default=[])
parser.add_argument("--reference", metavar="reference", help="Reference file", required=True)
parser.add_argument("--gaps", metavar="gaps", help="Gap bed file", required=False, default=None)
parser.add_argument("--filter_gaps", help="Filter out gaps", action="store_true", required=False)
parser.add_argument("--keep_standard_contigs", action="store_true", help="Keep only the major contigs + MT")
parser.add_argument("--wiggle", help="Wiggle for interval overlap", default=100, type=int, required=False)
parser.add_argument("--overlap_ratio", help="Reciprocal overlap ratio", default=0.5, type=float, required=False)
parser.add_argument("--workdir", help = "Scratch directory for working", default = "work", required = False)
parser.add_argument("--boost_ins", help = "Use soft-clips for improving insertion detection", action = "store_true")
parser.add_argument("--bam", help = "BAMs", type = file)
parser.add_argument("--chromosomes", help = "Chromosome list to process. If unspecified, then all chromosomes will be considered.", nargs = "+", default = [])
parser.add_argument("--num_threads", help = "Number of threads to use", type = int, default = 1)
parser.add_argument("--outdir", help = "Output directory", required = True)
parser.add_argument("--spades", help = "Path to SPAdes executable", required = True)
parser.add_argument("--age", help = "Path to AGE executable", required = True)
parser.add_argument("--disable_assembly", action = "store_true", help = "Disable assembly")

args = parser.parse_args()
mydir = os.path.dirname(os.path.realpath(__file__))

if not os.path.isdir(args.outdir):
  logger.info("Creating directory %s" % (args.outdir))
  os.makedirs(args.outdir)

if not os.path.isdir(args.workdir):
  logger.info("Creating directory %s" % (args.workdir))
  os.makedirs(args.workdir)

bedtools_tmpdir = os.path.join(args.workdir, "bedtools")

spades_tmpdir = os.path.join(args.workdir, "spades")
if not os.path.isdir(spades_tmpdir):
  logger.info("Creating directory %s" % (spades_tmpdir))
  os.makedirs(spades_tmpdir)

age_tmpdir = os.path.join(args.workdir, "age")
if not os.path.isdir(age_tmpdir):
  logger.info("Creating directory %s" % (age_tmpdir))
  os.makedirs(age_tmpdir)

# Initial checks
if (args.pindel_vcf is None) and (args.breakseq_vcf is None) and (args.breakdancer_vcf is None) and (args.cnvnator_vcf is None) and (args.pindel_native is None):
  logger.error("Nothing to do since no SV data specified")
  sys.exit(1)

if not os.path.isfile(args.reference + ".fai"):
  logger.error("Reference file %s is not indexed" % (args.reference))
  sys.exit(1)

fasta_handle = pysam.Fastafile(args.reference)
contigs = get_contigs(args.reference)

contig_whitelist = set(args.chromosomes) if args.chromosomes else set([contig.name for contig in contigs])
if args.keep_standard_contigs:
  contig_whitelist &= set([str(i) for i in xrange(1, 23)] + ["chr%d" % (i) for i in xrange(1, 23)] + ["X", "Y", "MT", "chrX", "chrY", "chrM"])
logger.info("Only SVs on the following contigs will be reported: %s" % (sorted(list(contig_whitelist))))

vcf_name_list = [("CNVnator", args.cnvnator_vcf), ("Pindel", args.pindel_vcf), ("BreakDancer", args.breakdancer_vcf), ("BreakSeq", args.breakseq_vcf), ("HaplotypeCaller", args.gatk_vcf)]

tools = []
intervals = {}
sv_types = set()

gap_intervals = []
if args.filter_gaps:
  if args.gaps is None:
    # Try to guess
    if "chr1" in contig_whitelist: gap_intervals = load_gap_intervals(os.path.join(mydir, "resources/hg19.gaps.bed"))
    elif "1" in contig_whitelist: gap_intervals = load_gap_intervals(os.path.join(mydir, "resources/b37.gaps.bed"))
    else: logger.error("Couldn't guess gaps file for reference. No gap filtering will be done")
  else:
    gap_intervals = load_gap_intervals(args.gaps)

pindel_lis = []
if args.pindel_native is not None:
  for pindel_native_file in args.pindel_native:
    for pindel_record in PindelReader(pindel_native_file, fasta_handle):
      if pindel_record.sv_type == "LI":
        interval = pindel_record.to_sv_interval()
        if not interval_overlaps_interval_list(interval, gap_intervals) and interval.chrom in contig_whitelist:
          pindel_lis.append(pindel_record.to_sv_interval())

if args.breakdancer_native is not None:
  for breakdancer_native_file in args.breakdancer_native:
    for breakdancer_record in BreakDancerReader(breakdancer_native_file):
      print repr(breakdancer_record.to_sv_interval())

if args.cnvnator_native is not None:
  for cnvnator_native_file in args.cnvnator_native:
    for cnvnator_record in CNVnatorReader(cnvnator_native_file):
      print cnvnator_record

for toolname, vcfname in vcf_name_list:
  if not vcfname: continue
  tools.append(toolname)
  intervals[toolname] = {}

  if toolname == "Pindel" and pindel_lis:
    intervals[toolname]["INS"] = pindel_lis
    sv_types |= set(["INS"])

  vcf_list = []
  for vcf in vcfname:
    if os.path.isdir(vcf):
      logger.info("Will load from per-chromosome VCFs from directory %s for tool %s" % (vcf, toolname))
      vcf_list += [os.path.join(vcf, "%s.vcf.gz" % (contig.name)) for contig in contigs if (not contig_whitelist or contig.name in contig_whitelist)]
    else:
      vcf_list.append(vcf)

  for vcf in vcf_list:
    load_intervals(vcf, intervals[toolname], gap_intervals, toolname, contig_whitelist, toolname == "HaplotypeCaller")
  sv_types |= set(intervals[toolname].keys())

logger.info("SV types are %s" % (str(sv_types)))
tool_merged_intervals = {}
final_intervals = []

bd_out = os.path.join(args.outdir, "breakdancer.vcf")
pindel_out = os.path.join(args.outdir, "pindel.vcf")
cnvnator_out = os.path.join(args.outdir, "cnvnator.vcf")
breakseq_out = os.path.join(args.outdir, "breakseq.vcf")
for toolname, tool_out in [("BreakDancer", bd_out), ("Pindel", pindel_out), ("CNVnator", cnvnator_out), ("BreakSeq", breakseq_out)]:
  if tool_out is not None and toolname in intervals:
    intervals_tool = []
    tool_out_fd = open(tool_out, "w")
    chr_intervals_tool = {contig.name: [] for contig in contigs}
    for sv_type in sv_types:
      if sv_type in intervals[toolname]:
        intervals_tool.extend([copy.deepcopy(interval) for interval in intervals[toolname][sv_type]])
    for interval in intervals_tool:
      interval.do_validation()
      interval.fix_pos()
      chr_intervals_tool[interval.chrom].append(interval)
    print_vcf_header(tool_out_fd, args.reference, contigs, args.sample)
    for contig in contigs:
      chr_intervals_tool[contig.name].sort()
      for interval in chr_intervals_tool[contig.name]:
        vcf_record = interval.to_vcf_record(fasta_handle)
        if vcf_record is not None:
          tool_out_fd.write("%s\n" % (vcf_record))
    tool_out_fd.close()
    pysam.tabix_index(tool_out, force=True, preset="vcf")

for sv_type in sv_types:
  logger.info("Processing SVs of type %s" % (sv_type))
  tool_merged_intervals[sv_type] = []
  for tool in tools:
    if sv_type not in intervals[tool]: continue
    logger.info("First level merging for %s for tool %s" % (sv_type, tool))
    tool_merged_intervals[sv_type] += sv_interval.merge_intervals(intervals[tool][sv_type])
  merged_intervals = sv_interval.merge_intervals(tool_merged_intervals[sv_type])
  intervals1 = []
  intervals2 = []
  for interval in tool_merged_intervals[sv_type]:
    if interval_overlaps_interval_list(interval, merged_intervals, 0.5, 0.5):
      intervals2.append(interval)
    else:
      intervals1.append(interval)

  final_intervals.extend(sv_interval.merge_intervals(intervals1) + sv_interval.merge_intervals(intervals2))

final_chr_intervals = {}
for contig in contigs: final_chr_intervals[contig.name] = []
for interval in final_intervals:
  interval.do_validation()
  interval.fix_pos()
  final_chr_intervals[interval.chrom].append(interval)

out_vcf = os.path.join(args.outdir, "metasv.vcf")
outfd = open(out_vcf, "w")
print_vcf_header(outfd, args.reference, contigs, args.sample)
final_stats = {}

bed_intervals = []
merged_bed = os.path.join(args.workdir, "metasv.bed")
for contig in contigs:
  final_chr_intervals[contig.name].sort()
  for interval in final_chr_intervals[contig.name]:
    vcf_record = interval.to_vcf_record(fasta_handle)
    if vcf_record is not None:
      key = (interval.sv_type, "PASS" if interval.is_validated else "LowQual", "PRECISE" if interval.is_precise else "IMPRECISE", tuple(sorted(list(interval.sources))))
      if key not in final_stats: final_stats[key] = 0
      final_stats[key] += 1
      outfd.write("%s\n" % (vcf_record))
    bed_interval = interval.to_bed_interval(args.sample)
    if bed_interval is not None:
      bed_intervals.append(bed_interval)

pybedtools.BedTool(bed_intervals).saveas(merged_bed)
outfd.close()
#pysam.tabix_index(out_vcf, force=True, preset="vcf")

for key in sorted(final_stats.keys()):
  logger.info(str(key) + ":" + str(final_stats[key]))

if not args.disable_assembly:
  assembly_bed = merged_bed
  #pybedtools.BedTool(bed_intervals).saveas(merged_bed)

  if args.boost_ins:
    logger.info("Generating intervals for insertions")
    assembly_bed = parallel_generate_sc_intervals([args.bam.name], list(contig_whitelist), merged_bed, args.workdir, num_threads = args.num_threads)
    logger.info("Generated intervals for assembly in %s" % (assembly_bed))

  logger.info("Will run assembly now")

  assembled_fasta, ignored_bed = run_spades_parallel(bam = args.bam.name, spades = args.spades, bed = assembly_bed, work = spades_tmpdir, pad = 500, nthreads = args.num_threads, chrs = list(contig_whitelist))
  breakpoints_bed = run_age_parallel(intervals_bed = assembly_bed, reference = args.reference, assembly = assembled_fasta, pad = 500, age = args.age, chrs = list(contig_whitelist), nthreads = args.num_threads, min_contig_len = 100, age_workdir = age_tmpdir)

  final_bed = os.path.join(args.workdir, "final.bed")
  if ignored_bed:
    pybedtools.BedTool(breakpoints_bed).cat(pybedtools.BedTool(ignored_bed), postmerge = False).sort().saveas(final_bed)
  else:
    pybedtools.BedTool(breakpoints_bed).saveas(final_bed)
else:
  final_bed = merged_bed

final_vcf = os.path.join(args.workdir, "final.vcf")
convert_metasv_bed_to_vcf(bedfile = final_bed, vcf_out = final_vcf, vcf_template = os.path.join(mydir, "resources/template.vcf"), sample = args.sample)

pybedtools.cleanup(remove_all = True)
