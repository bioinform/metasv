#!/net/kodiak/volumes/lake/shared/users/marghoob/my_env/bin/python

import sys
import os
import argparse
import subprocess
import shutil
import logging
import copy
from collections import defaultdict

import pysam
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
from generate_final_vcf import convert_metasv_bed_to_vcf


FORMAT = '%(levelname)s %(asctime)-15s %(name)-20s %(message)s'
logging.basicConfig(level=logging.INFO, format=FORMAT)
logger = logging.getLogger(__name__)

parser = argparse.ArgumentParser("Merge SVs from different tools",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("--sample", metavar="Sample", help="Sample name", required=True)
parser.add_argument("--pindel_vcf", nargs="+", metavar="pindel_vcf", help="VCF file or dir for Pindel VCFs",
                    required=False, default=[])
parser.add_argument("--pindel_native", nargs="+", metavar="File list", help="Pindel native files", required=False,
                    default=[])
parser.add_argument("--breakdancer_vcf", nargs="+", metavar="breakdancer_vcf",
                    help="VCF file or dir for BreakDancer VCFs", required=False, default=[])
parser.add_argument("--breakdancer_native", nargs="+", metavar="File list", help="BreakDancer native files",
                    required=False, default=[])
parser.add_argument("--breakseq_vcf", nargs="+", metavar="breakseq_vcf", help="VCF file or dir for BreakSeq VCFs",
                    required=False, default=[])
parser.add_argument("--cnvnator_vcf", nargs="+", metavar="cnvnator_vcf", help="VCF file or dir for CNVnator VCFs",
                    required=False, default=[])
parser.add_argument("--cnvnator_native", nargs="+", metavar="File list", help="CNVnator native files", required=False,
                    default=[])
parser.add_argument("--gatk_vcf", nargs="+", metavar="file", help="VCF file or dir for gatk VCFs", required=False,
                    default=[])
parser.add_argument("--reference", metavar="reference", help="Reference file", required=True)
parser.add_argument("--gaps", metavar="gaps", help="Gap bed file", required=False, default=None)
parser.add_argument("--filter_gaps", help="Filter out gaps", action="store_true", required=False)
parser.add_argument("--keep_standard_contigs", action="store_true", help="Keep only the major contigs + MT")
parser.add_argument("--wiggle", help="Wiggle for interval overlap", default=100, type=int, required=False)
parser.add_argument("--inswiggle", help="Wiggle for insertions, overides wiggle", default=100, type=int, required=False)
parser.add_argument("--minsvlen", help="Minimum length acceptable to be an SV", default=50, type=int, required=False)
parser.add_argument("--overlap_ratio", help="Reciprocal overlap ratio", default=0.5, type=float, required=False)
parser.add_argument("--workdir", help="Scratch directory for working", default="work", required=False)
parser.add_argument("--boost_ins", help="Use soft-clips for improving insertion detection", action="store_true")
parser.add_argument("--bam", help="BAM", type=file)
parser.add_argument("--chromosomes",
                    help="Chromosome list to process. If unspecified, then all chromosomes will be considered.",
                    nargs="+", default=[])
parser.add_argument("--num_threads", help="Number of threads to use", type=int, default=1)
parser.add_argument("--outdir", help="Output directory", required=True)
parser.add_argument("--spades", help="Path to SPAdes executable", required=False)
parser.add_argument("--age", help="Path to AGE executable", required=False)
parser.add_argument("--disable_assembly", action="store_true", help="Disable assembly")

args = parser.parse_args()
mydir = os.path.dirname(os.path.realpath(__file__))

if not os.path.isdir(args.outdir):
    logger.info("Creating directory %s" % (args.outdir))
    os.makedirs(args.outdir)

if not os.path.isdir(args.workdir):
    logger.info("Creating directory %s" % (args.workdir))
    os.makedirs(args.workdir)

bedtools_tmpdir = os.path.join(args.workdir, "bedtools")

# Initial checks
if (args.pindel_vcf is None) and (args.breakseq_vcf is None) and (args.breakdancer_vcf is None) and (
    args.cnvnator_vcf is None) and (args.pindel_native is None):
    logger.error("Nothing to do since no SV data specified")
    sys.exit(1)

if not os.path.isfile(args.reference + ".fai"):
    logger.error("Reference file %s is not indexed" % (args.reference))
    sys.exit(1)

# These are the included reference intervals
fasta_handle = pysam.Fastafile(args.reference)
contigs = get_contigs(args.reference)
include_intervals = sorted(
    [SVInterval(contig.name, 0, contig.length, contig.name, "include", length=contig.length) for contig in contigs])
logger.info(include_intervals)

# These are the contigs specified by the user to keep
contig_whitelist = set(args.chromosomes) if args.chromosomes else set([contig.name for contig in contigs])

if args.keep_standard_contigs:
    contig_whitelist &= set(
        [str(i) for i in xrange(1, 23)] + ["chr%d" % (i) for i in xrange(1, 23)] + ["X", "Y", "MT", "chrX", "chrY",
                                                                                    "chrM"])
logger.info("Only SVs on the following contigs will be reported: %s" % (sorted(list(contig_whitelist))))

vcf_name_list = [("CNVnator", args.cnvnator_vcf),
                 ("Pindel", args.pindel_vcf),
                 ("BreakDancer", args.breakdancer_vcf),
                 ("BreakSeq", args.breakseq_vcf),
                 ("HaplotypeCaller", args.gatk_vcf)]

native_name_list = [("CNVnator", args.cnvnator_native, CNVnatorReader),
                    ("Pindel", args.pindel_native, PindelReader),
                    ("BreakDancer", args.breakdancer_native, BreakDancerReader)]

tools = []
intervals = {}
sv_types = set()

# load gap files for filtering SVs that come from gaps in the reference
gap_intervals = []
if args.filter_gaps:
    if args.gaps is None:
        # Try to guess
        if "chr1" in contig_whitelist:
            gap_intervals = load_gap_intervals(os.path.join(mydir, "resources/hg19.gaps.bed"))
        elif "1" in contig_whitelist:
            gap_intervals = load_gap_intervals(os.path.join(mydir, "resources/b37.gaps.bed"))
        else:
            logger.error("Couldn't guess gaps file for reference. No gap filtering will be done")
    else:
        gap_intervals = sorted(load_gap_intervals(args.gaps))

# Handles native input
for toolname, nativename, svReader in native_name_list:
    # If no native file is given, ignore the tool
    if not nativename: continue

    tools.append(toolname)
    intervals[toolname] = defaultdict(list)

    for native_file in nativename:
        for record in svReader(native_file):
            interval = record.to_sv_interval()

            # Check length
            if interval.sv_len < args.minsvlen:
                continue

            # Set wiggle
            if interval.sv_type == "INS":
                interval.wiggle = max(args.inswiggle,args.wiggle)
            else:
                interval.wiggle = args.wiggle

            if not interval:
                # This is the case for SVs we want to skip
                continue
            if not interval_overlaps_interval_list(interval, gap_intervals) and interval.chrom in contig_whitelist:
                intervals[toolname][record.sv_type].append(interval)

    sv_types |= set(intervals[toolname].keys())

# Handles the VCF input cases, we will just deal with these cases
for toolname, vcfname in vcf_name_list:
    # If no VCF is given, ignore the tool
    if not vcfname: continue

    tools.append(toolname)
    intervals[toolname] = {}

    vcf_list = []
    for vcffile in vcfname:
        if os.path.isdir(vcffile):
            logger.info("Will load from per-chromosome VCFs from directory %s for tool %s" % (vcffile, toolname))
            vcf_list += [os.path.join(vcffile, "%s.vcf.gz" % (contig.name)) for contig in contigs if
                         (not contig_whitelist or contig.name in contig_whitelist)]
        else:
            vcf_list.append(vcffile)

    for vcffile in vcf_list:
        load_intervals(vcffile, intervals[toolname], gap_intervals, include_intervals, toolname, contig_whitelist,
                       minsvlen=args.minsvlen, wiggle=args.wiggle, inswiggle=args.inswiggle)
    sv_types |= set(intervals[toolname].keys())

logger.info("SV types are %s" % (str(sv_types)))
tool_merged_intervals = {}
final_intervals = []

bd_out = os.path.join(args.outdir, "breakdancer.vcf")
pindel_out = os.path.join(args.outdir, "pindel.vcf")
cnvnator_out = os.path.join(args.outdir, "cnvnator.vcf")
breakseq_out = os.path.join(args.outdir, "breakseq.vcf")

vcf_out_list = [("BreakDancer", bd_out),
                ("Pindel", pindel_out),
                ("CNVnator", cnvnator_out),
                ("BreakSeq", breakseq_out)]

# This will just output per-tool VCFs, no intra-tool merging is done yet
for toolname, tool_out in vcf_out_list:
    if tool_out is not None and toolname in intervals:
        logger.info("Outputting single tool VCF for %s" % (str(toolname)))
        vcf_template_reader = vcf.Reader(open(os.path.join(mydir, "resources/template.vcf"), "r"))
        vcf_template_reader.samples = [args.sample]

        intervals_tool = []
        tool_out_fd = open(tool_out, "w")
        vcf_writer = vcf.Writer(tool_out_fd, vcf_template_reader)
        chr_intervals_tool = {contig.name: [] for contig in contigs}
        for sv_type in sv_types:
            if sv_type in intervals[toolname]:
                intervals_tool.extend([copy.deepcopy(interval) for interval in intervals[toolname][sv_type]])
        for interval in intervals_tool:
            # Marghoob says that this is just to fill-in some metadata
            interval.do_validation(args.overlap_ratio)

            interval.fix_pos()
            chr_intervals_tool[interval.chrom].append(interval)

        for contig in contigs:
            chr_intervals_tool[contig.name].sort()
            for interval in chr_intervals_tool[contig.name]:
                vcf_record = interval.to_vcf_record(fasta_handle,args.sample)
                if vcf_record is not None:
                    vcf_writer.write_record(vcf_record)
        tool_out_fd.close()
        vcf_writer.close()
        logger.info("Indexing single tool VCF for %s" % (str(toolname)))
        pysam.tabix_index(tool_out, force=True, preset="vcf")

# Do merging here
for sv_type in sv_types:
    logger.info("Processing SVs of type %s" % (sv_type))
    tool_merged_intervals[sv_type] = []

    # Do the intra-tool merging
    for tool in tools:
        if sv_type not in intervals[tool]: continue
        logger.info("First level merging for %s for tool %s" % (sv_type, tool))
        tool_merged_intervals[sv_type] += sv_interval.merge_intervals(intervals[tool][sv_type])

    # Do the inter-tool merging
    merged_intervals = sv_interval.merge_intervals(tool_merged_intervals[sv_type])

    # Marghoob.... what does this do? :/ please add some comments
    intervals1 = []
    intervals2 = []
    for interval in tool_merged_intervals[sv_type]:
        if interval_overlaps_interval_list(interval, merged_intervals, args.overlap_ratio, args.overlap_ratio):
            intervals2.append(interval)
        else:
            intervals1.append(interval)

    final_intervals.extend(sv_interval.merge_intervals(intervals1) + sv_interval.merge_intervals(intervals2))

final_chr_intervals = {}
for contig in contigs: final_chr_intervals[contig.name] = []
for interval in final_intervals:
    interval.do_validation(args.overlap_ratio)  # [Marghoob..] I'm pretty sure it checks for 2 tools support here
    interval.fix_pos()
    final_chr_intervals[interval.chrom].append(interval)

# This is the merged VCF without assembly, ok for deletions at this point
vcf_template_reader = vcf.Reader(open(os.path.join(mydir, "resources/template.vcf"), "r"))
vcf_template_reader.samples = [args.sample]
out_vcf = os.path.join(args.outdir, "metasv.vcf")
vcf_fd = open(out_vcf, "w") if out_vcf is not None else sys.stdout
vcf_writer = vcf.Writer(vcf_fd, vcf_template_reader)

final_stats = {}

# This is the bedfile for the merged VCF, used for assembly
bed_intervals = []
merged_bed = os.path.join(args.workdir, "metasv.bed")

# output both VCF and BED here
for contig in contigs:
    final_chr_intervals[contig.name].sort()
    for interval in final_chr_intervals[contig.name]:
        vcf_record = interval.to_vcf_record(fasta_handle,args.sample)
        if vcf_record is not None:
            key = (interval.sv_type,
                   "PASS" if interval.is_validated else "LowQual",
                   "PRECISE" if interval.is_precise else "IMPRECISE",
                   tuple(sorted(list(interval.sources))))
            if key not in final_stats: final_stats[key] = 0
            final_stats[key] += 1
            vcf_writer.write_record(vcf_record)
        bed_interval = interval.to_bed_interval(args.sample)
        if bed_interval is not None:
            bed_intervals.append(bed_interval)

pybedtools.BedTool(bed_intervals).saveas(merged_bed)
vcf_fd.close()
vcf_writer.close()

for key in sorted(final_stats.keys()):
    logger.info(str(key) + ":" + str(final_stats[key]))

# Work on assembly after this point

if not args.disable_assembly:
    if args.spades is None:
        logger.error("Spades executable not specified")
        exit(1)

    if args.age is None:
        logger.error("AGE executable not specified")
        exit(1)

    spades_tmpdir = os.path.join(args.workdir, "spades")
    if not os.path.isdir(spades_tmpdir):
        logger.info("Creating directory %s" % (spades_tmpdir))
        os.makedirs(spades_tmpdir)

    age_tmpdir = os.path.join(args.workdir, "age")
    if not os.path.isdir(age_tmpdir):
        logger.info("Creating directory %s" % (age_tmpdir))
        os.makedirs(age_tmpdir)

    assembly_bed = merged_bed

    if args.boost_ins:
        logger.info("Generating intervals for insertions")
        assembly_bed = parallel_generate_sc_intervals([args.bam.name], list(contig_whitelist), merged_bed, args.workdir,
                                                      num_threads=args.num_threads)
        logger.info("Generated intervals for assembly in %s" % (assembly_bed))

    logger.info("Will run assembly now")

    assembled_fasta, ignored_bed = run_spades_parallel(bam=args.bam.name, spades=args.spades, bed=assembly_bed,
                                                       work=spades_tmpdir, pad=500, nthreads=args.num_threads,
                                                       chrs=list(contig_whitelist))
    breakpoints_bed = run_age_parallel(intervals_bed=assembly_bed, reference=args.reference, assembly=assembled_fasta,
                                       pad=500, age=args.age, chrs=list(contig_whitelist), nthreads=args.num_threads,
                                       min_contig_len=100, age_workdir=age_tmpdir)

    final_bed = os.path.join(args.workdir, "final.bed")
    if ignored_bed:
        pybedtools.BedTool(breakpoints_bed).cat(pybedtools.BedTool(ignored_bed), postmerge=False).sort().saveas(
            final_bed)
    else:
        pybedtools.BedTool(breakpoints_bed).saveas(final_bed)
else:
    final_bed = merged_bed

final_vcf = os.path.join(args.workdir, "final.vcf")
convert_metasv_bed_to_vcf(bedfile=final_bed, vcf_out=final_vcf,
                          vcf_template=os.path.join(mydir, "resources/template.vcf"), sample=args.sample)

pybedtools.cleanup(remove_all=True)
