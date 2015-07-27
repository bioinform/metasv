#!/usr/bin/env python

import sys
import argparse
from metasv.main import run_metasv
from metasv.defaults import *
from metasv._version import __version__

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Merge SVs from multiple tools for accurate SV calling",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    input_parser = parser.add_argument_group("Input data options")
    input_parser.add_argument("--sample", metavar="Sample", help="Sample name", required=True)
    input_parser.add_argument("--pindel_vcf", nargs="+", metavar="pindel_vcf", help="VCF file or dir for Pindel VCFs",
                              required=False, default=[])
    input_parser.add_argument("--pindel_native", nargs="+", metavar="File list", help="Pindel native files",
                              required=False,
                              default=[])
    input_parser.add_argument("--breakdancer_vcf", nargs="+", metavar="breakdancer_vcf",
                              help="VCF file or dir for BreakDancer VCFs", required=False, default=[])
    input_parser.add_argument("--breakdancer_native", nargs="+", metavar="File list", help="BreakDancer native files",
                              required=False, default=[])
    input_parser.add_argument("--breakseq_vcf", nargs="+", metavar="breakseq_vcf",
                              help="VCF file or dir for BreakSeq VCFs",
                              required=False, default=[])
    input_parser.add_argument("--breakseq_native", nargs="+", metavar="breakseq_native",
                              help="BreakSeq native GFF files",
                              required=False, default=[])
    input_parser.add_argument("--cnvnator_vcf", nargs="+", metavar="cnvnator_vcf",
                              help="VCF file or dir for CNVnator VCFs",
                              required=False, default=[])
    input_parser.add_argument("--cnvnator_native", nargs="+", metavar="File list", help="CNVnator native files",
                              required=False,
                              default=[])
    input_parser.add_argument("--gatk_vcf", nargs="+", metavar="file", help="VCF file or dir for gatk VCFs",
                              required=False,
                              default=[])

    reference_parser = parser.add_argument_group("Reference options")
    reference_parser.add_argument("--reference", metavar="reference", help="Reference file", required=True)
    reference_parser.add_argument("--chromosomes",
                                  help="Chromosome list to process. If unspecified, then all chromosomes will be considered.",
                                  nargs="+", default=[])
    reference_parser.add_argument("--gaps", metavar="gaps", help="Gap bed file", required=False, default=None)
    reference_parser.add_argument("--filter_gaps", help="Filter out gaps", action="store_true", required=False)
    reference_parser.add_argument("--keep_standard_contigs", action="store_true",
                                  help="Keep only the major contigs + MT")

    bam_parser = parser.add_argument_group("Input BAM options")
    bam_parser.add_argument("--bam", help="BAM", type=file)
    bam_parser.add_argument("--isize_mean", type=float, default=ISIZE_MEAN, help="Insert size mean")
    bam_parser.add_argument("--isize_sd", type=float, default=ISIZE_SD, help="Insert size standard deviation")

    merging_parser = parser.add_argument_group("Tool output merging options")
    merging_parser.add_argument("--wiggle", help="Wiggle for interval overlap", default=WIGGLE, type=int,
                                required=False)
    merging_parser.add_argument("--inswiggle", help="Wiggle for insertions, overides wiggle", default=INS_WIGGLE,
                                type=int,
                                required=False)
    merging_parser.add_argument("--minsvlen", help="Minimum length acceptable to be an SV", default=MIN_SV_LENGTH,
                                type=int,
                                required=False)
    merging_parser.add_argument("--maxsvlen", help="Maximum length SV to report", default=MAX_SV_LENGTH,
                                type=int,
                                required=False)
    merging_parser.add_argument("--overlap_ratio", help="Reciprocal overlap ratio", default=OVERLAP_RATIO, type=float,
                                required=False)

    insertion_parser = parser.add_argument_group("Insertion detection options")
    insertion_parser.add_argument("--boost_ins", help="Use soft-clips for improving insertion detection",
                                  action="store_true")
    insertion_parser.add_argument("--min_avg_base_qual", help="Minimum average base quality",
                                  default=SC_MIN_AVG_BASE_QUAL, type=int)
    insertion_parser.add_argument("--min_mapq", help="Minimum MAPQ", default=SC_MIN_MAPQ, type=int)
    insertion_parser.add_argument("--min_soft_clip", help="Minimum soft-clip", default=SC_MIN_SOFT_CLIP, type=int)
    insertion_parser.add_argument("--max_soft_clip", help="Maximum soft-clip", default=SC_MAX_SOFT_CLIP, type=int)
    insertion_parser.add_argument("--max_nm", help="Maximum number of edits", default=SC_MAX_NM, type=int)
    insertion_parser.add_argument("--min_matches", help="Mininum number of matches", default=SC_MIN_MATCHES, type=int)
    insertion_parser.add_argument("--min_ins_support",
                                  help="Minimum read support for calling insertions using soft-clips",
                                  type=int, default=MIN_SUPPORT)
    insertion_parser.add_argument("--min_ins_support_frac",
                                  help="Minimum fraction of reads supporting insertion using soft-clips", type=float,
                                  default=MIN_SUPPORT_FRAC)
    insertion_parser.add_argument("--max_ins_intervals", help="Maximum number of insertion intervals to generate",
                                  type=int,
                                  default=MAX_INTERVALS)

    as_parser = parser.add_argument_group("Assembly options")
    as_parser.add_argument("--spades", help="Path to SPAdes executable", required=False)
    as_parser.add_argument("--disable_assembly", action="store_true", help="Disable assembly")
    as_parser.add_argument("--svs_to_assemble", nargs="+", help="SVs to assemble", default=SVS_ASSEMBLY_SUPPORTED, choices=SVS_ASSEMBLY_SUPPORTED)
    as_parser.add_argument("--extraction_max_read_pairs", type=int, default=EXTRACTION_MAX_READ_PAIRS,
                           help="Maximum number of pairs to extract for assembly")
    as_parser.add_argument("--spades_max_interval_size", type=int, default=SPADES_MAX_INTERVAL_SIZE, help="Maximum SV length for assembly")
    as_parser.add_argument("--stop_spades_on_fail", action="store_true", help="Abort on SPAdes failure")
    as_parser.add_argument("--age", help="Path to AGE executable", required=False)

    gt_parser = parser.add_argument_group("Genotyping options")
    gt_parser.add_argument("--gt_window", type=int, default=GT_WINDOW, help="Window for genotyping")
    gt_parser.add_argument("--gt_normal_frac", type=float, default=GT_NORMAL_FRAC,
                           help="Min. fraction of reads supporting reference for genotyping")

    out_parser = parser.add_argument_group("Output options")
    out_parser.add_argument("--svs_to_report", nargs="+", help="SV types to report", default=SVS_SUPPORTED, choices=SVS_SUPPORTED)
    out_parser.add_argument("--enable_per_tool_output", action="store_true",
                            help="Enable output of merged SVs for individual tools")

    work_parser = parser.add_argument_group("Running environment options")
    work_parser.add_argument("--workdir", help="Scratch directory for working", default="work", required=False)
    work_parser.add_argument("--num_threads", help="Number of threads to use", type=int, default=1)
    work_parser.add_argument("--outdir", help="Output directory", required=True)

    other_parser = parser.add_argument_group("Other options")
    other_parser.add_argument('--version', action='version', version='%(prog)s ' + __version__)

    args = parser.parse_args()

    sys.exit(run_metasv(args.sample, args.reference, pindel_vcf=args.pindel_vcf, pindel_native=args.pindel_native,
                        breakdancer_vcf=args.breakdancer_vcf, breakdancer_native=args.breakdancer_native,
                        breakseq_vcf=args.breakseq_vcf, breakseq_native=args.breakseq_native,
                        cnvnator_vcf=args.cnvnator_vcf,
                        cnvnator_native=args.cnvnator_native, gatk_vcf=args.gatk_vcf,
                        gaps=args.gaps, filter_gaps=args.filter_gaps, keep_standard_contigs=args.keep_standard_contigs,
                        wiggle=args.wiggle, overlap_ratio=args.overlap_ratio,
                        workdir=args.workdir, outdir=args.outdir, boost_ins=args.boost_ins, bam=args.bam,
                        chromosomes=args.chromosomes, num_threads=args.num_threads, spades=args.spades, age=args.age,
                        disable_assembly=args.disable_assembly,
                        svs_to_assemble=args.svs_to_assemble,
                        asm_max_size=args.spades_max_interval_size,
                        minsvlen=args.minsvlen, maxsvlen=args.maxsvlen, inswiggle=args.inswiggle,
                        enable_per_tool_output=args.enable_per_tool_output, min_support=args.min_ins_support,
                        min_support_frac=args.min_ins_support_frac, max_intervals=args.max_ins_intervals,
                        stop_spades_on_fail=args.stop_spades_on_fail, gt_window=args.gt_window,
                        gt_normal_frac=args.gt_normal_frac, isize_mean=args.isize_mean, isize_sd=args.isize_sd,
                        extraction_max_read_pairs=args.extraction_max_read_pairs, svs_to_report=args.svs_to_report,
                        min_mapq=args.min_mapq, min_avg_base_qual=args.min_avg_base_qual,
                        min_soft_clip=args.min_soft_clip, max_soft_clip=args.max_soft_clip, max_nm=args.max_nm,
                        min_matches=args.min_matches))
