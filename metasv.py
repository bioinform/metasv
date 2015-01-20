#!/usr/bin/python

import sys
import argparse
from metasv.run_metasv import run_metasv

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
parser.add_argument("--bam", help = "BAM", type = file)
parser.add_argument("--chromosomes", help = "Chromosome list to process. If unspecified, then all chromosomes will be considered.", nargs = "+", default = [])
parser.add_argument("--num_threads", help = "Number of threads to use", type = int, default = 1)
parser.add_argument("--outdir", help = "Output directory", required = True)
parser.add_argument("--spades", help = "Path to SPAdes executable", required = False)
parser.add_argument("--age", help = "Path to AGE executable", required = False)
parser.add_argument("--disable_assembly", action = "store_true", help = "Disable assembly")

args = parser.parse_args()

sys.exit(run_metasv(args.sample, args.reference, pindel_vcf = args.pindel_vcf, pindel_native = args.pindel_native, breakdancer_vcf = args.breakdancer_vcf, breakdancer_native = args.breakdancer_native,
        breakseq_vcf = args.breakseq_vcf, cnvnator_vcf = args.cnvnator_vcf, cnvnator_native = args.cnvnator_native, gatk_vcf = args.gatk_vcf,
        gaps = args.gaps, filter_gaps = args.filter_gaps, keep_standard_contigs = args.keep_standard_contigs, wiggle = args.wiggle, overlap_ratio = args.overlap_ratio,
        workdir = args.workdir, outdir = args.outdir, boost_ins = args.boost_ins, bam = args.bam, chromosomes = args.chromosomes, num_threads = args.num_threads, spades = args.spades, age = args.age, disable_assembly = args.disable_assembly))
