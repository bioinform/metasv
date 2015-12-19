import traceback
import os
import argparse
import multiprocessing
import subprocess
import hashlib
from functools import partial
import json
import base64

import pybedtools

from age_parser import *
from defaults import *

from generate_final_vcf import convert_metasv_bed_to_vcf
from genotype import parallel_genotype_intervals

FORMAT = '%(levelname)s %(asctime)-15s %(name)-20s %(message)s'
logging.basicConfig(level=logging.INFO, format=FORMAT)
logger = logging.getLogger(__name__)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run AGE on files assembled under MetaSV.",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--reference", help="Reference FASTA", required=True, type=file)
    parser.add_argument("--breakpoints_bed", help="AGE bed output", required=True, type=file)
    parser.add_argument("--ignored_bed", help="SPADE ignored_bed output", type=file, default= None)
    parser.add_argument("--age", help="Path to AGE executable", required=True, type=file)
    parser.add_argument("--work", help="Work directory", default="work")
    parser.add_argument("--nthreads", help="Number of threads to use", type=int, default=1)
    parser.add_argument("--chrs", help="Chromosome list to process", nargs="+", default=[])
    parser.add_argument("--sample", metavar="Sample", help="Sample name", required=True)
    parser.add_argument("--isize_mean", type=float, default=ISIZE_MEAN, help="Insert size mean")
    parser.add_argument("--isize_sd", type=float, default=ISIZE_SD, help="Insert size standard deviation")
    parser.add_argument("--gt_normal_frac", type=float, default=GT_NORMAL_FRAC,
                           help="Min. fraction of reads supporting reference for genotyping")
    parser.add_argument("--bam", help="BAM", type=file)

    args = parser.parse_args()


    breakpoints_bed=args.breakpoints_bed.name
    if args.ignored_bed:
        ignored_bed=args.ignored_bed.name
    else:
        ignored_bed=None
    final_vcf = os.path.join(args.work, "variants.vcf")

    final_bed = os.path.join(args.work, "final.bed")
    if breakpoints_bed:
        if ignored_bed:
            pybedtools.BedTool(breakpoints_bed) \
                .cat(pybedtools.BedTool(ignored_bed), postmerge=False) \
                .sort().saveas(final_bed)
        else:
            pybedtools.BedTool(breakpoints_bed).saveas(final_bed)
    elif ignored_bed:
        pybedtools.BedTool(ignored_bed).sort().saveas(final_bed)
    else:
        final_bed = None

    genotyped_bed = parallel_genotype_intervals(final_bed, args.bam.name,
                                                workdir=os.path.join(args.work, "genotyping"),
                                                nthreads=args.nthreads, chromosomes=args.chrs,
                                                 isize_mean=args.isize_mean,
                                                isize_sd=args.isize_sd,
                                                normal_frac_threshold=args.gt_normal_frac)
    
    genotyped_bed = os.path.join(args.work, "genotyping","genotyped.bed")
    logger.info("Output final VCF file")

    convert_metasv_bed_to_vcf(bedfile=genotyped_bed, vcf_out=final_vcf, workdir=args.work, sample=args.sample, pass_calls=False)

    logger.info("Clean up pybedtools")

    pybedtools.cleanup(remove_all=True)

    logger.info("All Done!")
