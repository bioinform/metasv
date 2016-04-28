#!/usr/bin/env python

import argparse
import logging
import sys
from metasv.defaults import *
from metasv._version import __version__
from metasv.generate_sv_intervals import parallel_generate_sc_intervals

if __name__ == "__main__":
    log_format = '%(levelname)s %(asctime)-15s %(name)-20s %(message)s'
    logging.basicConfig(level=logging.INFO, format=log_format)
    logger = logging.getLogger(__name__)

    parser = argparse.ArgumentParser(
        description="Generate BED intervals for SVs using soft-clipped reads",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--bams", nargs="+", help="BAMs", required=True)
    parser.add_argument("--chromosomes", nargs="+", help="Chromosomes", default=[])
    parser.add_argument("--workdir", help="Working directory", default="work")
    parser.add_argument("--num_threads", help="Number of threads to use", default=1, type=int)
    parser.add_argument("--min_avg_base_qual", help="Minimum average base quality", default=SC_MIN_AVG_BASE_QUAL,
                        type=int)
    parser.add_argument("--min_mapq", help="Minimum MAPQ", default=SC_MIN_MAPQ, type=int)
    parser.add_argument("--min_soft_clip", help="Minimum soft-clip", default=SC_MIN_SOFT_CLIP, type=int)
    parser.add_argument("--max_nm", help="Maximum number of edits", default=SC_MAX_NM, type=int)
    parser.add_argument("--min_matches", help="Minimum number of matches", default=SC_MIN_MATCHES, type=int)
    parser.add_argument("--isize_mean", help="Insert-size mean", default=ISIZE_MEAN, type=float)
    parser.add_argument("--isize_sd", help="Insert-size s.d.", default=ISIZE_SD, type=float)
    parser.add_argument("--pad", help="Padding on both sides of the candidate locations", default=SC_PAD, type=int)
    parser.add_argument("--min_support_ins", help="Minimum read support for calling insertions using soft-clips (including neighbors)", default=MIN_SUPPORT_INS, type=int)
    parser.add_argument("--min_support_frac_ins", help="Minimum fraction of reads supporting insertion using soft-clips (including neighbors)",
                        default=MIN_SUPPORT_FRAC_INS, type=float)
    parser.add_argument("--skip_bed", help="BED regions with which no overlap should happen", type=file)
    parser.add_argument("--max_intervals",
                        help="Maximum number of intervals to process. Intervals are ranked by normalized read-support",
                        type=int, default=MAX_INTERVALS)
    parser.add_argument("--svs_to_softclip", nargs="+", help="SVs to perform soft-clip analysis on", default=SVS_SOFTCLIP_SUPPORTED,
                           choices=SVS_SOFTCLIP_SUPPORTED)
    parser.add_argument("--overlap_ratio", help="Reciprocal overlap ratio", default=OVERLAP_RATIO, type=float)
    parser.add_argument("--mean_read_length", type=float, default=MEAN_READ_LENGTH, help="Mean read length")
    parser.add_argument("--mean_read_coverage", type=float, default=MEAN_READ_COVERAGE, help="Mean read coverage")
    parser.add_argument("--min_ins_cov_frac", type=float, default=MIN_INS_COVERAGE_FRAC, help="Minimum read coverage around the insertion breakpoint.")
    parser.add_argument("--max_ins_cov_frac", type=float, default=MAX_INS_COVERAGE_FRAC, help="Maximum read coverage around the insertion breakpoint.")
    parser.add_argument("--assembly_max_tools", type=int, default=ASSEMBLY_MAX_TOOLS,
                           help="Skip assembly if more than this many tools support a call (default 1)")
    parser.add_argument("--other_scale", type=float, default=SC_OTHER_SCALE, help="Parameter to control none SV type resolution")
    parser.add_argument('--version', action='version', version='%(prog)s ' + __version__)

    args = parser.parse_args()

    logger.info("Command-line: " + " ".join(sys.argv))

    parallel_generate_sc_intervals(args.bams, args.chromosomes, args.skip_bed, args.workdir,
                                   num_threads=args.num_threads, min_avg_base_qual=args.min_avg_base_qual,
                                   min_mapq=args.min_mapq, min_soft_clip=args.min_soft_clip,
                                   pad=args.pad, min_support_ins=args.min_support_ins,
                                   min_support_frac_ins=args.min_support_frac_ins, max_intervals=args.max_intervals,
                                   max_nm=args.max_nm, min_matches=args.min_matches, isize_mean=args.isize_mean,
                                   isize_sd=args.isize_sd, svs_to_softclip=args.svs_to_softclip,
                                   overlap_ratio=args.overlap_ratio, mean_read_length=args.mean_read_length,
                                   mean_read_coverage=args.mean_read_coverage, min_ins_cov_frac=args.min_ins_cov_frac,
                                   max_ins_cov_frac=args.max_ins_cov_frac, assembly_max_tools=args.assembly_max_tools, other_scale=args.other_scale)
