import argparse
import logging
import sys
from metasv.defaults import *
from metasv.age import run_age_parallel

FORMAT = '%(levelname)s %(asctime)-15s %(name)-20s %(message)s'
logging.basicConfig(level=logging.INFO, format=FORMAT)
logger = logging.getLogger(__name__)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run AGE on files assembled under MetaSV.",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--reference", help="Reference FASTA", required=True, type=file)
    parser.add_argument("--assembly", help="Assembly FASTA", required=True, type=file)
    parser.add_argument("--age", help="Path to AGE executable", required=True, type=file)
    parser.add_argument("--work", help="Work directory", default="work")
    parser.add_argument("--pad", help="Padding to apply on both sides of the bed regions", type=int, default=AGE_PAD)
    parser.add_argument("--nthreads", help="Number of threads to use", type=int, default=1)
    parser.add_argument("--chrs", help="Chromosome list to process", nargs="+", default=[])
    parser.add_argument("--sv_types", help="SV types list to process (INS, DEL, INV)", nargs="+", default=[])
    parser.add_argument("--timeout", help="Max time for assembly processes to run", type=int, default=AGE_TIMEOUT)
    parser.add_argument("--keep_temp", help="Don't delete temporary files", action="store_true")
    parser.add_argument("--assembly_tool", help="Tool used for assembly", choices=["spades", "tigra"], default="spades")
    parser.add_argument("--min_contig_len", help="Minimum length of contig to consider", type=int,
                        default=AGE_MIN_CONTIG_LENGTH)
    parser.add_argument("--max_region_len", help="Maximum length of an SV interval", type=int,
                        default=AGE_MAX_REGION_LENGTH)
    parser.add_argument("--min_del_subalign_len", help="Minimum length of deletion sub-alginment", type=int,
                        default=MIN_DEL_SUBALIGN_LENGTH)
    parser.add_argument("--min_inv_subalign_len", help="Minimum length of inversion sub-alginment", type=int,
                        default=MIN_INV_SUBALIGN_LENGTH)
    parser.add_argument("--age_window", help="Window size for AGE to merge nearby breakpoints", type=int,
                        default=AGE_WINDOW_SIZE)
    parser.add_argument("--intervals_bed", help="BED file for assembly", type=file, required=True)

    args = parser.parse_args()

    logger.info("Command-line: {}".format(" ".join(sys.argv)))

    run_age_parallel(intervals_bed=args.intervals_bed.name, reference=args.reference.name, assembly=args.assembly.name,
                     pad=args.pad, age=args.age.name, age_workdir=args.work, timeout=args.timeout,
                     keep_temp=args.keep_temp, assembly_tool=args.assembly_tool, chrs=args.chrs, nthreads=args.nthreads,
                     min_contig_len=args.min_contig_len, max_region_len=args.max_region_len, sv_types=args.sv_types,
                     min_del_subalign_len=args.min_del_subalign_len, min_inv_subalign_len=args.min_inv_subalign_len,
                     age_window = args.age_window)