import os
import argparse
import logging
import multiprocessing
import subprocess
import fileinput
import traceback
from functools import partial, update_wrapper
import json
import base64
import time

import pysam
import pybedtools

import extract_pairs
from defaults import *

FORMAT = '%(levelname)s %(asctime)-15s %(name)-20s %(message)s'
logging.basicConfig(level=logging.INFO, format=FORMAT)
logger = logging.getLogger(__name__)

precise_methods = set(["AS", "SR", "JM"])


def run_cmd(cmd, logger, spades_log_fd):
    logger.info("Running command %s" % cmd)
    spades_log_fd.write("*************************************************\n")
    start_time = time.time()
    retcode = subprocess.call(cmd, shell=True, stderr=spades_log_fd, stdout=spades_log_fd)
    logger.info("Returned code %d (%g seconds)" % (retcode, time.time() - start_time))

    return retcode


def append_contigs(src, interval, dst_fd, fn_id=0, sv_type="INS"):
    with open(src, "r") as fd:
        for line in fd:
            if line.startswith(">"):
                dst_fd.write(
                    ">%s_%d_%d_%s_%d_%s" % (interval.chrom, interval.start, interval.end, sv_type, fn_id, line[1:]))
            else:
                dst_fd.write(line)


def run_spades_single(intervals=[], bam=None, spades=None, work=None, pad=SPADES_PAD, timeout=SPADES_TIMEOUT,
                      isize_min=ISIZE_MIN,
                      isize_max=ISIZE_MAX, stop_on_fail=False, max_read_pairs=EXTRACTION_MAX_READ_PAIRS):
    thread_logger = logging.getLogger("%s-%s" % (run_spades_single.__name__, multiprocessing.current_process()))

    if not os.path.isdir(work):
        thread_logger.info("Creating %s" % work)
        os.makedirs(work)

    merged_contigs = open(os.path.join(work, "merged.fa"), "w")
    spades_log_fd = open(os.path.join(work, "spades.log"), "w")

    extract_fns = [extract_pairs.all_pair, extract_pairs.non_perfect]

    try:
        for interval in intervals:
            region = "%s:%d-%d" % (interval.chrom, interval.start, interval.end)
            thread_logger.info("Processing interval %s" % (str(interval).strip()))

            sv_type = interval.name.split(",")[1]

            extraction_counts = extract_pairs.extract_read_pairs(bam, region, "%s/" % work, extract_fns, pad=pad,
                                                                 max_read_pairs=max_read_pairs, sv_type=sv_type)
            all_pair_count = extraction_counts[0][1]

            for fn_id, ((end1, end2), extracted_count) in enumerate(extraction_counts):
                extract_fn_name = extract_fns[fn_id].__name__

                if fn_id > 0 and extracted_count == all_pair_count:
                    thread_logger.info("Skipping assembly from %s since read count same as all_pairs" % extract_fn_name)
                    continue

                if extracted_count >= 5:
                    extra_opt = "--sc" if not fn_id == 0 else ""
                    spades_log_fd.write("Running spades for interval %s with extraction function %s\n" % (
                    str(interval).strip(), extract_fn_name))
                    retcode = run_cmd(
                        "bash -c \"timeout -k 1 %ds %s -1 %s -2 %s -o %s/spades_%s/ -m 4 -t 1 --phred-offset 33 %s\"" % (
                            timeout, spades, end1, end2, work, extract_fn_name, extra_opt), thread_logger,
                        spades_log_fd)
                    if retcode == 0:
                        append_contigs(os.path.join(work, "spades_%s/contigs.fasta") % extract_fn_name, interval,
                                       merged_contigs, fn_id, sv_type)
                    elif retcode == 1:
                        thread_logger.error("Spades failed for some reason")
                        if stop_on_fail:
                            thread_logger.error("Aborting!")
                            raise Exception("Spades failure on interval %s for extraction function %s\n" % (
                            str(interval).strip(), extract_fn_name))
                else:
                    thread_logger.info("Too few read pairs (%d) extracted. Skipping assembly." % extracted_count)
    except Exception as e:
        thread_logger.error('Caught exception in worker thread')

        # This prints the type, value, and stack trace of the
        # current exception being handled.
        traceback.print_exc()

        print()
        raise e

    merged_contigs.close()

    return os.path.abspath(merged_contigs.name)


def run_spades_single_callback(result, result_list):
    if result is not None:
        result_list.append(result)


def should_be_assembled(interval, max_interval_size=SPADES_MAX_INTERVAL_SIZE, svs_to_assemble=SVS_ASSEMBLY_SUPPORTED):   
    if interval.length > max_interval_size: 
        logger.info("Will not assemble (%s). Too large interval length: %d > %d" % (interval.name, interval.length,max_interval_size))
        return False
    # TODO: fix this later to make MetaSV do the right thing
    should_assemble = False
    for supported_sv in svs_to_assemble:
        if interval.name.find(supported_sv) >= 0:
            should_assemble = True
    if not should_assemble: return False

    name_fields = interval.name.split(",")
    methods = set(name_fields[3].split(";"))

    
    return len(methods) == 1 or not (methods & precise_methods)


def shouldnt_be_assembled(interval, max_interval_size=SPADES_MAX_INTERVAL_SIZE, svs_to_assemble=SVS_ASSEMBLY_SUPPORTED):
    return not should_be_assembled(interval, max_interval_size=max_interval_size,
                                   svs_to_assemble=svs_to_assemble)


def add_breakpoints(interval):
    fields = interval.fields
    name_fields = interval.name.split(",")
    methods = set(name_fields[3].split(";"))
    breakpoints = [-1, -1]
    if len(methods) > 1 and methods & precise_methods:
        breakpoints = [interval.start, interval.end]
    fields += map(str, breakpoints)
    fields += name_fields[2:3]
    fields.append(base64.b64encode(json.dumps(dict())))  # Does nothing, make sure the fields line up
    return pybedtools.create_interval_from_list(fields)


def run_spades_parallel(bam=None, spades=None, bed=None, work=None, pad=SPADES_PAD, nthreads=1, chrs=[],
                        max_interval_size=SPADES_MAX_INTERVAL_SIZE,
                        timeout=SPADES_TIMEOUT, isize_min=ISIZE_MIN, isize_max=ISIZE_MAX,
                        svs_to_assemble=SVS_ASSEMBLY_SUPPORTED,
                        stop_on_fail=False, max_read_pairs=EXTRACTION_MAX_READ_PAIRS):
    pybedtools.set_tempdir(work)

    logger.info("Running SPAdes on the intervals in %s" % bed)
    if not bed:
        logger.info("No BED file specified")
        return None, None

    bedtool = pybedtools.BedTool(bed)
    total = bedtool.count()

    chrs = set(chrs)
    all_intervals = [interval for interval in bedtool] if not chrs else [interval for interval in bedtool if
                                                                         interval.chrom in chrs]


    selected_intervals = filter(partial(should_be_assembled, max_interval_size=max_interval_size, svs_to_assemble=svs_to_assemble),
                                all_intervals)
    ignored_intervals = filter(partial(shouldnt_be_assembled, max_interval_size=max_interval_size, svs_to_assemble=svs_to_assemble),
                               all_intervals)

    logger.info("%d intervals selected" % len(selected_intervals))
    logger.info("%d intervals ignored" % len(ignored_intervals))

    pool = multiprocessing.Pool(nthreads)
    assembly_fastas = []
    for i in xrange(nthreads):
        intervals = [interval for (j, interval) in enumerate(selected_intervals) if (j % nthreads) == i]
        kwargs_dict = {"intervals": intervals, "bam": bam, "spades": spades, "work": "%s/%d" % (work, i), "pad": pad,
                       "timeout": timeout, "isize_min": isize_min, "isize_max": isize_max, "stop_on_fail": stop_on_fail,
                       "max_read_pairs": max_read_pairs}
        pool.apply_async(run_spades_single, kwds=kwargs_dict,
                         callback=partial(run_spades_single_callback, result_list=assembly_fastas))

    pool.close()
    pool.join()

    logger.info("Merging the contigs from %s" % (str(assembly_fastas)))
    assembled_fasta = os.path.join(work, "spades_assembled.fa")
    with open(assembled_fasta, "w") as assembled_fd:
        for line in fileinput.input(assembly_fastas):
            assembled_fd.write("%s\n" % (line.strip()))

    if os.path.getsize(assembled_fasta) > 0:
        logger.info("Indexing the assemblies")
        pysam.faidx(assembled_fasta)
    else:
        logger.error("No assembly generated")
        assembled_fasta = None

    ignored_bed = None
    if ignored_intervals:
        ignored_bed = os.path.join(work, "ignored.bed")
        pybedtools.BedTool(ignored_intervals).each(add_breakpoints).saveas(ignored_bed)

    pybedtools.cleanup(remove_all=True)

    return assembled_fasta, ignored_bed


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run spades on a bed file.",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--bam", help="BAM file to use reads from", required=True)
    parser.add_argument("--spades", help="Spades python executable", required=True)
    parser.add_argument("--work", help="Work directory", default="work")
    parser.add_argument("--bed", help="BED file for assembly regions", required=True)
    parser.add_argument("--pad", help="Padding to apply on both sides of the bed regions", type=int, default=SPADES_PAD)
    parser.add_argument("--nthreads", help="Number of threads to use", type=int, default=1)
    parser.add_argument("--chrs", help="Chromosome list to process", nargs="+", default=[])
    parser.add_argument("--timeout", help="Max time for assembly processes to run", type=int, default=SPADES_TIMEOUT)
    parser.add_argument("--max_interval_size", help="Maximum size of interval to process", type=int,
                        default=SPADES_MAX_INTERVAL_SIZE)
    parser.add_argument("--isize_min", help="Minimum insert size for normal pair", type=int, default=ISIZE_MIN)
    parser.add_argument("--isize_max", help="Maximum insert size for normal pair", type=int, default=ISIZE_MAX)
    parser.add_argument("--max_read_pairs", help="Maximum read pairs to use for assembly", type=int,
                        default=EXTRACTION_MAX_READ_PAIRS)

    args = parser.parse_args()

    if not os.path.isdir(args.work):
        logger.info("Creating %s" % args.work)
        os.makedirs(args.work)

    run_spades_parallel(bam=args.bam, spades=args.spades, bed=args.bed, work=args.work, pad=args.pad,
                        nthreads=args.nthreads, chrs=args.chrs, max_interval_size=args.max_interval_size,
                        isize_min=args.isize_min, isize_max=args.isize_max, max_read_pairs=args.max_read_pairs)
