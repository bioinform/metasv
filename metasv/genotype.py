import logging
import os
import multiprocessing
import traceback
import argparse
import sys
import time
import json
import base64
from functools import partial
import pybedtools
import pysam

from defaults import ISIZE_MEAN, ISIZE_SD, GT_WINDOW, GT_NORMAL_FRAC

GT_HET = "0/1"
GT_HOM = "1/1"
GT_REF = "0/0"
GT_UNK = "./."


def count_reads_supporting_ref(chrom, start, end, bam_handle, isize_min, isize_max, window):
    total_normal_reads = 0
    total_read_bases = 0
    total_reads = 0
    window_start = max(0, start - window)
    window_end = end + window
    for aln in bam_handle.fetch(chrom, window_start, window_end):
        if aln.is_duplicate or not aln.is_paired:
            continue
        total_reads += 1
        if aln.is_unmapped or aln.mate_is_unmapped:
            continue
        if aln.rnext != aln.tid: continue
        if aln.is_reverse:
            if not (aln.pnext < aln.pos and not aln.mate_is_reverse): continue
        else:
            if not (aln.pnext > aln.pos and aln.mate_is_reverse): continue
        if not (((aln.aend - end) >= 20 and (end - aln.pos) >= 20) or ((start - aln.pos) >= 20 and (aln.aend - start) >= 20)):
            continue
        tlen = abs(aln.tlen)
        if isize_min <= tlen <= isize_max:
            total_normal_reads += 1
            total_read_bases = total_read_bases + aln.qlen
    return total_normal_reads, total_read_bases, total_reads


def genotype_interval(chrom, start, end, sv_type, sv_length, bam_handle, isize_min, isize_max, window=GT_WINDOW, normal_frac_threshold=GT_NORMAL_FRAC):
    func_logger = logging.getLogger("%s-%s" % (genotype_interval.__name__, multiprocessing.current_process()))

    locations = [start, end] if sv_type != "INS" else [start]
    total_normal, total = 0, 0
    for location in locations:
        total_normal_, total_bases_, total_ = count_reads_supporting_ref(chrom, location, location, bam_handle, isize_min, isize_max, window)
        total_normal += total_normal_
        total += total_

    normal_frac = float(total_normal) / float(max(1, total))
    gt = GT_REF
    if normal_frac < 1 - normal_frac_threshold:
        gt = GT_HET if normal_frac >= normal_frac_threshold else GT_HOM

    func_logger.info("For interval %s:%d-%d %s counts are %d, %d and normal_frac is %g gt is %s" % (chrom, start, end, sv_type, total_normal, total, normal_frac, gt))
    return gt

def parse_interval(interval):
    chrom = interval.chrom
    pos = interval.start
    end = interval.end

    sub_names = interval.name.split(":")
    sub_lengths = map(lambda x: int(x.split(",")[2]), sub_names)

    sub_types = map(lambda x: x.split(",")[1], sub_names)
    sub_methods = [name.split(",")[3] for name in sub_names]
    try:
        info = json.loads(base64.b64decode(name.split(",")[0]))
    except TypeError:
        info = dict()
    if len(interval.fields) > 9:
        info.update(json.loads(base64.b64decode(interval.fields[9])))

    index_to_use = 0
    svlen = -1
    if "DEL" in sub_types:
        index_to_use = sub_types.index("DEL")
    elif "INV" in sub_types:
        index_to_use = sub_types.index("INV")
    elif "INS" in sub_types and "SC" in sub_methods:
        index_to_use = sub_methods.index("SC")
        pos = int(interval.fields[6])
        end = int(interval.fields[7])
        svlen = int(interval.fields[8])

    if svlen < 0: svlen = sub_lengths[index_to_use]
    if sub_types[index_to_use] == "DEL":
        svlen = -svlen
    sv_type = sub_types[index_to_use]

    return chrom, pos, end, sv_type, svlen


def genotype_intervals_callback(result, result_list):
    if result is not None:
        result_list.append(result)


def genotype_intervals(intervals_file=None, bam=None, workdir=None, window=GT_WINDOW, isize_mean=ISIZE_MEAN, isize_sd=ISIZE_SD, normal_frac_threshold=GT_NORMAL_FRAC):
    func_logger = logging.getLogger("%s-%s" % (genotype_intervals.__name__, multiprocessing.current_process()))

    if workdir and not os.path.isdir(workdir):
        os.makedirs(workdir)

    pybedtools.set_tempdir(workdir)

    genotyped_intervals = []
    start_time = time.time()

    isize_min = max(0, isize_mean - 3 * isize_sd)
    isize_max = isize_mean + 3 * isize_sd

    try:
        bam_handle = pysam.Samfile(bam, "rb")
        for interval in pybedtools.BedTool(intervals_file):
            chrom, start, end, sv_type, svlen = parse_interval(interval)
            genotype = genotype_interval(chrom, start, end, sv_type, svlen, bam_handle, isize_min, isize_max, window, normal_frac_threshold)
            fields = interval.fields + [genotype]
            genotyped_intervals.append(pybedtools.create_interval_from_list(fields))
        bedtool = pybedtools.BedTool(genotyped_intervals).moveto(os.path.join(workdir, "genotyped.bed"))
    except Exception as e:
        func_logger.error('Caught exception in worker thread')

        # This prints the type, value, and stack trace of the
        # current exception being handled.
        traceback.print_exc()

        print()
        raise e
    func_logger.info("Genotyped %d intervals in %g minutes" % (len(genotyped_intervals), (time.time() - start_time)/60.0))

    return bedtool.fn

def parallel_genotype_intervals(intervals_file=None, bam=None, workdir=None, nthreads=1, chromosomes=[], window=GT_WINDOW, isize_mean=ISIZE_MEAN, isize_sd=ISIZE_SD, normal_frac_threshold=GT_NORMAL_FRAC):
    func_logger = logging.getLogger("%s-%s" % (parallel_genotype_intervals.__name__, multiprocessing.current_process()))
    if workdir and not os.path.isdir(workdir):
        os.makedirs(workdir)

    chromosomes = set(chromosomes)

    start_time = time.time()

    bedtool = pybedtools.BedTool(intervals_file)
    selected_intervals = [interval for interval in bedtool if not chromosomes or interval.chrom in chromosomes]
    nthreads = min(len(selected_intervals), nthreads)
    intervals_per_process = (len(selected_intervals) + nthreads - 1) / nthreads

    pool = multiprocessing.Pool(nthreads)
    genotyped_beds = []
    for i in xrange(nthreads):
        process_workdir = os.path.join(workdir, str(i))
        if not os.path.isdir(process_workdir):
            os.makedirs(process_workdir)
        process_intervals = pybedtools.BedTool(selected_intervals[i*intervals_per_process: (i+1)*intervals_per_process]).saveas(os.path.join(process_workdir, "ungenotyped.bed"))
        kwargs_dict = {"intervals_file": process_intervals.fn, "bam": bam, "workdir": process_workdir, "window": window, "isize_mean": isize_mean, "isize_sd": isize_sd, "normal_frac_threshold": normal_frac_threshold}
        pool.apply_async(genotype_intervals, kwds=kwargs_dict, callback=partial(genotype_intervals_callback, result_list=genotyped_beds))

    pool.close()
    pool.join()

    func_logger.info("Following BED files will be merged: %s" % (str(genotyped_beds)))

    if not genotyped_beds:
        func_logger.warn("No intervals generated")
        return None

    pybedtools.set_tempdir(workdir)
    bedtool = pybedtools.BedTool(genotyped_beds[0])

    for bed_file in genotyped_beds[1:]:
        bedtool = bedtool.cat(pybedtools.BedTool(bed_file), postmerge=False)
    bedtool = bedtool.sort().moveto(os.path.join(workdir, "genotyped.bed"))

    func_logger.info("Finished parallel genotyping of %d intervals in %g minutes" % (len(selected_intervals), (time.time() - start_time)/60.0))

    return bedtool.fn


if __name__ == "__main__":
    FORMAT = '%(levelname)s %(asctime)-15s %(name)-20s %(message)s'
    logging.basicConfig(level=logging.INFO, format=FORMAT)
    logger = logging.getLogger(__name__)

    parser = argparse.ArgumentParser(
        description="Genotype final BED output from MetaSV assembly",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--bam", help="BAM", required=True, type=file)
    parser.add_argument("--chromosomes", nargs="+", help="Chromosomes to process. Leave unspecified to process all intervals.", default=[])
    parser.add_argument("--workdir", help="Working directory", default="work")
    parser.add_argument("--nthreads", help="Number of threads to use", default=1, type=int)
    parser.add_argument("--intervals_file", help="Final BED output from MetaSV assembly", required=True, type=file)
    parser.add_argument("--window", help="Window to use for genotyping", default=GT_WINDOW, type=int)
    parser.add_argument("--isize_mean", help="Insert size mean", default=ISIZE_MEAN, type=float)
    parser.add_argument("--isize_sd", help="Insert size standard deviation", default=ISIZE_SD, type=float)
    parser.add_argument("--normal_frac", help="Minimum fraction of normal reads to call heterozygous", default=GT_NORMAL_FRAC, type=float)

    args = parser.parse_args()

    logger.info("Command-line: " + " ".join(sys.argv))

    genotyped_bed = parallel_genotype_intervals(args.intervals_file.name, args.bam.name, args.workdir, args.nthreads, args.chromosomes, args.window, args.isize_mean, args.isize_sd, args.normal_frac)
    if genotyped_bed:
        logger.info("Generated genotyped BED as %s" % genotyped_bed)
        sys.exit(os.EX_OK)
    else:
        logger.error("No genotyped BED generated")
        sys.exit(os.EX_DATAERR)