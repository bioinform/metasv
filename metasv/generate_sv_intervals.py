#!/usr/bin/env python

import sys
import os
import argparse
import multiprocessing
import logging
import collections
import itertools
import traceback
from functools import partial
import json
import base64

import pysam
import pybedtools


DEFAULT_MIN_SUPPORT = 5
DEFAULT_MIN_SUPPORT_FRAC = 0.1


def concatenate_files(files, output):
    with open(output, 'w') as outfile:
        for fname in files:
            if not os.path.isfile(fname): continue
            with open(fname) as infile:
                outfile.write(infile.read())


def get_max_soft_clip(aln):
    if aln.cigar is None:
        return 0
    return max([length for (op, length) in aln.cigar if op == 4] + [0])


def is_good_soft_clip(aln):
    sc = get_max_soft_clip(aln)
    return 20 <= sc <= 50


def is_discordant(aln):
    return 200 > abs(aln.tlen) > 0


def is_good_candidate(aln, min_avg_base_qual=20, min_mapq=5, min_soft_clip=20, max_soft_clip=50):
    if aln.is_duplicate:
        return False
    if aln.is_unmapped:
        return False
    if aln.mapq < min_mapq:
        return False
    if aln.cigar is None:
        return False

    soft_clips = [length for (op, length) in aln.cigar if op == 4]
    if len(soft_clips) != 1:
        return False
    soft_clip = soft_clips[0]

    if not (min_soft_clip <= soft_clip <= max_soft_clip):
        return False

    if aln.cigar[0][0] == 4:
        avg_base_quality = float(sum(map(ord, aln.qual[:soft_clip]))) / soft_clip
    else:
        avg_base_quality = float(sum(map(ord, aln.qual[-soft_clip:]))) / soft_clip

    return avg_base_quality - 33 >= min_avg_base_qual


def get_interval(aln, pad=500):
    start = aln.pos
    end = aln.aend

    if aln.cigar[0][0] == 4:
        return start - pad, start + pad
    return end - pad, end + pad


def merged_interval_features(feature, bam_handle):
    support_list = feature.name.split(",")
    locations = sorted(map(int, support_list[0::2]))
    num_unique_locations = len(set(locations))
    count_str = ",".join(["%s,%s" % (i, c) for (i, c) in collections.Counter(locations).items()])
    plus_support = len([i for i in support_list[1::2] if i == "+"])
    minus_support = len(locations) - plus_support
    locations_span = max(locations) - min(locations)
    name = "%s,INS,0,SC,%d,%d,%d,%d,%s" % (
    base64.b64encode(json.dumps(dict())), plus_support, minus_support, locations_span, num_unique_locations, count_str)
    interval_readcount = bam_handle.count(reference=feature.chrom, start=feature.start, end=feature.end)

    return pybedtools.Interval(feature.chrom, feature.start, feature.end, name=name, score=feature.score, otherfields=[str(interval_readcount)])


def coverage_filter(feature, bam_handle, min_support_frac=DEFAULT_MIN_SUPPORT_FRAC):
    total_count = bam_handle.count(reference=feature.chrom, start=feature.start, end=feature.end)
    return float(feature.score) >= min_support_frac * float(total_count)


def generate_sc_intervals_callback(result, result_list):
    if result is not None:
        result_list.append(result)


def generate_sc_intervals(bam, chromosome, workdir, min_avg_base_qual=20, min_mapq=5, min_soft_clip=20,
                          max_soft_clip=50, pad=500, min_support=DEFAULT_MIN_SUPPORT, max_isize=1000000000,
                          min_support_frac=DEFAULT_MIN_SUPPORT_FRAC):
    func_logger = logging.getLogger("%s-%s" % (generate_sc_intervals.__name__, multiprocessing.current_process()))

    if not os.path.isdir(workdir):
        func_logger.error("Working directory %s doesn't exist" % workdir)
        return None

    func_logger.info("Generating candidate insertion intervals from %s for chromsome %s" % (bam, chromosome))
    pybedtools.set_tempdir(workdir)

    unmerged_intervals = []
    try:
        sam_file = pysam.Samfile(bam, "rb")
        for aln in sam_file.fetch(reference=chromosome):
            if abs(aln.tlen) > max_isize:
                continue
            if not is_good_candidate(aln, min_avg_base_qual=min_avg_base_qual, min_mapq=min_mapq,
                                     min_soft_clip=min_soft_clip, max_soft_clip=max_soft_clip): continue
            interval = get_interval(aln, pad=pad)
            soft_clip_location = sum(interval) / 2
            strand = "-" if aln.is_reverse else "+"
            name = "%d,%s" % (soft_clip_location, strand)
            unmerged_intervals.append(
                pybedtools.Interval(chromosome, interval[0], interval[1], name=name, score="1", strand=strand))
        sam_file.close()

        if not unmerged_intervals:
            sam_file.close()
            func_logger.warn("No intervals generated")
            return None

        unmerged_bed = os.path.join(workdir, "unmerged.bed")
        bedtool = pybedtools.BedTool(unmerged_intervals).sort().moveto(unmerged_bed)
        func_logger.info("%d candidate reads" % (bedtool.count()))

        merged_bed = os.path.join(workdir, "merged.bed")
        bedtool = bedtool.merge(c="4,5", o="collapse,sum", d=-500).moveto(merged_bed)
        func_logger.info("%d merged intervals" % (bedtool.count()))

        filtered_bed = os.path.join(workdir, "filtered_merged.bed")
        bedtool = bedtool.filter(lambda x: int(x.score) >= min_support).each(partial(merged_interval_features, bam_handle=sam_file)).moveto(
            filtered_bed)
        func_logger.info("%d filtered intervals" % (bedtool.count()))

        # Now filter based on coverage
        coverage_filtered_bed = os.path.join(workdir, "coverage_filtered_merged.bed")
        bedtool = bedtool.filter(lambda x: float(x.fields[6]) * min_support_frac <= float(x.score)).moveto(coverage_filtered_bed)
        func_logger.info("%d coverage filtered intervals" % (bedtool.count()))

        sam_file.close()
    except Exception as e:
        func_logger.error('Caught exception in worker thread')

        # This prints the type, value, and stack trace of the
        # current exception being handled.
        traceback.print_exc()

        print()
        raise e

    pybedtools.cleanup(remove_all=True)

    return coverage_filtered_bed


def parallel_generate_sc_intervals(bams, chromosomes, skip_bed, workdir, num_threads=1, min_avg_base_qual=20,
                                   min_mapq=5, min_soft_clip=20, max_soft_clip=50, pad=500,
                                   min_support=DEFAULT_MIN_SUPPORT, min_support_frac=DEFAULT_MIN_SUPPORT_FRAC):
    func_logger = logging.getLogger(
        "%s-%s" % (parallel_generate_sc_intervals.__name__, multiprocessing.current_process()))

    if not os.path.isdir(workdir):
        func_logger.info("Creating directory %s" % workdir)
        os.makedirs(workdir)

    if not chromosomes:
        func_logger.info("Chromosome list unspecified. Inferring from the BAMs")
        for bam in bams:
            bamfile = pysam.Samfile(bam, "rb")
            chromosomes += list(bamfile.references)
            bamfile.close()
        chromosomes = sorted(list(set(chromosomes)))
        func_logger.info("Chromosome list inferred as %s" % (str(chromosomes)))

    if not chromosomes:
        func_logger.error("Chromosome list empty")
        return None

    pool = multiprocessing.Pool(num_threads)

    bed_files = []
    for index, (bam, chromosome) in enumerate(itertools.product(bams, chromosomes)):
        process_workdir = os.path.join(workdir, str(index))
        if not os.path.isdir(process_workdir):
            os.makedirs(process_workdir)

        args_list = [bam, chromosome, process_workdir]
        kwargs_dict = {"min_avg_base_qual": min_avg_base_qual, "min_mapq": min_mapq, "min_soft_clip": min_soft_clip,
                       "max_soft_clip": max_soft_clip, "pad": pad, "min_support": min_support,
                       "min_support_frac": min_support_frac}
        pool.apply_async(generate_sc_intervals, args=args_list, kwds=kwargs_dict,
                         callback=partial(generate_sc_intervals_callback, result_list=bed_files))

    pool.close()
    pool.join()

    func_logger.info("Following BED files will be merged: %s" % (str(bed_files)))

    if not bed_files:
        func_logger.warn("No intervals generated")
        return None

    pybedtools.set_tempdir(workdir)
    bedtool = pybedtools.BedTool(bed_files[0])

    for bed_file in bed_files[1:]:
        bedtool = bedtool.cat(pybedtools.BedTool(bed_file), postmerge=False)

    if skip_bed:
        skip_bedtool = pybedtools.BedTool(skip_bed)
        func_logger.info(
            "Merging %d features with %d features from %s" % (bedtool.count(), skip_bedtool.count(), skip_bed))
        bedtool = skip_bedtool.cat(bedtool, postmerge=False).sort()
        func_logger.info("After merging with %s %d features" % (skip_bed, bedtool.count()))

    bedtool = bedtool.moveto(os.path.join(workdir, "intervals.bed"))
    pybedtools.cleanup(remove_all=True)

    return bedtool.fn


if __name__ == "__main__":
    FORMAT = '%(levelname)s %(asctime)-15s %(name)-20s %(message)s'
    logging.basicConfig(level=logging.INFO, format=FORMAT)
    logger = logging.getLogger(__name__)

    parser = argparse.ArgumentParser(
        description="Generate BED intervals for insertion detection using soft-clipped reads",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--bams", nargs="+", help="BAMs", required=True)
    parser.add_argument("--chromosomes", nargs="+", help="Chromosomes", default=[])
    parser.add_argument("--workdir", help="Working directory", default="work")
    parser.add_argument("--num_threads", help="Number of threads to use", default=1, type=int)
    parser.add_argument("--min_avg_base_qual", help="Minimum average base quality", default=20, type=int)
    parser.add_argument("--min_mapq", help="Minimum MAPQ", default=5, type=int)
    parser.add_argument("--min_soft_clip", help="Minimum soft-clip", default=20, type=int)
    parser.add_argument("--max_soft_clip", help="Maximum soft-clip", default=50, type=int)
    parser.add_argument("--pad", help="Padding on both sides of the candidate locations", default=500, type=int)
    parser.add_argument("--min_support", help="Minimum supporting reads", default=DEFAULT_MIN_SUPPORT, type=int)
    parser.add_argument("--min_support_frac", help="Minimum fraction of total reads for interval",
                        default=DEFAULT_MIN_SUPPORT_FRAC, type=float)
    parser.add_argument("--skip_bed", help="BED regions with which no overlap should happen", type=file)

    args = parser.parse_args()

    logger.info("Command-line: " + " ".join(sys.argv))

    parallel_generate_sc_intervals(args.bams, args.chromosomes, args.skip_bed, args.workdir,
                                   num_threads=args.num_threads, min_avg_base_qual=args.min_avg_base_qual,
                                   min_mapq=args.min_mapq, min_soft_clip=args.min_soft_clip,
                                   max_soft_clip=args.max_soft_clip, pad=args.pad, min_support=args.min_support,
                                   min_support_frac=args.min_support_frac)
