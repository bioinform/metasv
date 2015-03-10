from __future__ import print_function
import logging
import multiprocessing
from collections import defaultdict

import pybedtools


logger = logging.getLogger(__name__)


def pair_intervals(intervals, reference_length, min_interval_length=200, window=50):
    func_logger = logging.getLogger("%s-%s" % (pair_intervals.__name__, multiprocessing.current_process()))

    interval_pairs = []
    overlap_intervals = set()
    for interval1 in intervals:
        interval1_fixed = (min(interval1), max(interval1))
        if interval1_fixed[1] - interval1_fixed[0] < min_interval_length: continue
        for interval2 in intervals:
            if interval1 == interval2: continue
            interval2_fixed = (min(interval2), max(interval2))
            if interval1_fixed[0] >= interval2_fixed[0]: continue
            if interval2_fixed[1] - interval2_fixed[0] < min_interval_length: continue

            overlap_interval = (
                max(interval1_fixed[0], interval2_fixed[0]), min(interval1_fixed[1], interval2_fixed[1]))
            overlap_interval = (min(overlap_interval), max(overlap_interval))

            if overlap_interval[1] - overlap_interval[0] <= 20 and overlap_interval[0] >= window and overlap_interval[
                1] <= window + reference_length:
                interval_pairs.append((interval1_fixed, interval2_fixed))
                overlap_intervals.add(overlap_interval)

    interval_pairs = list(set(interval_pairs))
    func_logger.info("Paired %s as %s" % (str(intervals), str(interval_pairs)))
    func_logger.info("Overlap intervals %s" % (str(overlap_intervals)))

    if not overlap_intervals:
        return 0, -1, -1

    low = min([interval[0] for interval in overlap_intervals])
    high = max([interval[1] for interval in overlap_intervals])

    func_logger.info("low = %d, high = %d overlap_intervals = %s" % (low, high, str(overlap_intervals)))
    if high - low <= 20:
        func_logger.info("Pairing success! The overlap intervals are close enough.")
        return len(overlap_intervals), low, high
    return len(overlap_intervals), -1, -1


def get_insertion_breakpoints(age_records, intervals, window=20, start=0):
    func_logger = logging.getLogger("%s-%s" % (get_insertion_breakpoints.__name__, multiprocessing.current_process()))
    bedtools_intervals = [pybedtools.Interval("1", interval[0], interval[1]) for interval in sorted(intervals)]
    func_logger.info("bedtools_intervals %s" % (str(bedtools_intervals)))
    if not bedtools_intervals:
        return []

    potential_breakpoints = sorted(list(set(
        [interval.start for interval in bedtools_intervals] + [interval.end for interval in bedtools_intervals])))

    breakpoints = []
    for breakpoint in potential_breakpoints[1:-1]:
        func_logger.info("\tExamining potential breakpoint %d for support" % breakpoint)
        left_support = [interval[0] for interval in intervals if abs(interval[0] - breakpoint) <= window]
        right_support = [interval[1] for interval in intervals if abs(interval[1] - breakpoint) <= window]
        counter_examples = [age_record for age_record in age_records if age_record.has_long_ref_flanks() and (
            age_record.has_ref_deletion(window) or age_record.has_insertion(min_diff=1,
                                                                            max_diff=49)) and age_record.breakpoint_match(
            breakpoint, window)]
        if counter_examples:
            counter_example_ends = [age_record.start1_end1s for age_record in counter_examples]
            func_logger.info("\tSkipping breakpoint %d due to %s" % (breakpoint, str(counter_example_ends)))
            continue

        if left_support:
            func_logger.info("\tLeft support %s" % (str(left_support)))
        if right_support:
            func_logger.info("\tRight support %s" % (str(left_support)))

        if (left_support and right_support) and min(
                        [window + 1] + [abs(b[0] - breakpoint) for b in breakpoints]) > window:
            both_support = [age_record for age_record in age_records if
                            age_record.has_insertion(min_diff=50, max_diff=1000000000) and age_record.breakpoint_match(
                                breakpoint, window)]
            func_logger.info("\tboth_support = %s" % (str(both_support)))
            func_logger.info(
                "\tinsertion lengths = %s" % (str([age_record.insertion_length() for age_record in both_support])))
            insertion_length = max([0] + [age_record.insertion_length() for age_record in both_support])
            func_logger.info("\tInsertion length = %d" % insertion_length)
            breakpoints.append((breakpoint, insertion_length))

    func_logger.info("Gathered breakpoints as %s" % (str(breakpoints)))

    return [(start + b[0], b[1]) for b in breakpoints]


def get_deletion_breakpoints(age_records, window=20, min_flank_length=50, start=0):
    func_logger = logging.getLogger("%s-%s" % (get_deletion_breakpoints.__name__, multiprocessing.current_process()))

    potential_breakpoints = sorted(
        [age_record.start1_end1s[0][1] for age_record in age_records] + [age_record.start1_end1s[1][0] for age_record in
                                                                         age_records])
    breakpoints = []
    for breakpoint in potential_breakpoints:
        left_support = [age_record for age_record in age_records if
                        abs(age_record.start1_end1s[0][1] - breakpoint) < window]
        right_support = [age_record for age_record in age_records if
                         abs(age_record.start1_end1s[1][0] - breakpoint) < window]

        if (left_support or right_support) and min([window + 1] + [abs(b - breakpoint) for b in breakpoints]) >= window:
            breakpoints.append(breakpoint)

    func_logger.info("Gathered breakpoints as %s" % (str(breakpoints)))

    return [start + breakpoint for breakpoint in breakpoints]


def overlap(s1, e1, s2, e2):
    return min(e1, e2) - max(s1, s2) + 1


def are_positions_consistent(assemblies):
    if not assemblies:
        return 0
    low = min([assembly.start for assembly in assemblies])
    high = max([assembly.end for assembly in assemblies])

    return (len(assemblies) > 1) and (high - low) <= 20


def are_overlap_intervals_close(overlaps):
    low = min([overlap[0] for overlap in overlaps])
    high = max([overlap[1] for overlap in overlaps])

    return abs(high - low) <= 20


def get_overlap_interval(intervals):
    low = max([interval[0] for interval in intervals])
    high = min([interval[1] for interval in intervals])

    return min(low, high), max(low, high)


def pair_deletion_breakpoints(intervals, reference_length, window=50):
    interval_pairs = []
    overlap_intervals = set()
    for interval1 in intervals:
        interval1_fixed = (min(interval1), max(interval1))
        if interval1_fixed[1] - interval1_fixed[0] < 50: continue
        for interval2 in intervals:
            if interval1 == interval2: continue
            interval2_fixed = (min(interval2), max(interval2))
            if interval1_fixed[0] >= interval2_fixed[0]: continue
            if interval2_fixed[1] - interval2_fixed[0] < 50: continue

            overlap_interval = (
                max(interval1_fixed[0], interval2_fixed[0]), min(interval1_fixed[1], interval2_fixed[1]))
            overlap_interval = (min(overlap_interval), max(overlap_interval))

            if overlap_interval[1] - overlap_interval[0] <= 20 and overlap_interval[0] >= window and overlap_interval[
                1] <= window + reference_length:
                interval_pairs.append((interval1_fixed, interval2_fixed))
                overlap_intervals.add(overlap_interval)

    interval_pairs = list(set(interval_pairs))
    logger.info("Paired %s as %s" % (str(intervals), str(interval_pairs)))
    logger.info("Overlap intervals %s" % (str(overlap_intervals)))

    if not overlap_intervals:
        return 0, -1, -1

    low = min([interval[0] for interval in overlap_intervals])
    high = max([interval[1] for interval in overlap_intervals])

    logger.info("low = %d, high = %d overlap_intervals = %s" % (low, high, str(overlap_intervals)))
    if high - low <= 20:
        logger.info("Pairing success! The overlap intervals are close enough.")
        return len(overlap_intervals), low, high
    return len(overlap_intervals), -1, -1


def pair_intervals(intervals, reference_length, window=50):
    interval_pairs = []
    overlap_intervals = set()
    for interval1 in intervals:
        interval1_fixed = (min(interval1), max(interval1))
        if interval1_fixed[1] - interval1_fixed[0] < 50: continue
        for interval2 in intervals:
            if interval1 == interval2: continue
            interval2_fixed = (min(interval2), max(interval2))
            if interval1_fixed[0] >= interval2_fixed[0]: continue
            if interval2_fixed[1] - interval2_fixed[0] < 50: continue

            overlap_interval = (
                max(interval1_fixed[0], interval2_fixed[0]), min(interval1_fixed[1], interval2_fixed[1]))
            overlap_interval = (min(overlap_interval), max(overlap_interval))

            if overlap_interval[1] - overlap_interval[0] <= 20 and overlap_interval[0] >= window and overlap_interval[
                1] <= window + reference_length:
                interval_pairs.append((interval1_fixed, interval2_fixed))
                overlap_intervals.add(overlap_interval)

    interval_pairs = list(set(interval_pairs))
    logger.info("Paired %s as %s" % (str(intervals), str(interval_pairs)))
    logger.info("Overlap intervals %s" % (str(overlap_intervals)))

    if not overlap_intervals: return 0, -1, -1

    low = min([interval[0] for interval in overlap_intervals])
    high = max([interval[1] for interval in overlap_intervals])

    logger.info("low = %d, high = %d overlap_intervals = %s" % (low, high, str(overlap_intervals)))
    if high - low <= 20:
        logger.info("Pairing success! The overlap intervals are close enough.")
        return len(overlap_intervals), low, high
    return len(overlap_intervals), -1, -1


def get_good_assemblies(assemblies, min_len=25, window=50):
    # Pair-wise analysis of assemblies for long insertions

    good_assemblies = []
    for a1 in assemblies:
        if a1.used: continue
        interval1 = a1.start1_end1s[0] if a1.has_only_long_left_flank(min_len) else a1.start1_end1s[1]
        interval1 = [min(interval1), max(interval1)]
        for a2 in assemblies:
            if a2.used: continue
            if a1.tigra_contig.raw_name == a2.tigra_contig.raw_name:
                continue
            interval2 = a2.start1_end1s[0] if a2.has_only_long_left_flank(min_len) else a2.start1_end1s[1]
            interval2 = [min(interval2), max(interval2)]

            overlap_interval = [max(interval1[0], interval2[0]), min(interval1[1], interval2[1])]
            overlap_interval = [min(overlap_interval), max(overlap_interval)]
            if overlap_interval[1] - overlap_interval[0] <= 20:
                good_assemblies += [a1, a2]
                a1.used = True
                a2.used = True

                a1.start = a1.tigra_contig.sv_region.pos1 + overlap_interval[0] - 1 - window
                a1.end = a1.tigra_contig.sv_region.pos1 + overlap_interval[1] - window

                a2.start = a1.start
                a2.end = a1.end
                break

    return good_assemblies


def get_nonduplicate_assemblies(assemblies):
    for a1 in assemblies:
        if a1.duplicate: continue
        for a2 in assemblies:
            if a1.tigra_contig.raw_name == a2.tigra_contig.raw_name: continue
            if a2.duplicate: continue

            if a1.start1_end1s == a2.start1_end1s:
                a2.duplicate = True

    return [assembly for assembly in assemblies if not assembly.duplicate]


def get_reference_intervals(age_records, start=0, min_interval_len=100):
    intervals = []
    for age_record in age_records:
        intervals += map(lambda x: (min(x) + start - 1, max(x) + start - 1),
                         [interval for interval in age_record.start1_end1s if
                          abs(interval[0] - interval[1]) >= min_interval_len])

    return intervals


def process_age_records(age_records, sv_type="INS", ins_min_unaligned=10, min_interval_len=200, pad=500,
                        min_deletion_len=30):
    func_logger = logging.getLogger("%s-%s" % (process_age_records.__name__, multiprocessing.current_process()))

    good_age_records = age_records
    if sv_type == "INS":
        good_age_records = [age_record for age_record in good_age_records if
                            not age_record.almost_all_bases_aligned(ins_min_unaligned)]
        good_age_records = [age_record for age_record in good_age_records if not age_record.is_reference()]
    elif sv_type == "DEL":
        good_age_records = [age_record for age_record in good_age_records if
                            len(age_record.start1_end1s) == 2 and min(age_record.ref_flanking_regions) >= 50]
        good_age_records = [age_record for age_record in good_age_records if
                            abs(age_record.start1_end1s[0][1] - age_record.start1_end1s[1][0]) >= min_deletion_len]
        good_age_records = [age_record for age_record in good_age_records if
                            float(age_record.score) / sum(age_record.ref_flanking_regions) >= 0.7]
    elif sv_type == "INV":
        pass

    # Add some features to an info dict
    info = defaultdict(float)
    info["BA_NUM_GOOD_REC"] = len(good_age_records) if good_age_records else 0
    for rec in good_age_records:
        info["BA_FLANK_PERCENT"] = max(info["BA_FLANK_PERCENT"], rec.flank_percent)
        info["BA_NFRAGS"] = max(info["BA_NFRAGS"], rec.nfrags)
        info["BA_NUM_ALT"] = max(info["BA_NUM_ALT"], rec.n_alt)
        info["BA_PERCENT_MATCH"] = max(info["BA_PERCENT_MATCH"], rec.percent)

    if not good_age_records:
        func_logger.warning("No good records found for getting breakpoints")
        return [], dict(info)
    else:
        func_logger.info("Found %d good records for getting breakpoints" % (len(good_age_records)))
        func_logger.info("Good records")
        for age_record in good_age_records:
            func_logger.info(str(age_record))

    sv_region = good_age_records[0].contig.sv_region
    if sv_type == "DEL":
        breakpoints = get_deletion_breakpoints(good_age_records, start=sv_region.pos1 - pad)
    elif sv_type == "INS":
        reference_intervals = get_reference_intervals(good_age_records, start=1, min_interval_len=min_interval_len)

        func_logger.info("Gathered reference intervals as %s" % (str(reference_intervals)))
        breakpoints = get_insertion_breakpoints(good_age_records, reference_intervals, start=sv_region.pos1 - pad)
    else:
        return [], dict(info)

    func_logger.info("Detected breakpoints as %s" % (str(breakpoints)))

    # Add a few more features related to the breakpoints computed
    info["BA_NUM_BP"] = len(breakpoints)

    if sv_type == "DEL":
        if len(breakpoints) == 2:
            func_logger.info("True deletion interval %s" % (str(breakpoints)))
        else:
            func_logger.info("False deletion interval %s" % (str(breakpoints)))
            return [], dict(info)
    elif sv_type == "INS":
        if len(breakpoints) == 1:
            if sv_region.pos2 - sv_region.pos1 <= 20:
                info["BA_BP_SCORE"] = abs(breakpoints[0][0] - sv_region.pos1)
                if abs(breakpoints[0][0] - sv_region.pos1) > 20:
                    return [], dict(info)
            else:
                diff1 = breakpoints[0][0] - sv_region.pos1
                diff2 = sv_region.pos2 - breakpoints[0][0]
                info["BA_BP_SCORE"] = min(abs(diff1 - pad), abs(diff2 - pad))
                if not (pad - 10 <= diff1 <= pad + 10 or pad - 10 <= diff2 <= pad + 10):
                    return [], dict(info)
            func_logger.info("True insertion interval %s" % (str(breakpoints)))
        else:
            return [], dict(info)

    return breakpoints, dict(info)
