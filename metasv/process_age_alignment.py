from __future__ import print_function
import logging
import multiprocessing
from collections import defaultdict
import pybedtools

logger = logging.getLogger(__name__)

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
        # Check if the breakpoint is within window distance of a validated breakpoint
        if min([window + 1] + [abs(b[0] - breakpoint) for b in breakpoints]) <= window:
            continue
        func_logger.info("\tExamining potential breakpoint %d for support" % breakpoint)
        left_support = [interval[0] for interval in intervals if abs(interval[0] - breakpoint) <= window]
        right_support = [interval[1] for interval in intervals if abs(interval[1] - breakpoint) <= window]
        counter_examples = [age_record for age_record in age_records if age_record.has_long_ref_flanks() and (
            age_record.has_ref_deletion(window) or age_record.has_insertion(min_diff=1,
                                                                            max_diff=49)) and age_record.breakpoint_match(
            breakpoint, window)]
        if counter_examples:
            counter_example_ends = [age_record.start1_end1s for age_record in counter_examples]
            func_logger.info("\t\tSkipping breakpoint %d due to %s" % (breakpoint, str(counter_example_ends)))
            continue

        if left_support:
            func_logger.info("\t\tLeft support %s" % (str(left_support)))
        if right_support:
            func_logger.info("\t\tRight support %s" % (str(right_support)))

        if (left_support and right_support) and min(
                        [window + 1] + [abs(b[0] - breakpoint) for b in breakpoints]) > window:
            both_support = [age_record for age_record in age_records if
                            age_record.has_insertion(min_diff=50, max_diff=1000000000) and age_record.breakpoint_match(
                                breakpoint, window)]
            if both_support:
                func_logger.info("\t\tboth_support = %s" % (str(both_support)))
                func_logger.info("\t\tinsertion lengths = %s" % (
                str([age_record.insertion_length() for age_record in both_support])))
            insertion_length = max([0] + [age_record.insertion_length() for age_record in both_support])
            func_logger.info("\t\tInsertion length = %d" % insertion_length)
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



def get_inversion_breakpoints(age_records, window=20, min_endpoint_dist=0, start=0,pad=500,dist_to_expected_bp=400,min_interval_len_inv=100):
    func_logger = logging.getLogger("%s-%s" % (get_deletion_breakpoints.__name__, multiprocessing.current_process()))

    potential_breakpoints = []
    for age_record in age_records:
    
        polarities=[abs(age_record.polarities1[i]-age_record.polarities2[i]) for i in range(age_record.nfrags)]               
        good_intervals=[i for i in range(age_record.nfrags) if abs(age_record.start1_end1s[i][1]-age_record.start1_end1s[i][0]) > min_interval_len_inv and abs(age_record.start2_end2s[i][1]-age_record.start2_end2s[i][0]) > min_interval_len_inv]
        good_intervals=[i for i in good_intervals if abs(age_record.start1_end1s[i][1]-age_record.start1_end1s[i][0]) <= max(age_record.inputs[0].length-2*(pad-dist_to_expected_bp),pad+dist_to_expected_bp)]
        func_logger.info('Good intervals: %s'%str(good_intervals))
        if len(good_intervals)<2:
            func_logger.info('Not enough good interval for this age record: %s'%str(age_record))
            continue
        candidate_inv_intervals=[]
        far_from_endpoint_intervals=[]
        long_inversion=False
        for i in good_intervals:
            interval=age_record.start2_end2s[i]
            if min([interval[0],abs(interval[0]-age_record.inputs[1].length),
                    interval[1],abs(interval[1]-age_record.inputs[1].length)]) >= min_endpoint_dist:
                far_from_endpoint_intervals.append(i)
                if (abs(min(age_record.start1_end1s[i])-pad)<dist_to_expected_bp) and (abs(max(age_record.start1_end1s[i])-(age_record.inputs[0].length-pad))<dist_to_expected_bp): 
                    candidate_inv_intervals.append(i)
        
        if len(candidate_inv_intervals)>1:
            dist_to_exp_bps=map(lambda x: abs(min(age_record.start1_end1s[x])-pad)+abs(max(age_record.start1_end1s[x])-(age_record.inputs[0].length-pad)),candidate_inv_intervals)
            candidate_inv_intervals=[min(enumerate(dist_to_exp_bps),key=lambda x:x[1])[0]]
        
        #Potentially long inversion
        if not candidate_inv_intervals:
            for i in far_from_endpoint_intervals:
                interval_ends=sorted(age_record.start1_end1s[i])
                if (abs(min(interval_ends)-pad)<dist_to_expected_bp):
                    #if abs(interval_ends[0]-pad)<abs(interval_ends[1]-pad):
                    candidate_inv_intervals.append(i)
                    func_logger.info('Potentially long-inversion interval: %s'%i)
                    long_inversion=True
                elif (abs(max(interval_ends)-(age_record.inputs[0].length-pad))<dist_to_expected_bp): 
                    #if abs(interval_ends[0]-(age_record.inputs[0].length-pad))>abs(interval_ends[1]-(age_record.inputs[0].length-pad)):
                    candidate_inv_intervals.append(i)
                    func_logger.info('Potentially long-inversion interval: %s'%i)
                    long_inversion=True
                    
                        
            
            
        if not candidate_inv_intervals:
            func_logger.info('Not candidate inversion interval found for this age record: %s'%str(age_record))
            continue

        ref_polarities=[(i,polarities[i],abs(age_record.start1_end1s[i][0]-age_record.start1_end1s[i][1])) for i in good_intervals if i not in candidate_inv_intervals]

        func_logger.info('age_record: %s'%str(age_record))
        func_logger.info('candidate_inv_intervals: %s'%str(candidate_inv_intervals))
        func_logger.info('ref_polarities: %s'%str(ref_polarities))

        if not ref_polarities:
            func_logger.info('Cannot find reference polarity for this age record: %s'%str(age_record))
            continue
        ref_interval,ref_polarity,ref_length=max(ref_polarities,key=lambda x:x[2])
        
        candidate_inv_intervals=map(lambda i: (i,abs((age_record.start1_end1s[i][0]-age_record.start1_end1s[i][1]))),filter(lambda i: polarities[i]!=ref_polarity,candidate_inv_intervals))
        if candidate_inv_intervals:
            inv_interval=max(candidate_inv_intervals,key=lambda x:x[1])[0]
            
            s_inv=sorted(age_record.start1_end1s[inv_interval])
            s_ref=sorted(age_record.start1_end1s[ref_interval])
            if (s_ref[0]-s_inv[0])*(s_ref[1]-s_inv[1])<=0:
                func_logger.info('Bad intervals: %s'%str(age_record))
                continue
                
            
            if not long_inversion:
            
                interval=age_record.start2_end2s[inv_interval]            
                if min([interval[0],abs(interval[0]-age_record.inputs[1].length),
                        interval[1],abs(interval[1]-age_record.inputs[1].length)]) < 10:
                    func_logger.info('Inverted interval end points are too close to borders in Seq2: %s'%str(age_record))
                    continue
                        
    
            
                if ((s_ref[1]>s_inv[1]) and (abs(s_inv[1]-s_ref[0])>10)) or ((s_ref[0]<s_inv[0]) and (abs(s_inv[0]-s_ref[1])>10)): 
                    func_logger.info('Bad bp2  (does not match in seq1): %s'%str(age_record))
                    continue

                if ((s_ref[1]>s_inv[1]) and ((s_inv[1]-s_ref[0])>10)) or ((s_ref[0]<s_inv[0]) and ((s_ref[1]-s_inv[0])>10)): 
                    func_logger.info('Large overlaps between intervals in seq1: %s'%str(age_record))
                    continue



                if (s_ref[1]>s_inv[1]):
                    mbp_L_seq2=age_record.start2_end2s[inv_interval][filter(lambda x:age_record.start1_end1s[inv_interval][x]==s_inv[0],[0,1])[0]]
                    mbp_R_seq2=age_record.start2_end2s[ref_interval][filter(lambda x:age_record.start1_end1s[ref_interval][x]==s_ref[0],[0,1])[0]]
                else:
                    mbp_L_seq2=age_record.start2_end2s[ref_interval][filter(lambda x:age_record.start1_end1s[ref_interval][x]==s_ref[1],[0,1])[0]]
                    mbp_R_seq2=age_record.start2_end2s[inv_interval][filter(lambda x:age_record.start1_end1s[inv_interval][x]==s_inv[1],[0,1])[0]]
                
                if abs(mbp_L_seq2-mbp_R_seq2)>10:
                    func_logger.info('BPs do not match in seq2: %s'%str(age_record))
                    continue

                                
                potential_breakpoints += age_record.start1_end1s[inv_interval]
            else:
                if (abs(min(age_record.start1_end1s[inv_interval])-pad)<dist_to_expected_bp):
                    bp1=min(age_record.start1_end1s[inv_interval])
                    bp1_otherend=max(age_record.start1_end1s[inv_interval])
                else:
                    bp1=max(age_record.start1_end1s[inv_interval])
                    bp1_otherend=min(age_record.start1_end1s[inv_interval])
                bp1_seq2=age_record.start2_end2s[inv_interval][filter(lambda x:age_record.start1_end1s[inv_interval][x]==bp1,[0,1])[0]]
                
                if abs(age_record.start1_end1s[ref_interval][0]-bp1)<=abs(age_record.start1_end1s[ref_interval][1]-bp1):
                    bp2=age_record.start1_end1s[ref_interval][0]
                    bp2_otherend=age_record.start1_end1s[ref_interval][1]
                else:
                    bp2=age_record.start1_end1s[ref_interval][1]
                    bp2_otherend=age_record.start1_end1s[ref_interval][0]
                
                bp2_seq2=age_record.start2_end2s[ref_interval][filter(lambda x:age_record.start1_end1s[ref_interval][x]==bp2,[0,1])[0]]

                if not ((abs(bp2-pad)<dist_to_expected_bp) or (abs(bp2-(age_record.inputs[0].length-pad))<dist_to_expected_bp)):
                    func_logger.info('Bad bp2 in normal fragment: %s'%str(age_record))
                    continue


                if abs(bp1_seq2-bp2_seq2)>10:
                    func_logger.info('BPs do not match in seq2: %s'%str(age_record))
                    continue

                if (bp1_otherend-bp1)*(bp2_otherend-bp2)<0:
                    func_logger.info('BPs directions do not match in seq1: %s'%str(age_record))
                    continue
                    
                
                potential_breakpoints += [bp1,bp2]
            
    
    potential_breakpoints=sorted(potential_breakpoints)
    breakpoints = []
    for breakpoint in potential_breakpoints:
        if min([window + 1] + [abs(b - breakpoint) for b in breakpoints]) >= window:
            breakpoints.append(breakpoint)

    func_logger.info("Gathered breakpoints as %s" % (str(breakpoints)))

    return [start + breakpoint for breakpoint in breakpoints]


def get_reference_intervals(age_records, start=0, min_interval_len=100):
    intervals = []
    for age_record in age_records:
        intervals += map(lambda x: (min(x) + start - 1, max(x) + start - 1),
                         [interval for interval in age_record.start1_end1s if
                          abs(interval[0] - interval[1]) >= min_interval_len])

    return intervals


def process_age_records(age_records, sv_type="INS", ins_min_unaligned=10, min_interval_len=200, pad=500,
                        min_deletion_len=30,min_interval_len_inv=100):
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
        good_age_records = [age_record for age_record in good_age_records if
                            len(age_record.start1_end1s) >= 2 and min(map(lambda x:abs(x[1]-x[0]),age_record.start1_end1s)) >= min_interval_len_inv]
    else:
        pass
    # Add some features to an info dict
    info = defaultdict(float)
    info["BA_NUM_GOOD_REC"] = len(good_age_records)
    if not good_age_records:
        func_logger.warning("No good records found for getting breakpoints")
        return [], dict(info)

    for rec in good_age_records:
        info["BA_FLANK_PERCENT"] = max(info["BA_FLANK_PERCENT"], rec.flank_percent)
        info["BA_NFRAGS"] = max(info["BA_NFRAGS"], rec.nfrags)
        info["BA_NUM_ALT"] = max(info["BA_NUM_ALT"], rec.n_alt)
        info["BA_PERCENT_MATCH"] = max(info["BA_PERCENT_MATCH"], rec.percent)

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
    elif sv_type == "INV":
        breakpoints = get_inversion_breakpoints(good_age_records, start=sv_region.pos1 - pad ,pad=pad, min_interval_len_inv=min_interval_len_inv)
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
    elif sv_type == "INV":
        if len(breakpoints) == 2:
            func_logger.info("True inversion interval %s" % (str(breakpoints)))
        else:
            func_logger.info("False inversion interval %s" % (str(breakpoints)))
            return [], dict(info)

    return breakpoints, dict(info)
