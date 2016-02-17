from __future__ import print_function
import logging
import multiprocessing
from collections import defaultdict
import pybedtools

from defaults import MIN_INV_SUBALIGN_LENGTH, MIN_DEL_SUBALIGN_LENGTH,AGE_WINDOW_SIZE

logger = logging.getLogger(__name__)

def get_insertion_breakpoints(age_records, intervals, expected_bp_pos, window=AGE_WINDOW_SIZE, start=0, dist_to_expected_bp=50):
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
            age_record.has_ref_deletion(window) or age_record.has_insertion(min_diff=20,
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
            insertion_sequence = both_support[0].get_insertion_sequence if both_support else ""
            func_logger.info("\t\tInsertion length = %d %s" % (insertion_length, insertion_sequence))
            breakpoints.append((breakpoint, insertion_length, insertion_sequence))

    func_logger.info("Nonfiltered breakpoints as %s" % (str(breakpoints)))

    if len(breakpoints)>1:
        breakpoints=filter(lambda x: min(abs(x[0]-expected_bp_pos[0]),abs(expected_bp_pos[1]-x[0]))<dist_to_expected_bp,breakpoints)

    func_logger.info("Gathered breakpoints as %s" % (str(breakpoints)))

    return [(start + b[0], b[1], b[2]) for b in breakpoints]


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

def check_closeness_to_bp(pos,pad,dist_to_expected_bp,LR_bp,seq_length=0):
    if LR_bp == 'L':
        return abs(pos-pad)<dist_to_expected_bp
    else:
        return abs(pos-(seq_length-pad))<dist_to_expected_bp
    


def get_inversion_breakpoints(age_records, window=20, min_endpoint_dist=10, start=0, pad=500, dist_to_expected_bp=400, min_inv_subalign_len=MIN_INV_SUBALIGN_LENGTH):
    func_logger = logging.getLogger("%s-%s" % (get_deletion_breakpoints.__name__, multiprocessing.current_process()))

    potential_breakpoints = []
    for age_record in age_records:
    
        polarities=[abs(age_record.polarities1[i]-age_record.polarities2[i]) for i in range(age_record.nfrags)]               
        good_intervals=[i for i in range(age_record.nfrags) if abs(age_record.start1_end1s[i][1]-age_record.start1_end1s[i][0]) > min_inv_subalign_len and abs(age_record.start2_end2s[i][1]-age_record.start2_end2s[i][0]) > min_inv_subalign_len]
        good_intervals=[i for i in good_intervals if abs(age_record.start1_end1s[i][1]-age_record.start1_end1s[i][0]) <= max(age_record.inputs[0].length-2*(pad-dist_to_expected_bp),pad+dist_to_expected_bp)]
        func_logger.info('Good intervals: %s'%str(good_intervals))
        if len(good_intervals)<2:
            func_logger.info('Not enough good interval for this age record: %s'%str(age_record))
            continue
        candidate_inv_intervals=[]
        inv_interval=-1
        long_inversion=False
        
        left_end_near_l_bp=filter(lambda x: check_closeness_to_bp(min(age_record.start1_end1s[x]),pad,dist_to_expected_bp,"L"), good_intervals)
        right_end_near_r_bp=filter(lambda x: check_closeness_to_bp(max(age_record.start1_end1s[x]),pad,dist_to_expected_bp,"R",age_record.inputs[0].length), good_intervals)
        
        right_end_near_l_bp=filter(lambda x: check_closeness_to_bp(max(age_record.start1_end1s[x]),pad,dist_to_expected_bp,"L"), good_intervals)
        left_end_near_r_bp=filter(lambda x: check_closeness_to_bp(min(age_record.start1_end1s[x]),pad,dist_to_expected_bp,"R",age_record.inputs[0].length), good_intervals)

        candidate_inv_intervals=list(set(left_end_near_l_bp)&set(right_end_near_r_bp))
        candidate_norm_intervals=list(set(left_end_near_r_bp)|set(right_end_near_l_bp))
        
        if len(candidate_inv_intervals)>1 and len(candidate_norm_intervals)<=1: 
            candidate_inv_intervals=list(set(candidate_inv_intervals)-set(candidate_norm_intervals))
        
        if len(candidate_inv_intervals)>1:
            dist_to_exp_bps=map(lambda x: abs(min(age_record.start1_end1s[x])-pad)+abs(max(age_record.start1_end1s[x])-(age_record.inputs[0].length-pad)),candidate_inv_intervals)
            inv_interval=min(enumerate(dist_to_exp_bps),key=lambda x:x[1])[0]
        elif len(candidate_inv_intervals)==1 :
            inv_interval=candidate_inv_intervals[0]
        
        
        if inv_interval==-1:
            #Potentially long inversion
            candidate_inv_intervals=[i for i in left_end_near_l_bp if ((set(candidate_norm_intervals)&set(left_end_near_r_bp))-set([i]))] + \
                                    [i for i in right_end_near_r_bp if ((set(candidate_norm_intervals)&set(right_end_near_l_bp))-set([i]))]
            
            if len(candidate_inv_intervals)>1:
                candidate_inv_intervals=[i for i in set(candidate_inv_intervals)&set(left_end_near_l_bp) if (pad< (sum(age_record.start1_end1s[i])/2.0))] + \
                                        [i for i in set(candidate_inv_intervals)&set(right_end_near_r_bp) if ((age_record.inputs[0].length-pad) > (sum(age_record.start1_end1s[i])/2.0))]
            
            if candidate_inv_intervals:
                func_logger.info('Potentially long-inversion interval: %s'%candidate_inv_intervals)
                long_inversion=True
                if len(candidate_inv_intervals)>1:
                    dist_to_exp_bps=map(lambda x: abs(min(age_record.start1_end1s[x])-pad) if i in left_end_near_l_bp else abs(max(age_record.start1_end1s[x])-(age_record.inputs[0].length-pad)),candidate_inv_intervals)
                    inv_interval=min(enumerate(dist_to_exp_bps),key=lambda x:x[1])[0]
                else:
                    inv_interval=candidate_inv_intervals[0]
        elif age_record.inputs[0].length > ((2*pad+min_inv_subalign_len)):
            long_inversion=True

        if inv_interval==-1:
            func_logger.info('Not candidate inversion interval found for this age record: %s'%str(age_record))
            continue

        func_logger.info('age_record: %s'%str(age_record))
        func_logger.info('inverted interval: %s'%str(inv_interval))

        candidate_norm_intervals=filter(lambda x: polarities[x]!=polarities[inv_interval], set(candidate_norm_intervals)-set([inv_interval]))
        if long_inversion and (inv_interval not in set(left_end_near_l_bp) & set(right_end_near_r_bp)) :
            candidate_norm_intervals=list(set(candidate_norm_intervals)&set(left_end_near_r_bp if (inv_interval in left_end_near_l_bp) else right_end_near_l_bp))
    
        if not candidate_norm_intervals:
            func_logger.info('Cannot find the normal interval for this age record: %s'%str(age_record))
            continue
        
        if len(candidate_norm_intervals)>1:
            candidate_norm_intervals=map(lambda x: (x,abs(age_record.start1_end1s[x][0]-age_record.start1_end1s[x][1])),set(candidate_norm_intervals))
            norm_interval,norm_length=max(candidate_norm_intervals,key=lambda x:x[2])
        else:
            norm_interval=candidate_norm_intervals[0]

        func_logger.info('norm_interval: %s'%str(norm_interval))

        s_inv=sorted(age_record.start1_end1s[inv_interval])
        s_norm=sorted(age_record.start1_end1s[norm_interval])
        if (s_norm[0]-s_inv[0])*(s_norm[1]-s_inv[1])<=0:
            func_logger.info('Bad intervals (one fully covers the other): %s'%str(age_record))
            continue
        
        if not long_inversion:        
            interval=age_record.start2_end2s[inv_interval]            
            if min([interval[0],abs(interval[0]-age_record.inputs[1].length),
                    interval[1],abs(interval[1]-age_record.inputs[1].length)]) < min_endpoint_dist:
                func_logger.info('Inverted interval end points are too close to borders in Seq2: %s'%str(age_record))
                continue
            if (((s_norm[1]>s_inv[1]) and ((s_inv[1]-s_norm[0])>10)) or ((s_norm[0]<s_inv[0]) and ((s_norm[1]-s_inv[0])>10))):
                func_logger.info('Bad middle bp in seq1 (covers>10): %s'%str(age_record))
                continue    
            if (((s_norm[1]>s_inv[1]) and ((s_norm[0]-s_inv[1])>50)) or ((s_norm[0]<s_inv[0]) and ((s_inv[0]-s_norm[1])>50))):
                func_logger.info('Bad middle bp in seq1 (apart>50): %s'%str(age_record))
                continue

        bp_idx = 0 if (s_norm[1]>s_inv[1]) else 1
        bp1=s_inv[bp_idx]
        bp2=s_norm[bp_idx]
        bp1_seq2=age_record.start2_end2s[inv_interval][filter(lambda x:age_record.start1_end1s[inv_interval][x]==bp1,[0,1])[0]]
        bp2_seq2=age_record.start2_end2s[norm_interval][filter(lambda x:age_record.start1_end1s[norm_interval][x]==bp2,[0,1])[0]]

        if abs(bp1_seq2-bp2_seq2)>10:
            func_logger.info('BPs do not match in seq2: %s'%str(age_record))
            continue
        potential_breakpoints += [bp1,bp2]
            
    potential_breakpoints=sorted(potential_breakpoints)
    breakpoints = []
    for breakpoint in potential_breakpoints:
        if min([window + 1] + [abs(b - breakpoint) for b in breakpoints]) >= window:
            breakpoints.append(breakpoint)

    func_logger.info("Gathered breakpoints as %s" % (str(breakpoints)))

    return [start + breakpoint for breakpoint in breakpoints]


def get_duplication_breakpoints(age_records, window=20, max_endpoint_dist=10, start=0, pad=500, dist_to_expected_bp=400):
    func_logger = logging.getLogger("%s-%s" % (get_deletion_breakpoints.__name__, multiprocessing.current_process()))

    potential_breakpoints = []
    for age_record in age_records:
        left_end_near_l_bp=filter(lambda x: check_closeness_to_bp(min(age_record.start1_end1s[x]),pad,dist_to_expected_bp,"L"), [0,1])
        right_end_near_r_bp=filter(lambda x: check_closeness_to_bp(max(age_record.start1_end1s[x]),pad,dist_to_expected_bp,"R",age_record.inputs[0].length), [0,1])
        if (not left_end_near_l_bp) or (not right_end_near_r_bp):
            func_logger.info('Not close to expected BPs: %s'%str(age_record))        
            continue
        if len(left_end_near_l_bp)==2 and len(right_end_near_r_bp)==1:
            left_end_near_l_bp = list(set(left_end_near_l_bp)-set(right_end_near_r_bp))
        elif len(left_end_near_l_bp)==1 and len(right_end_near_r_bp)==2:
            right_end_near_r_bp = list(set(right_end_near_r_bp)-set(left_end_near_l_bp))
        elif len(left_end_near_l_bp)==2 and len(right_end_near_r_bp)==2:
            dist_to_exp_l_bp=map(lambda x: abs(min(age_record.start1_end1s[x])-pad),[0,1])
            dist_to_exp_r_bp=map(lambda x: abs(max(age_record.start1_end1s[x])-(age_record.inputs[0].length-pad)),[0,1])
            left_end_near_l_bp, right_end_near_r_bp = [[0],[1]]  if (dist_to_exp_l_bp[0]+dist_to_exp_r_bp[1]) < (dist_to_exp_l_bp[1]+dist_to_exp_r_bp[0]) else [[1],[0]]
        
        l_interval = left_end_near_l_bp[0]
        r_interval = right_end_near_r_bp[0]
        
        bp_idx_l = 0 if age_record.start1_end1s[l_interval][0]<age_record.start1_end1s[l_interval][1] else 1
        bp_idx_r = 1 if age_record.start1_end1s[r_interval][0]<age_record.start1_end1s[r_interval][1] else 0
        if abs(age_record.start2_end2s[l_interval][bp_idx_l]-age_record.start2_end2s[r_interval][bp_idx_r]) > 10:
            func_logger.info('BPs do not match in seq2: %s'%str(age_record))
            continue
        
        end_l_seq2 = age_record.start2_end2s[l_interval][1-bp_idx_l]
        end_r_seq2 = age_record.start2_end2s[r_interval][1-bp_idx_r]              
        if max(min(end_r_seq2,end_l_seq2),
               min(end_l_seq2-age_record.inputs[1].length,
                   end_r_seq2-age_record.inputs[1].length)) > max_endpoint_dist:
            func_logger.info('End points are too close to borders in Seq2: %s'%str(age_record))
            continue
        potential_breakpoints += [age_record.start1_end1s[l_interval][bp_idx_l],age_record.start1_end1s[r_interval][bp_idx_r]]

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
                        min_deletion_len=30, min_del_subalign_len=MIN_DEL_SUBALIGN_LENGTH, 
                        min_inv_subalign_len=MIN_INV_SUBALIGN_LENGTH, dist_to_expected_bp=400,
                        age_window=AGE_WINDOW_SIZE, pad_ins=0):
    func_logger = logging.getLogger("%s-%s" % (process_age_records.__name__, multiprocessing.current_process()))

    good_age_records = age_records
    if sv_type == "INS":
        good_age_records = [age_record for age_record in good_age_records if
                            not age_record.almost_all_bases_aligned(ins_min_unaligned)]
        good_age_records = [age_record for age_record in good_age_records if not age_record.is_reference()]
    elif sv_type == "DEL":
        good_age_records = [age_record for age_record in good_age_records if
                            len(age_record.start1_end1s) == 2 and min(age_record.ref_flanking_regions) >= min_del_subalign_len]
        good_age_records = [age_record for age_record in good_age_records if
                            abs(age_record.start1_end1s[0][1] - age_record.start1_end1s[1][0]) >= min_deletion_len]
        good_age_records = [age_record for age_record in good_age_records if
                            float(age_record.score) / sum(age_record.ref_flanking_regions) >= 0.7]

        good_age_records = [age_record for age_record in good_age_records if
                            abs(age_record.start2_end2s[0][1] - age_record.start2_end2s[1][0]) <= 50]
                                
        good_age_records = [age_record for age_record in good_age_records if 
                            check_closeness_to_bp(min(age_record.start1_end1s[0][1],
                                                  age_record.start1_end1s[1][0]),
                                                  pad,dist_to_expected_bp,"L") and 
                            check_closeness_to_bp(max(age_record.start1_end1s[0][1],
                                                  age_record.start1_end1s[1][0]),
                                                  pad,dist_to_expected_bp,"R",
                                                  age_record.inputs[0].length)]
    elif sv_type == "INV":
        good_age_records = [age_record for age_record in good_age_records if
                            len(age_record.start1_end1s) >= 2 and min(map(lambda x:abs(x[1]-x[0]),age_record.start1_end1s)) >= min_inv_subalign_len]
    elif sv_type == "DUP":
        good_age_records = [age_record for age_record in good_age_records if
                            len(age_record.start1_end1s) == 2 and min(age_record.ref_flanking_regions) >= 100]
    else:
        pass
    # Add some features to an info dict
    info = defaultdict(int)
    info["BA_NUM_GOOD_REC"] = len(good_age_records)
    if not good_age_records:
        func_logger.warning("No good records found for getting breakpoints")
        return [], dict(info)

    for rec in good_age_records:
        info["BA_FLANK_PERCENT"] = int(max(info["BA_FLANK_PERCENT"], rec.flank_percent))
        info["BA_NFRAGS"] = int(max(info["BA_NFRAGS"], rec.nfrags))
        info["BA_NUM_ALT"] = int(max(info["BA_NUM_ALT"], rec.n_alt))
        info["BA_PERCENT_MATCH"] = int(max(info["BA_PERCENT_MATCH"], rec.percent))

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
        breakpoints = get_insertion_breakpoints(good_age_records, reference_intervals, 
                                                expected_bp_pos=[pad+pad_ins,max((sv_region.pos2-sv_region.pos1)-pad_ins+pad,0)],
                                                window=age_window,
                                                start=sv_region.pos1 - pad)
    elif sv_type == "INV":
        breakpoints = get_inversion_breakpoints(good_age_records, start=sv_region.pos1 - pad ,pad=pad, min_inv_subalign_len=min_inv_subalign_len, dist_to_expected_bp=dist_to_expected_bp)
    elif sv_type == "DUP":
        breakpoints = get_duplication_breakpoints(good_age_records, start=sv_region.pos1 - pad ,pad=pad, dist_to_expected_bp=dist_to_expected_bp)
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
#             if sv_region.pos2 - sv_region.pos1 <= 20:
#                 info["BA_BP_SCORE"] = abs(breakpoints[0][0] - sv_region.pos1)
#                 if abs(breakpoints[0][0] - sv_region.pos1) > 20:
#                     return [], dict(info)
#             else:
            diff1 = breakpoints[0][0] - (sv_region.pos1+pad_ins)
            diff2 = (sv_region.pos2-pad_ins) - breakpoints[0][0]
            info["BA_BP_SCORE"] = min(abs(diff1), abs(diff2 ))
            if not  min(abs(diff1), abs(diff2 )) <=  100 :
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
    elif sv_type == "DUP":
        if len(breakpoints) == 2:
            func_logger.info("True duplication interval %s" % (str(breakpoints)))
        else:
            func_logger.info("False duplication interval %s" % (str(breakpoints)))
            return [], dict(info)
    return breakpoints, dict(info)
