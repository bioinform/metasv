import argparse
import logging
import multiprocessing
import time
from functools import partial, update_wrapper
from defaults import EXTRACTION_MAX_READ_PAIRS, EXTRACTION_MAX_NM, EXTRACTION_MAX_INTERVAL_TRUNCATION, EXTRACTION_TRUNCATION_PAD

import pysam

compl_table = [chr(i) for i in xrange(256)]
compl_table[ord('A')] = 'T'
compl_table[ord('C')] = 'G'
compl_table[ord('G')] = 'C'
compl_table[ord('T')] = 'A'


def compl(seq):
    return "".join([compl_table[ord(i)] for i in seq])


def get_sequence_quality(aln):
    if not aln.is_reverse:
        return aln.seq.upper(), aln.qual

    return compl(aln.seq.upper())[::-1], aln.qual[::-1]


def write_read(fd, aln):
    end_id = 1 if aln.is_read1 else 2

    sequence, quality = get_sequence_quality(aln)
    fd.write("@%s/%d\n%s\n+\n%s\n" % (aln.qname, end_id, sequence, quality))


def is_hq(aln, chr_tid, chr_start, chr_end):
    return aln.is_unmapped or aln.mapq>0 or (not (aln.tid==chr_tid and  chr_start<=aln.pos<=chr_end))

def all_pair(aln, mate, chr_tid, chr_start, chr_end):
    return True


def all_pair_hq(aln, mate, chr_tid, chr_start, chr_end):
    return  is_hq(aln, chr_tid, chr_start, chr_end) and is_hq(mate, chr_tid, chr_start, chr_end)

def get_nm(aln):
    nm_str = aln.opt("NM")
    return int(nm_str) if nm_str else 0

def perfect_aln(aln):
    return not aln.is_unmapped and aln.is_proper_pair and len(aln.cigar) == 1 and get_nm(aln) <= EXTRACTION_MAX_NM


def non_perfect(aln, mate, chr_tid, chr_start, chr_end):
    return not (perfect_aln(aln) and perfect_aln(mate))


def non_perfect_hq(aln, mate, chr_tid, chr_start, chr_end):
    return (not (perfect_aln(aln) and perfect_aln(mate))) and is_hq(aln, chr_tid, chr_start, chr_end) and is_hq(mate, chr_tid, chr_start, chr_end)



def discordant(aln, mate, chr_tid, chr_start, chr_end, isize_min=300, isize_max=400):
    if aln.tlen == 0: return True
    return not (isize_min <= abs(aln.tlen) <= isize_max)


def discordant_with_normal_orientation(aln, mate, chr_tid, chr_start, chr_end, isize_min=300, isize_max=400):
    if aln.tlen == 0: return True
    if aln.is_reverse and mate.is_reverse or not aln.is_reverse and not mate.is_reverse: return False
    return not (isize_min <= abs(aln.tlen) <= isize_max)

def get_mate(aln, bam_handles):
    mate = None
    for bam_handle in bam_handles:
        try:
            mate = bam_handle.mate(aln)
        except ValueError:
            pass
        if mate is not None:
            return mate
    return mate


def extract_read_pairs(bam_handles, region, prefix, extract_fns, pad=0, max_read_pairs = EXTRACTION_MAX_READ_PAIRS,
                       truncation_pad_read_extract = EXTRACTION_TRUNCATION_PAD,  
                       max_interval_len_truncation = EXTRACTION_MAX_INTERVAL_TRUNCATION, sv_type=''):
    logger = logging.getLogger("%s-%s" % (extract_read_pairs.__name__, multiprocessing.current_process()))

    extract_fn_names = [extract_fn.__name__ for extract_fn in extract_fns]
    logger.info("Extracting reads for region %s with padding %d using functions %s" % (
        region, pad, extract_fn_names))

    chr_name = str(region.split(':')[0])
    chr_start = int(region.split(':')[1].split("-")[0]) - pad
    chr_end = int(region.split(':')[1].split('-')[1]) + pad

    selected_pair_counts = [0] * len(extract_fn_names)
    start_time = time.time()

    if chr_start < 0:
        regions_to_extract = []
        logger.error("Skipping read extraction since interval too close to chromosome beginning")
    else:
        # Read alignments from the interval in memory and build a dictionary to get mate instead of calling bammate.mate() function
        regions_to_extract = [(chr_name, chr_start, chr_end)]
        if abs(chr_end-chr_start)>max_interval_len_truncation and sv_type in ["INV","DEL","DUP"]:
            # For large SVs, middle sequences has no effect on genotyping. So, we only extract reads around breakpoints to speed up
            truncate_start = chr_start + pad + truncation_pad_read_extract
            truncate_end = chr_end -  (pad + truncation_pad_read_extract)
            logger.info("Truncate the reads in [%d-%d] for %s_%d_%d" % (truncate_start,truncate_end,chr_name,chr_start,chr_end))
            regions_to_extract = [(chr_name, chr_start, truncate_start-1), (chr_name, truncate_end+1, chr_end)]

    aln_list = [aln for (chr_, start_, end_) in regions_to_extract for bam_handle in bam_handles for aln in bam_handle.fetch(chr_, start=start_, end=end_) if not aln.is_secondary]

    aln_dict = {}
    for aln in aln_list:
        if aln.qname not in aln_dict:
            aln_dict[aln.qname] = [None, None]
        aln_dict[aln.qname][0 if aln.is_read1 else 1] = aln

    aln_pairs = []
    if len(aln_dict) <= max_read_pairs:
        logger.info("Building mate dictionary from %d reads" % len(aln_list))
        for aln_pair in aln_dict.values():
            missing_index = 0 if aln_pair[0] is None else (1 if aln_pair[1] is None else 2)
            if missing_index < 2:
                mate = get_mate(aln_pair[1 - missing_index], bam_handles)
                if mate is not None:
                    aln_pair[missing_index] = mate
                    aln_pairs.append(aln_pair)
            else:
                aln_pairs.append(aln_pair)
    else:
        logger.info("Too many reads encountered for %s. Skipping read extraction. (%d >%d)"%(region, len(aln_dict),max_read_pairs))

    ends = [(open("%s_%s_1.fq" % (prefix, name), "w"), open("%s_%s_2.fq" % (prefix, name), "w")) for name in
            extract_fn_names]

    chr_tid = bam_handles[0].gettid(chr_name) if bam_handles else -1
    for first, second in aln_pairs:
        for fn_index, extract_fn in enumerate(extract_fns):
            if extract_fn(first, second,chr_tid,chr_start,chr_end):
                write_read(ends[fn_index][0], first)
                write_read(ends[fn_index][1], second)

                selected_pair_counts[fn_index] += 1

    for end1, end2 in ends:
        end1.close()
        end2.close()

    logger.info("Examined %d pairs in %g seconds" % (len(aln_pairs), time.time() - start_time))
    logger.info("Extraction counts %s" % (zip(extract_fn_names, selected_pair_counts)))

    return zip([(end[0].name, end[1].name) for end in ends], selected_pair_counts)


if __name__ == "__main__":
    FORMAT = '%(levelname)s %(asctime)-15s %(name)-20s %(message)s'
    logging.basicConfig(level=logging.INFO, format=FORMAT)

    parser = argparse.ArgumentParser(description="Extract reads and mates from a region for spades assembly",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--bams", nargs='+', help="BAM files to extract reads from", required=True, default=[])
    parser.add_argument("--region", help="Samtools region string", required=True)
    parser.add_argument("--prefix", help="Output FASTQ prefix", required=True)
    parser.add_argument("--extract_fn", help="Extraction function", choices=["all_pair", "non_perfect", "discordant"],
                        default="all_pair")
    parser.add_argument("--pad", help="Padding to apply on both sides of the interval", type=int, default=0)
    parser.add_argument("--isize_min", help="Minimum insert size", default=200, type=int)
    parser.add_argument("--isize_max", help="Maximum insert size", default=500, type=int)
    parser.add_argument("--max_read_pairs", help="Maximum read pairs to extract for an interval",
                        default=EXTRACTION_MAX_READ_PAIRS, type=int)

    args = parser.parse_args()

    if args.extract_fn == 'all_pair':
        extract_fn = all_pair
    elif args.extract_fn == 'non_perfect':
        extract_fn = non_perfect
    else:
        extract_fn = partial(discordant, isize_min=args.isize_min, isize_max=args.isize_max)
        update_wrapper(extract_fn, discordant)

    bam_handles = [pysam.Samfile(bam, "rb") for bam in args.bams]

    extract_read_pairs(bam_handles, args.region, args.prefix, [extract_fn], pad=args.pad,
                       max_read_pairs=args.max_read_pairs)

    for bam_handle in bam_handles:
        bam_handle.close()
