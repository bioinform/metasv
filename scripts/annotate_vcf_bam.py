#!/net/kodiak/volumes/lake/shared/users/marghoob/my_env/bin/python
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
import vcf
import random
import math

import pysam

def annotate_vcfs(bam, chromosomes, workdir, num_threads, vcfs):
    func_logger = logging.getLogger("%s-%s" % (annotate_vcfs.__name__, multiprocessing.current_process()))
    random.seed(0)

    # Load indexed BAM file
    sam_file = pysam.Samfile(bam.name, "rb")

    if not chromosomes:
        func_logger.info("Chromosome list unspecified. Inferring from the BAMs")
        chromosomes += list(sam_file.references)
        chromosomes = sorted(list(set(chromosomes)))
        func_logger.info("Chromosome list inferred as %s" % (str(chromosomes)))

    if not chromosomes or len(chromosomes) == 0:
        func_logger.error("Chromosome list empty")
        return None

    # Read through samfile and get some statistics
    # hard code this for now
    read_limit = 1000

    # this is temporary, needs to read the reference to be sensible
    num_read = 0.0
    cover_sum = 0.0
    template_list = list()
    first_chr = sam_file.getrname(0)
    for i in xrange(0, read_limit):
        loc = random.randint(0, 30000000)
        alignments = sam_file.fetch(first_chr, loc, loc+1)

        curr_num = 0
        for aln in alignments:
            if aln.mapq < 18:
                continue
            curr_num += 1
            cover_sum += 1
            template_list.append(abs(aln.tlen))

        if curr_num > 0:
            num_read += 1

    template_list.sort()
    num_template = float(len(template_list))
    low_bound = int(math.floor(num_template*0.05))
    upp_bound = int(math.ceil(num_template*0.95))

    insert_count = 0
    insert_sum = 0.0
    insert_sq_sum = 0.0

    for i in xrange(low_bound, upp_bound):
        insert_count += 1
        insert_sum += template_list[i]
        insert_sq_sum += template_list[i] * template_list[i]

    mean_coverage = cover_sum/num_read
    mean_insert_size = insert_sum/insert_count
    sd_insert_size = math.sqrt((insert_sq_sum/insert_count) - (mean_insert_size * mean_insert_size))
    func_logger.info("Estimated coverage mean:      {0:.2f}".format(mean_coverage))
    func_logger.info("Estimated template size mean: {0:.2f}".format(mean_insert_size))
    func_logger.info("Estimated template size sd:   {0:.2f}".format(sd_insert_size))
    func_logger.info("Estimated template size Q5:   {0:.2f}".format(template_list[low_bound]))
    func_logger.info("Estimated template size Q95:  {0:.2f}".format(template_list[upp_bound - 1]))

    # Read though VCF one line at a time

    for inVCF in vcfs:
        vcf_reader = vcf.Reader(open(inVCF.name))
        vcf_template_reader = vcf.Reader(open(inVCF.name))
        vcf_writer = vcf.Writer(open("anno_"+inVCF.name, 'w'), vcf_template_reader)
        num_processed = 0
        for vcf_record in vcf_reader:
            if vcf_record.CHROM not in chromosomes:
                continue
            num_processed += 1
            if num_processed % 100 == 0:
                func_logger.info("{0} read from {1}".format(num_processed,inVCF.name))

            # get the interval that corresponds to the SV
            breakpoints = (vcf_record.start, vcf_record.INFO['END'])

            process_variant = True
            if breakpoints[1] - breakpoints[0] > 1000000:
                process_variant = False

            if process_variant:
                # get reads between breakpoints
                # sample with replacement 100 points
                unique_coverage = 0.0
                total_coverage = 0.0

                num_forward = 0.0

                num_repeat = 10
                for i in xrange(0, num_repeat):
                    loc = random.randint(breakpoints[0], breakpoints[1])
                    alignments = sam_file.fetch(vcf_record.CHROM, loc, loc + 1)
                    for rec in alignments:
                        if rec.mapq >= 18:
                            unique_coverage += 1
                            if not rec.is_reverse:
                                num_forward += 1
                        total_coverage += 1

                # get coverage between the breakpoints
                vcf_record.INFO["AA_UNIQ_COV"] = (unique_coverage/num_repeat)/mean_coverage
                vcf_record.INFO["AA_TOTAL_COV"] = (total_coverage/num_repeat)/mean_coverage

                if unique_coverage > 0:
                    vcf_record.INFO["AA_TOTAL_STRAND"] = (num_forward/unique_coverage - 0.5) ** 2

            # get strand bias

            # get mapping quality stats

            # get clipped reads stats

            # get discordant reads stats

            # get supplementary alignment stats

            vcf_writer.write_record(vcf_record)

        vcf_reader.close()
        vcf_template_reader.close()
        vcf_writer.close()




if __name__ == "__main__":
    FORMAT = '%(levelname)s %(asctime)-15s %(name)-20s %(message)s'
    logging.basicConfig(level=logging.INFO, format=FORMAT)
    logger = logging.getLogger(__name__)

    parser = argparse.ArgumentParser(
        description="Annotate VCF with additional useful features",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--bam", help="BAM file", required=True, type=file)
    parser.add_argument("--chromosomes", nargs="+", help="Chromosomes", default=[])
    parser.add_argument("--workdir", help="Working directory", default="work")
    parser.add_argument("--num_threads", help="Number of threads to use", default=1, type=int)
    parser.add_argument("--vcfs", nargs="+", help="Input VCF files", type=file)

    args = parser.parse_args()

    logger.info("Command-line: " + " ".join(sys.argv))

    annotate_vcfs(args.bam, args.chromosomes, args.workdir, args.num_threads, args.vcfs)
