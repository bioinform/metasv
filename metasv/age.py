import traceback
import os
import multiprocessing
import subprocess
import hashlib
from functools import partial
import json
import base64
from external_cmd import TimedExternalCmd

import pysam
import pybedtools

from spades_contig import SpadesContig
from tigra_contig import TigraContig
from svregion import SVRegion
from age_parser import *
from process_age_alignment import process_age_records
from defaults import *

FORMAT = '%(levelname)s %(asctime)-15s %(name)-20s %(message)s'
logging.basicConfig(level=logging.INFO, format=FORMAT)
logger = logging.getLogger(__name__)


def get_age_file_prefix(contig):
    return hashlib.md5(contig.raw_name).hexdigest()


def run_cmd(cmd, logger, out, err):
    logger.info("Running command %s" % cmd)
    retcode = subprocess.call(cmd, shell=True, stderr=err, stdout=out)
    logger.info("Returned code %d" % retcode)

    return retcode


def run_age_single(intervals_bed=None, region_list=[], contig_dict={}, reference=None, assembly=None, pad=AGE_PAD,
                   age=None, truncation_pad_read_age = AGE_TRUNCATION_PAD,
                   max_interval_len_truncation_age = AGE_MAX_INTERVAL_TRUNCATION,
                   dist_to_expected_bp = AGE_DIST_TO_BP, min_del_subalign_len = MIN_DEL_SUBALIGN_LENGTH, 
                   min_inv_subalign_len = MIN_INV_SUBALIGN_LENGTH, age_window = AGE_WINDOW_SIZE,
                   age_workdir=None, timeout=AGE_TIMEOUT, keep_temp=False, myid=0):
    thread_logger = logging.getLogger("%s-%s" % (run_age_single.__name__, multiprocessing.current_process()))

    bedtools_intervals = []
    intervals_bedtool = pybedtools.BedTool(intervals_bed)

    assembly_fasta = pysam.Fastafile(assembly) if assembly else None
    reference_fasta = pysam.Fastafile(reference)

    breakpoints_bed = None

    thread_logger.info("Will process %d intervals" % (len(region_list)))

    try:
        for region in region_list:
            bedtools_interval = pybedtools.Interval(region[0], region[1], region[3])
            matching_intervals = [interval for interval in intervals_bedtool if (
                interval.start == bedtools_interval.start and interval.end == bedtools_interval.end and interval.chrom == bedtools_interval.chrom)]
            if not matching_intervals:
                thread_logger.info("Matching interval not found for %s" % (str(bedtools_interval)))
                matching_interval = bedtools_interval
            else:
                matching_interval = matching_intervals[0]
            thread_logger.info("Matching interval %s" % (str(matching_interval)))
            sc_locations = []
            try:
                sc_locations = map(int, json.loads(base64.b64decode(matching_interval.name.split(",")[0]))["SC_LOCATIONS"].split(","))
            except:
                pass

            if region not in contig_dict:
                continue
            if not contig_dict[region]:
                continue

            region_object = SVRegion(region[0], region[1], region[2], region[3])
            if region_object.pos1 - pad < 0:
                thread_logger.error("Region too close to start of chromosome. Skipping.")
                continue

            reference_sequence = reference_fasta.fetch(reference=region_object.chrom1, start=region_object.pos1 - pad,
                                                       end=region_object.pos2 + pad)
            region_name = "%s.%d.%d" % (region_object.chrom1, region_object.pos1, region_object.pos2)
            ref_name = os.path.join(age_workdir, "%s.ref.fa" % region_name)

            thread_logger.info("Writing the ref sequence for region %s" % region_name)
            with open(ref_name, "w") as file_handle:
                file_handle.write(">{}.ref\n{}".format(region_name, reference_sequence))


            

            age_records = []
            thread_logger.info("Processing %d contigs for region %s" % (len(contig_dict[region]), str(region_object)))
            for contig in contig_dict[region]:
                thread_logger.info(
                    "Writing the assembeled sequence %s of length %s" % (contig.raw_name, contig.sequence_len))
                
                tr_region=[]
                if region_object.length()>max_interval_len_truncation_age and contig.sv_type in ["INV","DEL","DUP"]:
                    # For large SVs, middle sequences has no effect on genotyping. So, we truncate middle region of reference to speed up
                    thread_logger.info("Truncate the reference sequence.")
                    

                    truncate_start = pad + dist_to_expected_bp + truncation_pad_read_age +1
                    truncate_end = len(reference_sequence) -  (pad + dist_to_expected_bp + truncation_pad_read_age)
                    reference_sequence_tr=reference_sequence[0:truncate_start-1]+reference_sequence[truncate_end:]
                    region_name_tr = "%s.%d.%d.tr_%d_%d" % (region_object.chrom1, region_object.pos1, region_object.pos2,truncate_start,truncate_end)
                    ref_name_tr = os.path.join(age_workdir, "%s.ref.fa" % region_name_tr)

                    thread_logger.info("Writing the truncated ref sequence for region %s, contig %s" % (region_name_tr, contig.raw_name))
                    with open(ref_name_tr, "w") as file_handle:
                        file_handle.write(">{}.ref\n{}".format(region_name_tr, reference_sequence_tr))
                        
                    ref_len = len(reference_sequence_tr)
                    ref_f_name = ref_name_tr
                    tr_region = [truncate_start,truncate_end-truncate_start+1]
                    
                else:
                    ref_len = region_object.length()
                    ref_f_name = ref_name
                    
                if contig.sequence_len * ref_len >= 100000000:
                    thread_logger.info("Skipping contig because AGE problem is large (contig_len = %d , ref_len= %d)"%(contig.sequence_len, ref_len))
                    continue

                contig_sequence = assembly_fasta.fetch(contig.raw_name)

                prefix = get_age_file_prefix(contig)
                asm_name = os.path.join(age_workdir, "%s.as.fa" % prefix)
                out = os.path.join(age_workdir, "%s.age.out" % prefix)
                err = os.path.join(age_workdir, "%s.age.err" % prefix)
                fd_out = open(out, "w")
                fd_err = open(err, "w")

                with open(asm_name, "w") as file_handle:
                    file_handle.write(">{}.as\n{}".format(region_name, contig_sequence))

                age_cmd = "%s %s -both -go=-6 %s %s" % (
                    age, "-inv" if contig.sv_type == "INV" else "-tdup" if contig.sv_type == "DUP" else "-indel",
                    ref_f_name, asm_name)
                cmd_runner = TimedExternalCmd(age_cmd, thread_logger)
                retcode = cmd_runner.run(timeout=timeout, cmd_log_fd_out=fd_out, cmd_log_fd_err=fd_err)
                fd_out.close()
                fd_err.close()

                if retcode == 0:
                    age_record = AgeRecord(out,tr_region_1=tr_region)
                    if len(age_record.inputs) == 2:
                        age_record.contig = contig
                        age_record.set_assembly_contig(contig_sequence)
                        age_records.append(age_record)
                    else:
                        thread_logger.error("Number of inputs != 2 in age output file %s. Skipping." % out)

                if not keep_temp:
                    os.remove(asm_name)
                    os.remove(err)
                    if tr_region:
                        os.remove(ref_name_tr)

            unique_age_records = get_unique_age_records(age_records)

            thread_logger.info("Unique %d AGE records for region %s" % (len(unique_age_records), str(region_object)))
            for age_record in unique_age_records:
                thread_logger.info(str(age_record))

            sv_types = list(set([age_record.contig.sv_type for age_record in unique_age_records]))
            if len(sv_types) != 1:
                thread_logger.error("Some problem. Mixed SV types for this interval %s" % (str(sv_types)))
            else:
                sv_type = sv_types[0]
                thread_logger.info("Processing region of type %s" % sv_type)
                breakpoints, info_dict = process_age_records(unique_age_records, sv_type=sv_type, 
                                                             pad=pad, dist_to_expected_bp=dist_to_expected_bp,
                                                             min_del_subalign_len=min_del_subalign_len,
                                                             min_inv_subalign_len=min_inv_subalign_len,
                                                             age_window=age_window, sc_locations=sc_locations)
                bedtools_fields = matching_interval.fields
                if len(breakpoints) == 1 and sv_type == "INS":
                    bedtools_fields += map(str, [breakpoints[0][0], breakpoints[0][0] + 1, breakpoints[0][1], breakpoints[0][2]])
                elif len(breakpoints) == 2 and (sv_type in ["DEL","INV","DUP"]):
                    bedtools_fields += map(str, breakpoints + [breakpoints[1] - breakpoints[0]] + ["."])
                else:
                    bedtools_fields += map(str, [bedtools_fields[1], bedtools_fields[2], -1, "."])
                bedtools_fields[3] += ";AS"
                bedtools_fields.append(base64.b64encode(json.dumps(info_dict)))
                thread_logger.info("Writing out fields %s" % (str(bedtools_fields)))
                bedtools_intervals.append(pybedtools.create_interval_from_list(bedtools_fields))

            if not keep_temp:
                os.remove(ref_name)
    except Exception as e:
        thread_logger.error('Caught exception in worker thread')

        # This prints the type, value, and stack trace of the
        # current exception being handled.
        traceback.print_exc()

        print()
        raise e

    if assembly_fasta:
        assembly_fasta.close()
    reference_fasta.close()

    thread_logger.info("Writing %d intervals" % (len(bedtools_intervals)))
    if bedtools_intervals:
        breakpoints_bed = os.path.join(age_workdir, "%d_breakpoints.bed" % myid)
        pybedtools.BedTool(bedtools_intervals).saveas(breakpoints_bed)

    return breakpoints_bed


def run_age_single_callback(result, result_list):
    if result is not None:
        result_list.append(result)


def run_age_parallel(intervals_bed=None, reference=None, assembly=None, pad=AGE_PAD, age=None, age_workdir=None,
                     timeout=AGE_TIMEOUT, keep_temp=False, assembly_tool="spades", chrs=[], nthreads=1,
                     min_contig_len=AGE_MIN_CONTIG_LENGTH,
                     max_region_len=AGE_MAX_REGION_LENGTH, sv_types=[], 
                     min_del_subalign_len=MIN_DEL_SUBALIGN_LENGTH, min_inv_subalign_len=MIN_INV_SUBALIGN_LENGTH,
                     age_window = AGE_WINDOW_SIZE):
    func_logger = logging.getLogger("%s-%s" % (run_age_parallel.__name__, multiprocessing.current_process()))

    if not os.path.isdir(age_workdir):
        func_logger.info("Creating %s" % age_workdir)
        os.makedirs(age_workdir)

    if assembly:
        if not os.path.isfile("%s.fai" % assembly):
            func_logger.info("Assembly FASTA wasn't indexed. Will attempt to index now.")
            pysam.faidx(assembly)

        func_logger.info("Loading assembly contigs from %s" % assembly)
        with open(assembly) as assembly_fd:
            if assembly_tool == "spades":
                contigs = [SpadesContig(line[1:]) for line in assembly_fd if line[0] == '>']
            elif assembly_tool == "tigra":
                contigs = [TigraContig(line[1:]) for line in assembly_fd if line[0] == '>']
    else:
        contigs = []

    chrs = set(chrs)
    sv_types = set(sv_types)
    contig_dict = {contig.sv_region.to_tuple(): [] for contig in contigs if (len(
        chrs) == 0 or contig.sv_region.chrom1 in chrs) and contig.sequence_len >= min_contig_len and contig.sv_region.length() <= max_region_len and (
                       len(sv_types) == 0 or contig.sv_type in sv_types)}

    func_logger.info("Generating the contig dictionary for parallel execution")
    small_contigs_count = 0
    for contig in contigs:
        if contig.sv_region.length() > max_region_len: 
            func_logger.info("Too large SV region length: %d > %d" % (contig.sv_region.length(),max_region_len))
            continue
        if (len(chrs) == 0 or contig.sv_region.chrom1 in chrs) and (len(sv_types) == 0 or contig.sv_type in sv_types):
            if contig.sequence_len >= min_contig_len:
                contig_dict[contig.sv_region.to_tuple()].append(contig)
            else:
                small_contigs_count += 1

    region_list = sorted(contig_dict.keys())
    nthreads = min(nthreads, len(region_list))

    if nthreads == 0:
        func_logger.warning("AGE not run since no contigs found")
        return None

    func_logger.info("Will process %d regions with %d contigs (%d small contigs ignored) using %d threads" % (
        len(region_list), sum([len(value) for value in contig_dict.values()]), small_contigs_count, nthreads))

    pybedtools.set_tempdir(age_workdir)
    pool = multiprocessing.Pool(nthreads)

    breakpoints_beds = []
    for i in xrange(nthreads):
        region_sublist = [region for (j, region) in enumerate(region_list) if (j % nthreads) == i]
        kwargs_dict = {"intervals_bed": intervals_bed, "region_list": region_sublist, "contig_dict": contig_dict,
                       "reference": reference, "assembly": assembly, "pad": pad, "age": age, "age_workdir": age_workdir,
                       "timeout": timeout, "keep_temp": keep_temp, "myid": i, 
                       "min_del_subalign_len": min_del_subalign_len, "min_inv_subalign_len": min_inv_subalign_len,
                       "age_window" : age_window}
        pool.apply_async(run_age_single, args=[], kwds=kwargs_dict,
                         callback=partial(run_age_single_callback, result_list=breakpoints_beds))

    pool.close()
    pool.join()

    func_logger.info("Finished parallel execution")

    func_logger.info("Will merge the following breakpoints beds %s" % (str(breakpoints_beds)))

    pybedtools.cleanup(remove_all=True)

    if not breakpoints_beds:
        return None

    bedtool = pybedtools.BedTool(breakpoints_beds[0])
    for bed_file in breakpoints_beds[1:]:
        bedtool = bedtool.cat(pybedtools.BedTool(bed_file), postmerge=False)

    bedtool = bedtool.moveto(os.path.join(age_workdir, "breakpoints_unsorted.bed"))
    merged_bed = os.path.join(age_workdir, "breakpoints.bed")
    bedtool.sort().saveas(merged_bed)

    return merged_bed