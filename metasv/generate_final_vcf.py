import argparse
import sys
import datetime
import os
from collections import OrderedDict
import json
import base64
import logging

import pybedtools
import pysam
import vcf

import fasta_utils

mydir = os.path.dirname(os.path.realpath(__file__))
vcf_template = os.path.join(mydir, "resources/template.vcf")


def check_duplicates(interval1,interval2,max_dist=10):
    if interval1.chrom != interval2.chrom or \
       abs(interval1.start-interval2.start)>max_dist or \
       abs(interval1.end-interval2.end)>max_dist or \
       interval1.fields[3].split(",")[1] != interval2.fields[3].split(",")[1]:
        return None
    info1=json.loads(base64.b64decode(interval1.fields[3].split(",")[0]))
    info2=json.loads(base64.b64decode(interval2.fields[3].split(",")[0]))
    svmethods= sorted(list(set(info1["SVMETHOD"]+info2["SVMETHOD"])))
    sources = []
    if "SOURCES" in info1:
        sources.append(info1["SOURCES"])
    if "SOURCES" in info2:
        sources.append(info2["SOURCES"])
    sources=",".join(sources)
    if sources:
        info1["SOURCES"]=sources
    if "PASS" in [interval1.fields[7], interval2.fields[7]] or ("AS" not in svmethods and len(set(svmethods) - {"SC","AS"})>1):
        sv_filter = "PASS"
    else:
        sv_filter = "LowQual"
            
    end = max(interval1.end,interval2.end)
    start = min(interval1.start,interval2.start)
    info1.update(
        {"END": end, "SVMETHOD": svmethods, "NUM_SVMETHODS": len(svmethods)})
    return pybedtools.Interval(interval1.chrom, start,end, name="%s,%s,%d,%s" % (
                               base64.b64encode(json.dumps(info1)), 
                               info1["SVTYPE"], end-start,
                               ";".join(svmethods)), 
                               score = interval1.score, 
                               otherfields=[interval1.fields[6],sv_filter])



def get_interval_info(feature,pass_calls):
    func_logger = logging.getLogger("%s" % (get_interval_info.__name__))
    
    pos = feature.start
    end = feature.end
    genotype = "./." if len(feature.fields) < 11 else feature.fields[10]
    if genotype == "0/0":
        func_logger.info("Skipping homozygous reference %s" % str(feature))
        return None

    sub_names = feature.name.split(":")
    sub_lengths = map(lambda x: int(x.split(",")[2]), sub_names)

    sub_types = map(lambda x: x.split(",")[1], sub_names)
    sub_methods = [name.split(",")[3] for name in sub_names]
    svmethods = (";".join([name.split(",")[3] for name in sub_names])).split(";")
    try:
        info = json.loads(base64.b64decode(name.split(",")[0]))
    except TypeError:
        info = dict()
    if len(feature.fields) > 9:
        info.update(json.loads(base64.b64decode(feature.fields[9])))

    index_to_use = 0
    is_pass = False
    svlen = -1
    if "DEL" in sub_types:
        index_to_use = sub_types.index("DEL")
        svmethods_s = set(svmethods) - {"SC","AS"}
        is_pass = len(svmethods_s) > 1
        if "AS" in svmethods:
            pos, end, svlen = map(int, feature.fields[6:9])
            is_pass = svlen >= 100
    elif "INV" in sub_types:
        index_to_use = sub_types.index("INV")
        svmethods_s = set(svmethods) - {"SC","AS"}
        is_pass = len(svmethods_s) > 1
        if "AS" in svmethods:
            pos, end, svlen = map(int, feature.fields[6:9])
            is_pass = svlen >= 100               
    elif "INS" in sub_types and ("SC" in svmethods or "AS" in svmethods):
        # TODO: I think it should be sub_types.index
        #index_to_use = [i for i,methods in enumerate(sub_methods) if ("SC" in methods) or ("AS" in svmethods)][0]
        index_to_use = sub_types.index("INS")
        pos, end, svlen = map(int, feature.fields[6:9])
    elif "ITX" in sub_types:
        index_to_use = sub_types.index("ITX")
        svmethods_s = set(svmethods) - {"SC","AS"}
        is_pass = len(svmethods_s) > 1
        end = info["END"]
    elif "CTX" in sub_types:
        index_to_use = sub_types.index("CTX")
        svmethods_s = set(svmethods) - {"SC","AS"}
        is_pass = len(svmethods_s) > 1
        end = info["END"]
        
    if pos < 1:
        func_logger.info("Variant with pos < 1 encountered. Skipping! %s" % str(feature))
        return None

    if svlen < 0: svlen = sub_lengths[index_to_use]
    if sub_types[index_to_use] == "DEL":
        svlen = -svlen

    sv_type = sub_types[index_to_use]
    if sv_type == "INS":
        if pass_calls and end != pos + 1:
            return None
        end = pos
        is_pass = (int(feature.fields[8]) != -1) and (svlen == 0 or svlen >= 100)
    info.update(
        {"END": end, "SVLEN": svlen, "SVTYPE": sv_type, "SVMETHOD": svmethods, "NUM_SVMETHODS": len(svmethods)})
    sv_filter = "PASS" if is_pass else "LowQual"
    interval_info={"pos": pos, "end": end, "info": info, "sv_type": sv_type, "genotype": genotype,
                   "sv_length": abs(svlen), "svmethods": svmethods, "sv_filter": sv_filter}    
    return interval_info

def filter_confused_INS_calls(nonfilterd_bed, filterd_bed, wiggle=20):

    nonins_intervals=[]
    bedtool = pybedtools.BedTool(nonfilterd_bed)
    bedtool_INS = bedtool.filter(lambda x: "INS" in x.name.split(",")[1]).sort()
    bedtool_others = bedtool.filter(lambda x: "INS" not in x.name.split(",")[1]).sort()
    bedtool_good_nonINS = bedtool_others.filter(lambda x: ("DEL" in x.name.split(",")[1] or "INV" in x.name.split(",")[1]) and x.fields[7] != "LowQual").saveas()
    nonINS_bp_intervals=[]
    for interval in bedtool_good_nonINS:
        start=interval.start
        end=interval.end
        nonINS_bp_intervals.append(pybedtools.Interval(interval.chrom,max(start-wiggle,0),start+wiggle))
        nonINS_bp_intervals.append(pybedtools.Interval(interval.chrom,max(end-wiggle,0),end+wiggle))

    bedtool_bp_nonINS=pybedtools.BedTool(nonINS_bp_intervals) 
    bad_INS=bedtool_INS.window(bedtool_bp_nonINS,w=wiggle)

    bedtool_filtered=bedtool_INS.window(bad_INS,w=wiggle,v=True).saveas(filterd_bed)
    bedtool_filtered =bedtool_filtered.cat(bedtool_others,postmerge=False).sort().saveas(filterd_bed)
    return filterd_bed

def convert_metasv_bed_to_vcf(bedfile=None, vcf_out=None, workdir=None, vcf_template_file=vcf_template, sample=None, reference=None,
                              pass_calls=True):
    func_logger = logging.getLogger("%s" % (convert_metasv_bed_to_vcf.__name__))

    if bedfile:
        intervals = []
        for interval in pybedtools.BedTool(bedfile):
            interval_info = get_interval_info(interval,pass_calls)            
            if interval_info:
                updated_interval = pybedtools.Interval(interval.chrom, interval_info["pos"], 
                                                       interval_info["end"], name="%s,%s,%d,%s" % (
                                                       base64.b64encode(json.dumps(interval_info["info"])), 
                                                       interval_info["sv_type"], interval_info["sv_length"],
                                                       ";".join(interval_info["svmethods"])), 
                                                       score = interval.score, 
                                                       otherfields=[interval_info["genotype"]
                                                                    , interval_info["sv_filter"]])
                if not intervals:
                    intervals.append(updated_interval)
                else:
                    merged_interval=check_duplicates(updated_interval,intervals[-1])
                    if merged_interval:
                        func_logger.info("Merging intervals: %s and %s" % (updated_interval,intervals[-1]))
                        intervals.pop()
                        intervals.append(merged_interval)
                    else:
                        intervals.append(updated_interval)
            else: 
                func_logger.info("Skip interval: %s" % (interval))

    nonfilterd_bed = os.path.join(workdir, "final_nonfilterd.bed")
    filterd_bed = os.path.join(workdir, "final_filterd.bed")
    bedtool = pybedtools.BedTool(intervals).sort().moveto(nonfilterd_bed)
    filterd_bed = filter_confused_INS_calls(nonfilterd_bed,filterd_bed)    

    vcf_template_reader = vcf.Reader(open(vcf_template_file, "r"))
    # The following are hacks to ensure sample name and contig names are put in the VCF header
    vcf_template_reader.samples = [sample]
    contigs = []
    if reference:
        contigs = fasta_utils.get_contigs(reference)
        contigs_order_dict = {contig.name: index for (index, contig) in enumerate(contigs)}
        vcf_template_reader.contigs = OrderedDict([(contig.name, (contig.name, contig.length)) for contig in contigs])
        vcf_template_reader.metadata["reference"] = reference

    vcf_template_reader.metadata["fileDate"] = str(datetime.date.today())
    vcf_template_reader.metadata["source"] = [" ".join(sys.argv)]
    vcf_writer = vcf.Writer(open(vcf_out, "w"), vcf_template_reader)
    vcf_records = []
    if filterd_bed:
        bedtool = pybedtools.BedTool(filterd_bed)
        for interval in bedtool:
            name_split=interval.name.split(",")
            info = json.loads(base64.b64decode(name_split[0]))
            sv_type = name_split[1]
            sv_id = "."
            ref = "."
            alt = [vcf.model._SV(sv_type)]
            qual = "."
            sv_filter = [interval.fields[7]]
            genotype = interval.fields[6]
            sv_format = "GT"
            sample_indexes = [0]
            vcf_record = vcf.model._Record(interval.chrom, interval.start, sv_id, ref, alt, qual,
                                           sv_filter, info, sv_format, sample_indexes)
            vcf_record.samples = vcf_template_reader._parse_samples([genotype], "GT", vcf_record)
            vcf_records.append(vcf_record)
            
    if contigs:
        vcf_records.sort(key=lambda x: (contigs_order_dict[x.CHROM], x.POS))
    else:
        vcf_records.sort(key=lambda x: (x.CHROM, x.POS))

    for vcf_record in vcf_records:
        vcf_writer.write_record(vcf_record)
    vcf_writer.close()

    func_logger.info("Tabix compressing and indexing %s" % vcf_out)
    pysam.tabix_index(vcf_out, force=True, preset="vcf")


if __name__ == "__main__":
    FORMAT = '%(levelname)s %(asctime)-15s %(name)-20s %(message)s'
    logging.basicConfig(level=logging.INFO, format=FORMAT)
    logger = logging.getLogger(__name__)

    parser = argparse.ArgumentParser(description="Convert MetaSV final BED to VCF",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("--sample", help="Sample name", required=True)
    parser.add_argument("--bed", help="MetaSV final BED", required=True)
    parser.add_argument("--vcf", help="Final VCF to output", required=True)
    parser.add_argument("--vcf_template", help="VCF template", default=vcf_template)
    parser.add_argument("--reference", help="Reference FASTA", required=False)
    parser.add_argument("--work", help="Work directory", default="work")
    parser.add_argument("--pass_only", action="store_true", help="Output only PASS calls")

    args = parser.parse_args()

    convert_metasv_bed_to_vcf(bedfile=args.bed, vcf_out=args.vcf, workdir=args.work, vcf_template_file=args.vcf_template,
                              sample=args.sample,
                              reference=args.reference, pass_calls=args.pass_only)
