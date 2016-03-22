import argparse
import sys
import datetime
import os
from collections import OrderedDict
import json
import base64
import logging
from functools import partial

import pybedtools
import pysam
import vcf
import fasta_utils

mydir = os.path.dirname(os.path.realpath(__file__))
VCF_TEMPLATE = os.path.join(mydir, "resources/template.vcf")


def check_duplicates(interval1, interval2, max_dist=10):
    if interval1.chrom != interval2.chrom or \
                    abs(interval1.start - interval2.start) > max_dist or \
                    abs(interval1.end - interval2.end) > max_dist or \
                    interval1.fields[3].split(",")[1] != \
                    interval2.fields[3].split(",")[1]:
        return None
    info1 = json.loads(base64.b64decode(interval1.fields[3].split(",")[0]))
    info2 = json.loads(base64.b64decode(interval2.fields[3].split(",")[0]))
    svmethods = sorted(list(set(info1["SVMETHOD"] + info2["SVMETHOD"])))
    sources = []
    if "SOURCES" in info1:
        sources.append(info1["SOURCES"])
    if "SOURCES" in info2:
        sources.append(info2["SOURCES"])
    sources = ",".join(sources)
    if sources:
        info1["SOURCES"] = sources
    if "PASS" in [interval1.fields[7], interval2.fields[7]] or (
            "AS" not in svmethods and len(set(svmethods) - {"SC", "AS"}) > 1):
        sv_filter = "PASS"
    else:
        sv_filter = "LowQual"

    end = max(interval1.end, interval2.end)
    start = min(interval1.start, interval2.start)
    info1.update(
        {"END": end, "SVMETHOD": svmethods, "NUM_SVMETHODS": len(svmethods)})
    return pybedtools.Interval(interval1.chrom, start, end,
                               name="%s,%s,%d,%s" % (
                                   base64.b64encode(json.dumps(info1)),
                                   info1["SVTYPE"], end - start,
                                   ";".join(svmethods)),
                               score=interval1.score,
                               otherfields=[interval1.fields[6], sv_filter])


def get_interval_info(feature, pass_calls):
    func_logger = logging.getLogger("%s" % get_interval_info.__name__)

    pos = feature.start
    end = feature.end
    genotype = "./." if len(feature.fields) < 12 else feature.fields[11]
    if genotype == "0/0":
        func_logger.info("Skipping homozygous reference %s" % str(feature))
        return None

    sub_names = feature.name.split(":")
    sub_lengths = map(lambda x: int(x.split(",")[2]), sub_names)

    sub_types = map(lambda x: x.split(",")[1], sub_names)
    sub_methods = [name.split(",")[3] for name in sub_names]
    svmethods = (";".join([name.split(",")[3] for name in sub_names])).split(
        ";")
    try:
        info = json.loads(base64.b64decode(name.split(",")[0]))
    except TypeError:
        info = dict()
    if len(feature.fields) > 10:
        info.update(json.loads(base64.b64decode(feature.fields[10])))

    index_to_use = 0
    is_pass = False
    svlen = -1

    sv_type = None
    for sub_type in ["DEL", "INV", "DUP", "ITX", "CTX", "INS"]:
        if sub_type in sub_types:
            index_to_use = sub_types.index(sub_type)
            sv_type = sub_type
            break

    if not sv_type:
        func_logger.info("Unknown SV type! %s" % str(feature))
        return None

    if sv_type != "INS":
        svmethods_s = set(svmethods) - {"SC", "AS"}
        is_pass = len(svmethods_s) > 1
        if "AS" in svmethods:
            pos, end, svlen = map(int, feature.fields[6:9])
            is_pass = svlen >= 100
        if sv_type in ["ITX", "CTX"]:
            end = info["END"]

        if svlen < 0: svlen = sub_lengths[index_to_use]
        if sv_type == "DEL":
            svlen = -svlen
    else:
        if "SC" in svmethods or "AS" in svmethods:
            # TODO: I think it should be sub_types.index
            # index_to_use = [i for i,methods in enumerate(sub_methods) if ("SC" in methods) or ("AS" in svmethods)][0]
            pos, end, svlen = map(int, feature.fields[6:9])
        if svlen < 0: svlen = sub_lengths[index_to_use]
        if pass_calls and end != pos + 1:
            return None
        end = pos
        is_pass = (int(feature.fields[8]) != -1) and (
        svlen == 0 or svlen >= 100) and (
                  ("SC" in svmethods or "AS" in svmethods) or (
                  len(set(svmethods) - {"SC", "AS"}) > 1))

    if pos < 1:
        func_logger.info(
            "Variant with pos < 1 encountered. Skipping! %s" % str(feature))
        return None

    info.update(
        {"END": end, "SVLEN": svlen, "SVTYPE": sv_type, "SVMETHOD": svmethods,
         "NUM_SVMETHODS": len(svmethods)})
    if "IMPRECISE" in info:
        info["IMPRECISE"] = True
    sv_filter = "PASS" if is_pass else "LowQual"
    interval_info = {"pos": pos, "end": end, "info": info, "sv_type": sv_type,
                     "genotype": genotype,
                     "sv_length": abs(svlen), "svmethods": svmethods,
                     "sv_filter": sv_filter}
    return interval_info


def filter_confused_INS_calls(nonfilterd_bed, filterd_bed, wiggle=20):
    nonins_intervals = []
    bedtool = pybedtools.BedTool(nonfilterd_bed)
    bedtool_INS = bedtool.filter(lambda x: "INS" in x.name.split(",")[1]).sort()
    bedtool_others = bedtool.filter(
        lambda x: "INS" not in x.name.split(",")[1]).sort()
    bedtool_good_nonINS = bedtool_others.filter(lambda x: ("DEL" in
                                                           x.name.split(",")[
                                                               1] or "INV" in
                                                           x.name.split(",")[
                                                               1]) and x.fields[
                                                                           7] != "LowQual").saveas()
    nonINS_bp_intervals = []
    for interval in bedtool_good_nonINS:
        start = interval.start
        end = interval.end
        nonINS_bp_intervals.append(
            pybedtools.Interval(interval.chrom, max(start - wiggle, 0),
                                start + wiggle))
        nonINS_bp_intervals.append(
            pybedtools.Interval(interval.chrom, max(end - wiggle, 0),
                                end + wiggle))

    bedtool_bp_nonINS = pybedtools.BedTool(nonINS_bp_intervals)
    bad_INS = bedtool_INS.window(bedtool_bp_nonINS, w=wiggle)

    bedtool_filtered = bedtool_INS.window(bad_INS, w=wiggle, v=True).saveas()
    if len(bedtool_filtered) == 0:
        bedtool_others.saveas(filterd_bed)
    elif len(bedtool_others) == 0:
        bedtool_filtered.saveas(filterd_bed)
    else:
        bedtool_filtered = bedtool_filtered.cat(bedtool_others,
                                                postmerge=False).sort().saveas(
            filterd_bed)

    return filterd_bed


def find_idp(feature, wiggle):
    n = len(feature.fields) / 2
    if feature.chrom != feature.fields[n]:
        return None
    start_dup = feature.start
    end_dup = feature.end
    start_del = int(feature.fields[n + 1])
    end_del = int(feature.fields[n + 2])
    if abs(start_del - end_del) > (abs(start_dup - end_dup) - wiggle):
        return None
    dist_ends = [abs(start_del - start_dup), abs(end_del - end_dup)]
    if min(dist_ends) > wiggle:
        return None
    del_pos = start_del if dist_ends[0] > dist_ends[1] else end_del
    name = "%s,%s" % (feature.name, feature.fields[n + 3])
    score = "%s,%s" % (feature.score, feature.fields[n + 4])
    return pybedtools.Interval(feature.chrom, feature.start, feature.end,
                               name=name, score=score,
                               otherfields=["%d" % del_pos,
                                            "%d-%d" % (start_del, end_del)])


def find_itx(feature, wiggle):
    n = len(feature.fields) / 2
    start_idp1 = feature.start
    end_idp1 = feature.end
    start_idp2 = int(feature.fields[n + 1])
    end_idp2 = int(feature.fields[n + 2])
    dist_ends = [abs(start_idp1 - start_idp2), abs(end_idp1 - end_idp2)]
    if min(dist_ends) > wiggle:
        return None
    del_pos1 = int(feature.fields[6])
    del_pos2 = int(feature.fields[n + 6])
    if abs(del_pos1 - del_pos2) > wiggle:
        return None

    del_interval1 = map(int, feature.fields[7].split("-"))
    del_interval2 = map(int, feature.fields[n + 7].split("-"))
    lr_1 = 1 if abs(del_pos1 - del_interval1[0]) < abs(
        del_pos1 - del_interval1[1]) else 0
    lr_2 = 1 if abs(del_pos2 - del_interval2[0]) < abs(
        del_pos2 - del_interval2[1]) else 0
    if lr_1 == lr_2 or lr_2 < lr_1:
        return None

    del_id_2 = feature.name.split(",")[-1]
    del_filter_2 = feature.score.split(",")[-1]
    name = "%s,%s" % (feature.name, del_id_2)
    score = "%s,%s" % (feature.score, del_filter_2)

    return pybedtools.Interval(feature.chrom, feature.start, feature.end,
                               name=name, score=score,
                               otherfields=["%d" % ((del_pos1 + del_pos2) / 2),
                                            "%d-%d,%d-%d" % (
                                            del_interval1[0], del_interval1[1],
                                            del_interval2[0],
                                            del_interval2[1])])


def build_chr2_ins(feature, thr_top=0.15):
    sc_chr2_str = feature.fields[6]
    if sc_chr2_str == ".":
        return []
    sub_str = map(lambda x: [x.split(";")[0], map(int, x.split(";")[1:])],
                  sc_chr2_str.split(","))
    chr2_dict = {}
    for chr2, poses in sub_str:
        if chr2 not in chr2_dict:
            chr2_dict[chr2] = []
        chr2_dict[chr2].append(poses)

    chr2_dict = {k: [sum(map(lambda x: x[0], v)), min(map(lambda x: x[1], v)),
                     max(map(lambda x: x[2], v))] for k, v in
                 chr2_dict.iteritems()}
    sorted_chr2 = sorted(chr2_dict.items(), key=lambda x: x[1][0], reverse=True)
    n_reads = sum(map(lambda x: x[1][0], sorted_chr2))
    top_chr2s = filter(
        lambda x: x[1][0] > (thr_top * n_reads) and x[0] not in ["-1",
                                                                 feature.chrom],
        sorted_chr2)
    if not top_chr2s:
        return []
    ctx_intervals = []
    for chr2, [cnt, start, end] in top_chr2s:
        ctx_intervals.append(pybedtools.Interval(chr2, start, end,
                                                 name=feature.name,
                                                 score=feature.score))
    return ctx_intervals


def find_ctx(feature, overlap_ratio=0.9):
    n = len(feature.fields) / 2
    start_del_ins = int(feature.fields[n + 1])
    end_del_ins = int(feature.fields[n + 2])
    name = "%s,%s" % (feature.name, feature.fields[n + 3])
    score = "%s,%s" % (feature.score, feature.fields[n + 4])
    return pybedtools.Interval(feature.chrom, feature.start, feature.end,
                               name=name, score=score,
                               otherfields=[".", "%d-%d" % (
                               start_del_ins, end_del_ins)])


def extract_del_interval(feature):
    start, end = map(int, feature.fields[7].split("-"))
    return pybedtools.Interval(feature.chrom, start, end)


def filter_itxs(feature):
    n = len(feature.fields) / 2
    del_interval_idp = map(int, feature.fields[7].split("-"))
    del_interval_itx_1 = map(int,
                             feature.fields[n + 7].split(",")[0].split("-"))
    del_interval_itx_2 = map(int,
                             feature.fields[n + 7].split(",")[1].split("-"))
    if filter(lambda x: abs(x[0] - del_interval_idp[0]) + abs(
                    x[1] - del_interval_idp[1]) == 0,
              [del_interval_itx_1, del_interval_itx_2]) and "LowQual" not in \
            feature.fields[n + 4]:
        return None
    return pybedtools.Interval(feature.chrom, feature.start, feature.end,
                               name=feature.name,
                               score=feature.score,
                               otherfields=feature.fields[6:n])


def merge_idp_itx(fasta_file, record_dup, records_del, del_pos, del_interval,
                  score, svtype):
    info = {}
    info.update(record_dup.INFO)
    start = int(record_dup.POS)
    end = info["END"]
    dup_interval = "%s-%d-%d" % (record_dup.CHROM, start, end)
    if svtype == "IDP":
        del_interval_ends = map(int, del_interval.split("-"))
        if abs(del_pos - del_interval_ends[0]) < abs(
                        del_pos - del_interval_ends[1]):
            pos = start
            info["END"] = max(del_pos, pos)
            info["POS2"] = end
        else:
            pos = del_pos
            info["END"] = max(end, pos)
            info["POS2"] = start
    elif svtype == "ITX":
        pos = start
        info["END"] = max(del_pos, pos)
        info["POS2"] = end

    info["CHR2"] = record_dup.CHROM
    info["SVLEN"] = max(info["END"] - pos, 0)
    info["%s_INTERVALS" % svtype] = "DUP-%s," % dup_interval + ",".join(
        map(lambda x: "DEL-%s-%s" % (record_dup.CHROM, x),
            del_interval.split(",")))
    info["SVTYPE"] = svtype
    info["SVMETHOD"] = list(
        set(reduce(lambda y, z: y + z, map(lambda x: x.INFO["SVMETHOD"]
                                           , [record_dup] + records_del))))
    info["SOURCES"] = ",".join(
        map(lambda x: x.INFO["SOURCES"], [record_dup] + records_del))
    info["NUM_SVMETHODS"] = len(info["SVMETHOD"])
    info["NUM_SVTOOLS"] = len(
        set(map(lambda x: x.split('-')[-1], info["SOURCES"].split(','))))
    sv_id = "."
    ref = fasta_file.fetch(record_dup.CHROM, pos,
                           pos + 1) if fasta_file else "."
    alt = [vcf.model._SV(svtype)]
    qual = "."
    sv_filter = ["PASS"] if "LowQual" not in score else ["LowQual"]
    sv_format = "GT"
    sample_indexes = [0]
    vcf_record = vcf.model._Record(record_dup.CHROM, pos, sv_id, ref, alt, qual,
                                   sv_filter, info, sv_format, sample_indexes)
    vcf_record.samples = record_dup.samples

    return vcf_record


def merge_ctx(fasta_file, record_del, record_ins, score):
    info = {}
    info.update(record_del.INFO)
    start = int(record_del.POS)
    end = info["END"]
    del_interval = "%s-%d-%d" % (record_del.CHROM, start, end)
    ins_interval = "%s-%d-%d" % (
    record_ins.CHROM, record_ins.POS, record_ins.POS)
    pos = start
    info["POS2"] = record_ins.POS
    info["CHR2"] = record_ins.CHROM
    info["CTX_INTERVALS"] = "DEL-%s,INS-%s" % (del_interval, ins_interval)
    info["SVTYPE"] = "CTX"
    info["SVMETHOD"] = list(
        set(reduce(lambda y, z: y + z, map(lambda x: x.INFO["SVMETHOD"]
                                           , [record_del, record_ins]))))
    info["SOURCES"] = ",".join(
        map(lambda x: x.INFO["SOURCES"], [record_del, record_ins]))
    info["NUM_SVMETHODS"] = len(info["SVMETHOD"])
    info["NUM_SVTOOLS"] = len(
        set(map(lambda x: x.split('-')[-1], info["SOURCES"].split(','))))
    sv_id = "."
    ref = fasta_file.fetch(record_del.CHROM, pos,
                           pos + 1) if fasta_file else "."
    alt = [vcf.model._SV("CTX")]
    qual = "."
    sv_filter = ["PASS"] if "LowQual" not in score else ["LowQual"]
    sv_format = "GT"
    sample_indexes = [0]
    vcf_record = vcf.model._Record(record_del.CHROM, pos, sv_id, ref, alt, qual,
                                   sv_filter, info, sv_format, sample_indexes)
    vcf_record.samples = record_del.samples

    return vcf_record


def remove_info_fields(record, fields):
    info = {}
    info.update(record.INFO)
    for field in fields:
        if field in info:
            del info[field]
    sample_indexes = [0]
    vcf_record = vcf.model._Record(record.CHROM, record.POS, record.ID,
                                   record.REF, record.ALT, record.QUAL,
                                   record.FILTER, info, record.FORMAT,
                                   sample_indexes)
    vcf_record.samples = record.samples
    return vcf_record


def resolve_for_IDP_ITX_CTX(vcf_records, fasta_file, pad=0, wiggle=10,
                            overlap_ratio=0.9):
    del_records = filter(lambda x: (x.INFO["SVTYPE"] == "DEL"), vcf_records)
    dup_records = filter(lambda x: (x.INFO["SVTYPE"] == "DUP"), vcf_records)
    ins_records = filter(lambda x: (x.INFO["SVTYPE"] == "INS"), vcf_records)
    other_records = filter(
        lambda x: (x.INFO["SVTYPE"] not in ["DEL", "DUP", "INS"]), vcf_records)
    del_bedtool = pybedtools.BedTool(
        [pybedtools.Interval(x.CHROM, x.POS, (x.POS + abs(x.INFO["SVLEN"])),
                             name="DEL_%d" % i, score=x.FILTER[0]) for i, x in
         enumerate(del_records)])
    dup_bedtool = pybedtools.BedTool(
        [pybedtools.Interval(x.CHROM, x.POS, (x.POS + abs(x.INFO["SVLEN"])),
                             name="DUP_%d" % i, score=x.FILTER[0]) for i, x in
         enumerate(dup_records)])
    ins_bedtool = pybedtools.BedTool(
        [pybedtools.Interval(x.CHROM, x.POS, (x.POS + 1),
                             name="INS_%d" % i, score=x.FILTER[0],
                             otherfields=[x.INFO["SC_CHR2_STR"] if
                                          "SC_CHR2_STR" in x.INFO else "."])
         for i, x in enumerate(ins_records)])
    chr2_intervals = []
    for interval in ins_bedtool:
        chr2_intervals.extend(build_chr2_ins(interval))

    chr2_ins_bedtool = pybedtools.BedTool(chr2_intervals).sort()

    idp_bedtool = pybedtools.BedTool([])
    remained_dup_bedtool = pybedtools.BedTool([])
    if len(dup_bedtool):
        idp_bedtool = dup_bedtool.window(del_bedtool, w=wiggle).each(
            partial(find_idp, wiggle=wiggle))
        remained_dup_bedtool = dup_bedtool.intersect(idp_bedtool, f=0.95, r=True,
                                                     wa=True, v=True)
    if len(idp_bedtool):
        idp_bedtool = idp_bedtool.sort()
    if len(remained_dup_bedtool):
        remained_dup_bedtool = remained_dup_bedtool.sort()

    remained_del_bedtool = pybedtools.BedTool([])
    if len(del_bedtool):
        remained_del_bedtool = del_bedtool.intersect(
            idp_bedtool.each(partial(extract_del_interval)).sort(), f=0.95, r=True,
            wa=True, v=True)

    itx_bedtool = pybedtools.BedTool([])
    remained_idp_bedtool_1 = pybedtools.BedTool([])
    remained_idp_bedtool_2 = pybedtools.BedTool([])
    if len(idp_bedtool):
        itx_bedtool = idp_bedtool.window(idp_bedtool, w=wiggle).each(
            partial(find_itx, wiggle=wiggle))
        remained_idp_bedtool_1 = idp_bedtool.window(itx_bedtool, w=wiggle).each(
            partial(filter_itxs))
        remained_idp_bedtool_2 = idp_bedtool.window(itx_bedtool, w=wiggle,
                                                    c=True).filter(
            lambda x: x.fields[-1] == "0")
    if len(itx_bedtool):
        itx_bedtool = itx_bedtool.sort()
    if len(remained_idp_bedtool_1):
        remained_idp_bedtool_1 = remained_idp_bedtool_1.sort()
    if len(remained_idp_bedtool_2):
        remained_idp_bedtool_2 = remained_idp_bedtool_2.sort()

    ctx_bedtool = pybedtools.BedTool([])
    if len(remained_del_bedtool):
        ctx_bedtool = remained_del_bedtool.intersect(chr2_ins_bedtool, r=True,
                                                     f=overlap_ratio, wa=True,
                                                     wb=True).each(
            partial(find_ctx, overlap_ratio=overlap_ratio))
        remained_del_bedtool = remained_del_bedtool.intersect(ctx_bedtool, f=0.95,
                                                              r=True, wa=True,
                                                              v=True)
    if len(ctx_bedtool):
        ctx_bedtool = ctx_bedtool.sort()
    if len(remained_del_bedtool):
        remained_del_bedtool = remained_del_bedtool.sort()

    if len(remained_idp_bedtool_2):
        remained_idp_bedtool_2 = remained_idp_bedtool_2.cut(
            range(idp_bedtool.field_count())).sort()

    recoverd_pass_del_dup_ins = []
    removed_pass_del_dup_ins = []
    for bed in remained_idp_bedtool_1, remained_idp_bedtool_2, itx_bedtool, ctx_bedtool:
        recoverd_pass_del_dup_ins.append(",".join(
            map(lambda y: y.name, filter(lambda x: "LowQual" in x.score, bed))))
        removed_pass_del_dup_ins.append(",".join(map(lambda y: y.name, filter(
            lambda x: "LowQual" not in x.score, bed))))

    recoverd_pass_del_dup_ins = set(
        (",".join(recoverd_pass_del_dup_ins)).split(",")) - set([''])
    removed_pass_del_dup_ins = set(
        (",".join(removed_pass_del_dup_ins)).split(",")) - set([''])
    recoverd_pass_del_dup_ins = recoverd_pass_del_dup_ins - removed_pass_del_dup_ins

    recoverd_dups = list(set([x.name for x in remained_dup_bedtool]) | set(
        filter(lambda x: "DUP" in x, recoverd_pass_del_dup_ins)))
    recoverd_dels = list(set([x.name for x in remained_del_bedtool]) | set(
        filter(lambda x: "DEL" in x, recoverd_pass_del_dup_ins)))
    recoverd_inss = list(set([x.name for x in ins_bedtool]) - (
    set(filter(lambda x: "INS" in x, removed_pass_del_dup_ins))))

    vcf_records = other_records + [dup_records[int(x.split("_")[-1])] for x in
                                   recoverd_dups] + \
                  [del_records[int(x.split("_")[-1])] for x in recoverd_dels] + \
                  [ins_records[int(x.split("_")[-1])] for x in recoverd_inss] + \
                  [merge_idp_itx(fasta_file, dup_records[
                      int(x.name.split(",")[0].split("_")[-1])],
                                 [del_records[int(
                                     x.name.split(",")[1].split("_")[-1])]],
                                 int(x.fields[6]), x.fields[7], x.score, "IDP")
                   for x in remained_idp_bedtool_1] + \
                  [merge_idp_itx(fasta_file, dup_records[
                      int(x.name.split(",")[0].split("_")[-1])],
                                 [del_records[int(
                                     x.name.split(",")[1].split("_")[-1])]],
                                 int(x.fields[6]), x.fields[7], x.score, "IDP")
                   for x in remained_idp_bedtool_2] + \
                  [merge_idp_itx(fasta_file, dup_records[
                      int(x.name.split(",")[0].split("_")[-1])],
                                 [del_records[
                                      int(x.name.split(",")[1].split("_")[-1])],
                                  del_records[int(
                                      x.name.split(",")[2].split("_")[-1])]],
                                 int(x.fields[6]), x.fields[7], x.score, "ITX")
                   for x in itx_bedtool] + \
                  [merge_ctx(fasta_file, del_records[
                      int(x.name.split(",")[0].split("_")[-1])],
                             ins_records[
                                 int(x.name.split(",")[1].split("_")[-1])],
                             x.score) for x in ctx_bedtool]

    vcf_records = sorted(
        map(lambda x: remove_info_fields(x, ["SC_CHR2_STR"]), vcf_records),
        key=lambda x: (x.CHROM, x.POS))
    return vcf_records


def convert_metasv_bed_to_vcf(bedfile=None, vcf_out=None, workdir=None,
                              vcf_template_file=VCF_TEMPLATE, sample=None,
                              reference=None,
                              pass_calls=True):
    func_logger = logging.getLogger("%s" % convert_metasv_bed_to_vcf.__name__)
    if not os.path.exists(workdir):
        os.makedirs(workdir)

    intervals = []
    if bedfile:

        for interval in pybedtools.BedTool(bedfile):
            interval_info = get_interval_info(interval, pass_calls)
            if interval_info:
                updated_interval = pybedtools.Interval(interval.chrom,
                                                       interval_info["pos"],
                                                       interval_info["end"],
                                                       name="%s,%s,%d,%s" % (
                                                           base64.b64encode(
                                                               json.dumps(
                                                                   interval_info[
                                                                       "info"])),
                                                           interval_info[
                                                               "sv_type"],
                                                           interval_info[
                                                               "sv_length"],
                                                           ";".join(
                                                               interval_info[
                                                                   "svmethods"])),
                                                       score=interval.score,
                                                       otherfields=[
                                                           interval_info[
                                                               "genotype"]
                                                           , interval_info[
                                                               "sv_filter"]])
                if not intervals:
                    intervals.append(updated_interval)
                else:
                    merged_interval = check_duplicates(updated_interval,
                                                       intervals[-1])
                    if merged_interval:
                        func_logger.info("Merging intervals: %s and %s" % (
                        updated_interval, intervals[-1]))
                        intervals.pop()
                        intervals.append(merged_interval)
                    else:
                        intervals.append(updated_interval)
            else:
                func_logger.info("Skip interval: %s" % (interval))

    nonfilterd_bed = os.path.join(workdir, "final_nonfilterd.bed")
    filterd_bed = os.path.join(workdir, "final_filterd.bed")
    bedtool = pybedtools.BedTool(intervals).sort().moveto(nonfilterd_bed)
    filterd_bed = filter_confused_INS_calls(nonfilterd_bed, filterd_bed)

    vcf_template_reader = vcf.Reader(open(vcf_template_file, "r"))
    # The following are hacks to ensure sample name and contig names are put in the VCF header
    vcf_template_reader.samples = [sample]
    contigs = []
    fasta_file = None
    if reference:
        contigs = fasta_utils.get_contigs(reference)
        contigs_order_dict = {contig.name: index for (index, contig) in
                              enumerate(contigs)}
        vcf_template_reader.contigs = OrderedDict(
            [(contig.name, (contig.name, contig.length)) for contig in contigs])
        vcf_template_reader.metadata["reference"] = reference
        fasta_file = pysam.Fastafile(reference)

    vcf_template_reader.metadata["fileDate"] = str(datetime.date.today())
    vcf_template_reader.metadata["source"] = [" ".join(sys.argv)]
    vcf_writer = vcf.Writer(open(vcf_out, "w"), vcf_template_reader)
    vcf_records = []
    if filterd_bed:
        bedtool = pybedtools.BedTool(filterd_bed)
        for interval in bedtool:
            name_split = interval.name.split(",")
            info = json.loads(base64.b64decode(name_split[0]))
            # Fix info
            if "INSERTION_SEQUENCE" in info and not info["INSERTION_SEQUENCE"]:
                del info["INSERTION_SEQUENCE"]
            sv_type = name_split[1]
            sv_id = "."
            ref = fasta_file.fetch(str(interval.chrom), interval.start,
                                   interval.start + 1) if fasta_file else "."
            alt = [vcf.model._SV(sv_type)]
            qual = "."
            sv_filter = [interval.fields[7]]
            genotype = interval.fields[6]
            sv_format = "GT"
            sample_indexes = [0]
            vcf_record = vcf.model._Record(interval.chrom, interval.start,
                                           sv_id, ref, alt, qual,
                                           sv_filter, info, sv_format,
                                           sample_indexes)
            vcf_record.samples = vcf_template_reader._parse_samples([genotype],
                                                                    "GT",
                                                                    vcf_record)
            vcf_records.append(vcf_record)

    if contigs:
        vcf_records.sort(key=lambda x: (contigs_order_dict[x.CHROM], x.POS))
    else:
        vcf_records.sort(key=lambda x: (x.CHROM, x.POS))

    resolved_vcf_records = resolve_for_IDP_ITX_CTX(vcf_records, fasta_file)

    for vcf_record in resolved_vcf_records:
        vcf_writer.write_record(vcf_record)
    vcf_writer.close()

    func_logger.info("Tabix compressing and indexing %s" % vcf_out)
    pysam.tabix_index(vcf_out, force=True, preset="vcf")
