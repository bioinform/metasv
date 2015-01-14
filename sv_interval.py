#!/net/kodiak/volumes/lake/shared/users/marghoob/my_env/bin/python

import logging

logger = logging.getLogger(__name__)

import copy
import bisect
import pysam
import pybedtools

svs_of_interest = ["DEL", "INS", "DUP", "DUP:TANDEM", "INV"]
sv_sources = ["Pindel", "BreakSeq", "HaplotypeCaller", "BreakDancer", "CNVnator"]
precise_sv_sources = ["Pindel", "BreakSeq", "HaplotypeCaller"]
sv_sources_to_type = {"Pindel": "SR", "BreakSeq": "JM", "BreakDancer": "RP", "CNVnator": "RD", "HaplotypeCaller": "AS"}

class SVInterval:
  def __init__(self, chrom=None, start=0, end=0, name=None, sv_type=None, length=0, sources=set(), gt="./1", wiggle=0, info=None, cipos=[], ciend=[], native_sv=None):
    self.chrom = chrom
    self.start = start
    self.end = end
    self.name = name
    self.sv_type = sv_type
    self.length = length
    self.info = info
    self.sources = sources
    self.gt = gt
    self.wiggle = wiggle
    self.sub_intervals = []
    self.is_precise = False
    self.is_validated = False
    self.validating_interval = None
    self.info_dict = {}
    self.cipos = cipos
    self.ciend = ciend
    self.native_sv = native_sv

  def merge(self, interval):
    self.start = min(self.start, interval.start - interval.wiggle)
    self.end = max(self.end, interval.end + interval.wiggle)
    self.length = max(self.length, interval.length)
    self.name = self.name + "," + interval.name
    self.sub_intervals.append(interval)
    self.sources |= interval.sources

  def set_merged(self, interval1, interval2):
    self.chrom = interval1.chrom
    self.start = min(interval1.start - interval1.wiggle, interval2.start - interval2.wiggle)
    self.end = max(interval1.end + interval1.wiggle, interval2.end + interval2.wiggle)
    self.length = max(interval1.length, interval2.length)
    self.name = interval1.name + "," + interval2.name
    self.sv_type = interval1.sv_type
    self.info = "merged"
    self.sub_intervals = [interval1, interval2]
    self.sources = interval1.sources | interval2.sources
    self.gt = interval1.gt

  def __lt__(self, other):
    if self.chrom != other.chrom: return self.chrom < other.chrom
    if self.start != other.start: return self.start < other.start
    return self.end < other.end

  def __str__(self):
    if self.sub_intervals:
      return ",".join([str(interval) for interval in self.sub_intervals])
    return "%s-%d-%d-%d-%s" % (self.chrom, self.start, self.end, self.length, ",".join(list(self.sources)))

  def __repr__(self):
    return "<" + self.__class__.__name__ + " " + str(self.__dict__) + ">"

  def overlaps(self, other, min_fraction_self=1e-9, min_fraction_other=1e-9, min_overlap_length_self=1, min_overlap_length_other=1):
    if self.chrom != other.chrom: return False
    if max(self.start - self.wiggle, other.start - other.wiggle) >= min(self.end + self.wiggle, other.end + other.wiggle): return False

    self_length = float(self.end - self.start + 2 * self.wiggle)
    other_length = float(other.end - other.start + 2 * other.wiggle)
    overlap_length = min(self.end + self.wiggle, other.end + other.wiggle) - max(self.start - self.wiggle, other.start - other.wiggle)
    return float(overlap_length) >= max(min_fraction_self * self_length, min_fraction_other * other_length) and overlap_length >= max(min_fraction_self, min_overlap_length_other)

  def is_adjacent(self, other, gap=0):
    if self.chrom != other.chrom: return False
    return (self.end + gap >= other.start and self.end + gap < other.end) or (other.end + gap >= self.start and other.end + gap < self.end)

  def get_start(self):
    if not self.sub_intervals:
      return self.start

    return min([interval.get_start() for interval in self.sub_intervals])

  def get_end(self):
    if not self.sub_intervals:
      return self.end

    return max([interval.get_end() for interval in self.sub_intervals])

  def do_validation(self, overlap_ratio = 0.5):
    self.start = self.get_start()
    self.end = self.get_end()

    if self.sv_type not in svs_of_interest:
      return

    if not self.sub_intervals:
      self.is_precise = list(self.sources)[0] in precise_sv_sources
      #if self.sv_type == "INV" and self.length >= 1000: self.is_validated = True
      return

    lists = {source: [] for source in sv_sources}

    precise_merged = merge_intervals([interval for interval in self.sub_intervals if list(interval.sources)[0] in precise_sv_sources])

    for interval in self.sub_intervals:
      #logger.debug("Examining subinterval " + str(interval))
      lists[list(interval.sources)[0]].append(interval)

    for source in sv_sources:
      if len(lists[source]) == 1:
        if self.overlaps(lists[source][0], overlap_ratio, overlap_ratio) \
                or interval_overlaps_interval_list(lists[source][0], precise_merged, overlap_ratio, overlap_ratio):
          self.is_validated = len(self.sources) > 1
          self.validating_interval = lists[source][0]
          self.is_precise = source in precise_sv_sources
          break

    if not self.is_validated:
      return

    if self.is_precise and self.validating_interval is not None:
      self.start = self.validating_interval.start
      self.end = self.validating_interval.end
      self.length = self.validating_interval.length
      self.gt = self.validating_interval.gt

      if self.length == 0:
        for source in sv_sources:
          if source in self.sources:
            if lists[source] and lists[source][0].length > 0:
              self.length = lists[source][0].length
              break

  def fix_pos(self):
    if self.sv_type == "INS":
      mid = (self.start + self.end) / 2
      if not self.is_precise:
        self.end = self.start
        self.cipos = ["0", "%d" % (self.end - self.start)]
        #self.cipos = ["%d" % (mid - self.start), "%d" % (mid - self.end)]
      else:
        self.start = mid
        self.end = mid

  def to_vcf_record(self, fasta_handle):
    if self.start <= 0: return None
    if self.sv_type not in svs_of_interest: return None
    if not self.sub_intervals and list(self.sources)[0] == "HaplotypeCaller": return None
    if len(self.sources) == 1 and list(self.sources)[0] == "HaplotypeCaller": return None

    #logger.debug("Converting interval %s to VCF record" % (str(self)))

    vcf_record_ref = fasta_handle.fetch(self.chrom, self.start - 1, self.start)
    vcf_record_alt = "<%s>" % (self.sv_type)
    vcf_record_filter = "PASS" if self.is_validated else "LowQual"

    vcf_record_info = "END=%d;SVLEN=%d;SVTYPE=%s;VT=SV;SVTOOL=%s;NUM_SVMETHODS=%d;SOURCES=%s" % (self.end, -self.length if self.sv_type == "DEL" else self.length, self.sv_type, "MetaSVMerge", len(self.sources), str(self))
    if self.cipos: vcf_record_info += ";CIPOS=%s" % (",".join(self.cipos))
    if not self.is_precise: vcf_record_info += ";IMPRECISE"
    svmethods = [sv_sources_to_type[tool] for tool in self.sources]
    svmethods.sort()
    vcf_record_info += ";SVMETHOD=" + ",".join(svmethods)
    return "%s\t%d\t.\t%s\t%s\t.\t%s\t%s\tGT\t%s" % (self.chrom, self.start, vcf_record_ref, vcf_record_alt, vcf_record_filter, vcf_record_info, self.gt)

    return vcf_record   

  def to_bed_interval(self, sample_name):
    if self.start <= 0: return None
    if self.sv_type not in ["DEL", "INS", "INV"]: return None
    #if not self.sub_intervals and list(self.sources)[0] == "HaplotypeCaller": return None
    #if len(self.sources) == 1 and list(self.sources)[0] == "HaplotypeCaller": return None
    end = self.end if self.sv_type != "INS" else (self.end + 1)

    return pybedtools.Interval(self.chrom, self.start, end, name = "%s,%d,%s" % (self.sv_type, self.length, ";".join(sorted([sv_sources_to_type[tool] for tool in self.sources]))), score = str(len(self.sources)))

  def to_svp_record(self, sample_name, id_num):
    if self.start <= 0: return None
    if self.sv_type not in ["DEL", "INS", "INV"]: return None
    if not self.sub_intervals and list(self.sources)[0] == "HaplotypeCaller": return None
    if len(self.sources) == 1 and list(self.sources)[0] == "HaplotypeCaller": return None

    start_outer = self.start
    start_inner = self.start
    end_inner = self.end
    end_outer = self.end
    type_of_computational_approach = ",".join(sorted([sv_sources_to_type[tool] for tool in self.sources]))

    return "%s\t%d\t%d\t%d\t%d\t%s\t%d\t%s\t%s\t%s\t%s\tMetaSV\t%d" % (self.chrom, start_outer, start_inner, end_inner, end_outer, self.sv_type, self.length, "BWA", "Illumina", sample_name, type_of_computational_approach, id_num)

def interval_overlaps_interval_list(interval, interval_list, min_fraction_self=1e-9, min_fraction_other=1e-9):
    index = bisect.bisect_left(interval_list, interval)
    if index > 0 and interval.overlaps(interval_list[index-1], min_fraction_self, min_fraction_other, min_overlap_length_self=1, min_overlap_length_other=1): return True
    if index < len(interval_list) and interval.overlaps(interval_list[index], min_fraction_self, min_fraction_other, min_overlap_length_self=1, min_overlap_length_other=1): return True
    return False

def merge_intervals(interval_list):
  interval_list.sort()
  merged_intervals = []
  if not interval_list: return []

  current_merged_interval = copy.deepcopy(interval_list[0])
  for i in xrange(len(interval_list)-1):
    next_interval = interval_list[i+1]
    if current_merged_interval.overlaps(next_interval) or current_merged_interval.is_adjacent(next_interval):
      if current_merged_interval.sub_intervals:
        current_merged_interval.merge(next_interval)
      else:
        new_merged_interval = SVInterval()
        #logger.debug("Merging %s with %s" % (repr(current_merged_interval), repr(next_interval)))
        new_merged_interval.set_merged(current_merged_interval, next_interval)
        current_merged_interval = new_merged_interval
    else:
      merged_intervals.append(current_merged_interval)
      current_merged_interval = copy.deepcopy(next_interval)
  merged_intervals.append(current_merged_interval)
  merged_intervals.sort()
  return merged_intervals
