#!/usr/bin/env python2.7

import sys
import os
import argparse
import subprocess
import pysam

parser = argparse.ArgumentParser("Compute zygosity for SV tools output", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--bams", metavar="bams", nargs="+", help="BAMs", required=True)
parser.add_argument("--input", metavar="input", help="Breakdancer out file or pindel out prefix", required=True)
parser.add_argument("--tool", metavar="tool", help="SV tool name", required=True, choices = ["pindel", "breakdancer"])
parser.add_argument("--mean", metavar="mean", help="Insert size mean", required=True, type=float)
parser.add_argument("--sd", metavar="sd", help="Standard deviation", required=True, type=float)
parser.add_argument("--dg", action = "store_true", help = "Disable genotyping")

args = parser.parse_args()

svtypes = ["D", "LI", "SI", "TD", "INV"]

if args.tool == "pindel":
  input_handles = [open(args.input + "_" + svtype) for svtype in svtypes]
  tool = "Pindel"
elif args.tool == "breakdancer":
  input_handles = [sys.stdin if args.input == "-" else open(args.input)]
  tool = "Breakdancer"

bam_handles = []
for bam in args.bams:
  bam_handles.append(pysam.Samfile(bam, "rb"))

def get_num_reads_supporting_ref(bam_handles, min_isize, max_isize, chromosome, start, end):
  total_normal_reads = 0
  total_read_bases = 0
  for bam_handle in bam_handles:
    itr = bam_handle.fetch(chromosome, start, end)
    for aln in itr:
      if aln.is_duplicate: continue
      if not aln.is_paired: continue
      if aln.is_unmapped or aln.mate_is_unmapped: continue
      if aln.rnext != aln.tid: continue
      if aln.is_reverse:
        if not (aln.pnext < aln.pos and not aln.mate_is_reverse): continue
      else:
        if not (aln.pnext > aln.pos and aln.mate_is_reverse): continue
      tlen = aln.tlen if aln.tlen > 0 else (-aln.tlen)
      if tlen >= min_isize and tlen <= max_isize:
        total_normal_reads = total_normal_reads + 1
        total_read_bases = total_read_bases + aln.qlen
  return total_normal_reads, total_read_bases

pindel_sv_type_dict = {"D": "DEL", "I": "INS", "TD": "DUP:TANDEM", "LI": "LI", "INV": "INV"};

def parse_line(line, tool):
  line_stripped = line.strip()
  if tool == "pindel" and line.find("ChrID") == -1: return None
  if tool == "breakdancer" and line.startswith("#"): return None

  fields = line.split()
  if tool == "pindel":
    sv_type = pindel_sv_type_dict[fields[1]]
    if sv_type != "LI":
      size, chr1, pos1, pos2, num_reads, score = int(fields[2]), fields[7], int(fields[9]) - 1, int(fields[10]) - 1, int(fields[16]), int(fields[24])
    else:
      chr1, pos1, pos2, num_reads = fields[3], int(fields[4]) - 1, int(fields[7]) - 1, int(fields[6]) + int(fields[9])
      size = 0
      pos1 = min(pos1, pos2)
      pos2 = pos1 + 1
      score = 99
      sv_type = "INS"
  if tool == "breakdancer":
    chr1, pos1, orientation1, chr2, pos2, orientation2, sv_type, size, score, num_reads = fields[0:10]
    pos1, pos2, size, num_reads = int(pos1), int(pos2), abs(int(size)), int(num_reads)
    score = float(score)
    if sv_type == "ITX" or sv_type == "CTX": return None
  return sv_type, chr1, pos1, pos2, size, num_reads, score

min_isize = max(args.mean - 3 * args.sd, 0)
max_isize = args.mean + 3 * args.sd

for input_handle in input_handles:
  for line in input_handle.readlines():
    record = parse_line(line, args.tool)
    if not record: continue
    sv_type, chr1, pos1, pos2, size, num_reads, score = record

    pos2 = max(pos2, pos1+1)
    if args.dg:
    	normal_read_count = 1
    	normal_read_bases = 1
    	normal_coverage = 1
    	gt = "1/1"
    else:
    	normal_read_count, normal_read_bases = get_num_reads_supporting_ref(bam_handles, min_isize, max_isize, chr1, pos1, pos2)
    	normal_coverage = float(normal_read_bases) / max(1, pos2 - pos1)
    	gt = "0/1" if normal_coverage / float(max(1, num_reads)) > 0.2 else "1/1"
    print "%s\t%s\t%d\t%d\t%d\t%d\t%d\t%g\t%s\t%s" % (sv_type, chr1, pos1+1, pos2, size, normal_read_count, num_reads, score, gt, tool)
    sys.stdout.flush()
  input_handle.close()

for bam_handle in bam_handles:
  bam_handle.close()
