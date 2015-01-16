#!/usr/bin/env python2.7

import sys
import os
import argparse
import subprocess
import pysam

parser = argparse.ArgumentParser("Convert genotyped BreakDancer output to VCF", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--sv_file", metavar="sv_file", help="SV file", required=False, default="-")
parser.add_argument("--reference", metavar="reference", help="Reference file", required=True)
parser.add_argument("--sample", metavar="Sample", help="Sample name", required=True)
parser.add_argument("--sort", action="store_true", help="Sort the input")

args = parser.parse_args()

input_handle = sys.stdin if args.sv_file == "-" else open(args.sv_file)

fasta_handle = pysam.Fastafile(args.reference)

def get_contigs(fai_filename):
  fai_file = open(fai_filename)
  contigs = {}
  contigs_order = {}
  linenum = 0
  for line in fai_file.readlines():
    line = line.strip()
    line_items = line.split("\t")
    name, length = line_items[0:2]
    name = name.split(" ")[0]
    contigs[name] = int(length)
    contigs_order[name] = linenum
    linenum += 1
  fai_file.close()
  return contigs, contigs_order

def line_to_tuple(line):
  line = line.strip()
  fields = line.split("\t")
  return tuple(fields[0:2]) + tuple([int(i) for i in fields[2:7]]) + (float(fields[7]),) + tuple(fields[8:10])

contigs, contigs_order = get_contigs(args.reference + ".fai")
contig_names = contigs.keys()
contig_names.sort(key = lambda tup: contigs_order[tup])

contig_str = ""
for contig_name in contig_names:
  contig_str += "##contig=<ID=%s,length=%d>\n" % (contig_name, contigs[contig_name])

sys.stdout.write("""##fileformat=VCFv4.1
##reference=%s
##source=%s
##INFO=<ID=BKPTID,Number=.,Type=String,Description=\"ID of the assembled alternate allele in the assembly file\">
##INFO=<ID=CIEND,Number=2,Type=Integer,Description=\"Confidence interval around END for imprecise variants\">
##INFO=<ID=CIPOS,Number=2,Type=Integer,Description=\"Confidence interval around POS for imprecise variants\">
##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">
##INFO=<ID=HOMLEN,Number=.,Type=Integer,Description=\"Length of base pair identical micro-homology at event breakpoints\">
##INFO=<ID=HOMSEQ,Number=.,Type=String,Description=\"Sequence of base pair identical micro-homology at event breakpoints\">
##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Imprecise structural variation\">
##INFO=<ID=MEINFO,Number=4,Type=String,Description=\"Mobile element info of the form NAME,START,END,POLARITY\">
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">
##INFO=<ID=NORMAL_COUNT,Number=1,Type=Integer,Description=\"Number of normal reads supporting reference\">
##INFO=<ID=SCORE,Number=1,Type=Integer,Description=\"Score from tool\">
##INFO=<ID=NUM_READS,Number=1,Type=Integer,Description=\"Number of reads supporting event\">
##INFO=<ID=TOOLNAME,Number=1,Type=String,Description=\"Tool used to generate the SV\">
##FILTER=<ID=LowQual,Description=\"Low Quality\">
%s##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
##ALT=<ID=DEL,Description=\"Deletion\">
##ALT=<ID=DEL:ME:ALU,Description=\"Deletion of ALU element\">
##ALT=<ID=DEL:ME:L1,Description=\"Deletion of L1 element\">
##ALT=<ID=DUP,Description=\"Duplication\">
##ALT=<ID=DUP:TANDEM,Description=\"Tandem Duplication\">
##ALT=<ID=INS,Description=\"Insertion of novel sequence\">
##ALT=<ID=INS:ME:ALU,Description=\"Insertion of ALU element\">
##ALT=<ID=INS:ME:L1,Description=\"Insertion of L1 element\">
##ALT=<ID=INV,Description=\"Inversion\">
##ALT=<ID=CNV,Description=\"Copy number variable region\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n""" % (args.reference, sys.argv[0], contig_str, args.sample))

records = map(line_to_tuple, input_handle.readlines())
if args.sort:
  records.sort(key = lambda tup: (contigs_order[tup[1]], tup[2], tup[3], tup[4]))

for record in records:
  sv_type, chr1, pos1, pos2, size, normal_read_count, num_reads, score, gt, tool = record
  alt_allele = "<%s>" % (sv_type)

  ref_allele = fasta_handle.fetch(chr1, pos1-1, pos1)

  if sv_type in ["DEL"]: size = -size
  info = "TOOLNAME=%s;SVLEN=%d;SVTYPE=%s;END=%d;IMPRECISE;NORMAL_COUNT=%d;NUM_READS=%d;SCORE=%g" % (tool, size, sv_type, pos2, normal_read_count, num_reads, score)

  sys.stdout.write("%s\t%d\t.\t%s\t%s\t.\tPASS\t%s\tGT\t%s\n" % (chr1, pos1, ref_allele, alt_allele, info, gt))

fasta_handle.close()
