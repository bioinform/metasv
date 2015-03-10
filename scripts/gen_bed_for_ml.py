#!/net/kodiak/volumes/lake/shared/users/marghoob/my_env/bin/python

from __future__ import print_function
import argparse
import sys
import logging

FORMAT = '%(levelname)s %(asctime)-15s %(message)s'
logging.basicConfig(level=logging.INFO, format=FORMAT)
logger = logging.getLogger(__name__)

parser = argparse.ArgumentParser("Take all and false MetaSV VCF (from VarSim). Make BED with features out of them.\n",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--in_all_vcf", help="VCF file of all variants", required=True)
parser.add_argument("--in_false_vcf", help="VCF file of false positive variants", required=True)
parser.add_argument("--in_unknown_vcf", help="VCF file of unknown status variants", required=True)
parser.add_argument("--out_bed", help="Output BED file")
parser.add_argument("--ignore_len", help="Ignore the length when matching false positives", action="store_true")


args = parser.parse_args()


class VcfRec:
    def __init__(self, line):
        ll = line.split("\t")
        self.chrom = ll[0]
        self.loc = int(ll[1])
        self.ref = ll[3]
        self.alt = ll[4]
        self.qual = ll[5]
        self.filter = ll[6]
        info = ll[7]
        form = ll[8]
        data = ll[9]
        self.info = {}
        self.data = {}
        self.line = line

        info_ll = info.split(";")
        for field in info_ll:
            fll = field.split("=")
            if len(fll) == 2:
                self.info[fll[0]] = fll[1]
            else:
                self.info[fll[0]] = True

        form_ll = form.split(":")
        data_ll = data.split(":")

        if len(form_ll) != len(data_ll):
            raise Exception("format not equal data", ll[0], ll[1], form, data)

        for i in xrange(0, len(form_ll)):
            self.data[form_ll[i]] = data_ll[i]


def read_all_vcf(vcf_file):
    f = open(vcf_file, 'r')
    vcf_list = []

    for line in f:
        line = line.strip()
        if len(line) == 0:
            continue
        if line[0] == "#":
            continue
        vcf_list.append(VcfRec(line))

    return vcf_list

def gen_false_tuple(var):
    if args.ignore_len:
        return "-".join([line.chrom, str(line.loc), line.alt])
    else:
        return "-".join([line.chrom, str(line.loc), line.alt, line.info["SVLEN"]])

# read in the false VCF and store in a dictionary
false_vcf = read_all_vcf(args.in_false_vcf)
false_vcf_dict = {}  # this just stores some things that will make the VCF record unique (chr,loc,alt,SVLEN)
for line in false_vcf:
    false_vcf_dict[gen_false_tuple(line)] = 1

# read in the unknown VCF and store in a dictionary
unknown_vcf = read_all_vcf(args.in_unknown_vcf)
unknown_vcf_dict = {}  # this just stores some things that will make the VCF record unique (chr,loc,alt,SVLEN)
for line in unknown_vcf:
    unknown_vcf_dict[gen_false_tuple(line)] = 1

# open the output file
outf = sys.stdout
if not args.out_bed is None:
    outf = open(args.out_bed, 'w')

"""
# Example lines

14      106880906       .       T       <DEL>   .       PASS    END=107175027;SVLEN=-242290;SVTYPE=DEL;
DP=58;SVTOOL=MetaSVMerge;SOURCES=14-106550910-106924720-373810-BreakDancer,14-106791601-106792400-799-CNVnator,
14-106805001-106814100-9099-CNVnator,14-106857201-106858000-799-CNVnator,14-106876601-106880600-3999-CNVnator,
14-106881006-107065514-184507-BreakSeq,14-106932636-107174927-242290-BreakSeq,14-106883201-106918800-35599-CNVnator,
14-106925577-106925749-172-BreakDancer,14-106932701-106935400-2699-CNVnator,14-106932820-107175108-242288-BreakDancer,
14-106986860-106986920-60-HaplotypeCaller,14-107174301-107174900-599-CNVnator;NUM_SVMETHODS=4;VT=SV;
SVMETHOD=AS,JM,RD,RP;AC=1;AF=0.500;AN=2;BD_CHR1=14;BD_CHR2=14;BD_ORI1=25+0-;BD_ORI2=0+25-;BD_POS1=106932820;
BD_POS2=107174928;BD_SCORE=99.0;BD_SUPPORTING_READ_PAIRS=25;BaseQRankSum=-0.380;CN_EVAL1=11.4296;
CN_EVAL2=8.9711e-12;CN_EVAL3=1.0;CN_EVAL4=1.0;CN_NORMALIZED_RD=0.311274;CN_Q0=0.0208333;FS=2.391;
MLEAC=1;MLEAF=0.500;MQ=59.74;MQ0=0;MQRankSum=0.181;QD=11.86;ReadPosRankSum=1.354;
VQSLOD=-1.060e-01;culprit=DP   GT      1/1

1       953160  .       C       <DEL>   .       PASS    END=953319;SVLEN=-158;SVTYPE=DEL;SVTOOL=MetaSVMerge;
SOURCES=1-953028-953326-165-BreakDancer,1-953161-953319-158-Pindel;NUM_SVMETHODS=2;VT=SV;SVMETHOD=RP,SR;
BD_CHR1=1;BD_CHR2=1;BD_ORI1=3+0-;BD_ORI2=0+3-;BD_POS1=953027;BD_POS2=953326;BD_SCORE=43.0;BD_SUPPORTING_READ_PAIRS=3;
PD_BP_RANGE_END=953344;PD_BP_RANGE_START=953161;PD_DOWN_READ_SUPP=5;PD_DOWN_UNIQ_READ_SUPP=5;PD_HOMLEN=25;
PD_HOMSEQ=;PD_NT_ADDED=;PD_NUM_NT_ADDED=0;PD_NUM_SAMPLE=1;PD_NUM_SAMPLE_SUPP=1;PD_NUM_SAMPLE_UNIQ_SUPP=1;
PD_READ_SUPP=8;PD_SIMPLE_SCORE=24;PD_SUM_MAPQ=440;PD_UNIQ_READ_SUPP=8;PD_UP_READ_SUPP=3;
PD_UP_UNIQ_READ_SUPP=3      GT      1/1

"""
# read in all VCF and flag the ones that are false
# at the same time generate the features
# and also output the BED file

# The order matters here
methods = ["RP", "RD", "SR", "JM", "AS", "SC"]
features = ["CHROM", "START", "END", "LEN", "TRUE", "NUM_METHODS"] + methods

# List of tuples, feature name and default value
details = [("CN_EVAL1", -1),
           ("CN_EVAL2", -1),
           ("CN_EVAL3", -1),
           ("CN_EVAL4", -1),
           ("CN_NORMALIZED_RD", -1),
           ("CN_Q0", -1),
           ("MQRankSum", -1),
           ("QD", -1),
           ("MQ0", -1),
           ("MQ", -1),
           ("ReadPosRankSum", -1),
           ("BD_SCORE", -1),
           ("BD_SUPPORTING_READ_PAIRS", -1),
           ("BaseQRankSum", -1),
           ("PD_DOWN_UNIQ_READ_SUPP", -1),
           ("PD_UP_UNIQ_READ_SUPP", -1),
           ("PD_UNIQ_READ_SUPP", -1),
           ("PD_SIMPLE_SCORE", -1),
           ("PD_SUM_MAPQ", -1),
           ("PD_HOMLEN", -1),
           ("AA_UNIQ_COV", -1),
           ("AA_TOTAL_COV", -1),
           ("AA_TOTAL_STRAND", -1),
           ("AA_PROP_REPEAT", -1),
           ("AA_PROP_ALIGNED", -1),
           ("AA_DISCORDANT_HIGH", -1),
           ("AA_DISCORDANT_LOW", -1),
           ("AA_END_PROP_ALIGNED", -1),
           ("BA_FLANK_PERCENT", -1),
           ("BA_NFRAGS", -1),
           ("BA_NUM_ALT", -1),
           ("BA_NUM_BP", -1),
           ("BA_NUM_GOOD_REC", -1),
           ("BA_PERCENT_MATCH", -1),
           ("BA_BP_SCORE", 10000)

]

outf.write("#" + "\t".join(features))

for d in details:
    outf.write('\t' + d[0])
outf.write('\n')

all_vcf = read_all_vcf(args.in_all_vcf)

for line in all_vcf:
    false_vec = gen_false_tuple(line)

    if false_vec in unknown_vcf_dict:
        # ignore the unknown state records for analysis
        continue

    true_var = True
    if false_vec in false_vcf_dict:
        true_var = False

    method_ll = {i: 1 for i in line.info["SVMETHOD"].split(',')}

    out_vec = []

    # CHROM
    out_vec.append(line.chrom)

    # START
    out_vec.append(str(line.loc + 1))

    # END
    out_vec.append(str(line.loc + abs(int(line.info["SVLEN"]))))

    # LEN
    out_vec.append(str(abs(int(line.info["SVLEN"]))))

    # TRUE
    out_vec.append(str(true_var))

    # NUM_METHODS
    out_vec.append(str(len(method_ll)))

    #RP,RD,SR,JT,AS

    for m in methods:
        if m in method_ll:
            out_vec.append("1")
        else:
            out_vec.append("0")

    for d in details:
        if d[0] in line.info:
            out_vec.append(str(line.info[d[0]]))
        else:
            out_vec.append(str(d[1]))

    outf.write('\t'.join(out_vec) + '\n')

outf.close()




