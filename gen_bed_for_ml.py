#!/net/kodiak/volumes/lake/shared/users/marghoob/my_env/bin/python

from __future__ import print_function
import argparse
import sys
import logging

FORMAT = '%(levelname)s %(asctime)-15s %(message)s'
logging.basicConfig(level=logging.INFO, format=FORMAT)
logger = logging.getLogger(__name__)

parser = argparse.ArgumentParser("Take all and false MetaSV VCF. Make BED with features out of them.\n",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--in_all_vcf", help="VCF file of all variants", required=True)
parser.add_argument("--in_false_vcf", help="VCF file of all variants", required=True)
parser.add_argument("--out_bed", help="Output BED file")

args = parser.parse_args()


class VcfRec:
    def __init__(self, line):
        ll = line.split("\t");
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
        ll = line.split('\t')
        vcf_list.append(VcfRec(line))

    return vcf_list

# read in the false VCF and store in a dictionary
false_vcf = read_all_vcf(args.in_false_vcf)
false_vcf_dict = {}  # this just stores some things that will make the VCF record unique (chr,loc,alt,SVLEN)
for line in false_vcf:
    false_vcf_dict["-".join([line.chrom, str(line.loc), line.alt, line.info["SVLEN"]])] = 1


# open the output file
outf = sys.stdout
if not args.out_bed == None:
    outf = open(args.out_bed, 'w')

# read in all VCF and flag the ones that are false
# at the same time generate the features
# and also output the BED file

methods = ["RP", "RD", "SR", "JM", "AS"]
features = ["CHROM", "START", "END", "LEN", "TRUE", "NUM_METHODS"] + methods

outf.write("#" + "\t".join(features) + "\n")

all_vcf = read_all_vcf(args.in_all_vcf)

for line in all_vcf:
    false_vec = "-".join([line.chrom, str(line.loc), line.alt, line.info["SVLEN"]])
    true_var = True
    if false_vec in false_vcf_dict:
        true_var = False

    method_ll = {i: 1 for i in line.info["SVMETHOD"].split(',')}

    out_vec = []

    # CHROM
    out_vec.append(line.chrom)

    #START
    out_vec.append(str(line.loc + 1))

    #END
    out_vec.append(str(line.loc + abs(int(line.info["SVLEN"]))))

    #LEN
    out_vec.append(str(abs(int(line.info["SVLEN"]))))

    #TRUE
    out_vec.append(str(true_var))

    #NUM_METHODS
    out_vec.append(str(len(method_ll)))

    #RP,RD,SR,JT,AS

    for m in methods:
        if m in method_ll:
            out_vec.append("1")
        else:
            out_vec.append("0")

    outf.write('\t'.join(out_vec) + '\n')

outf.close()




