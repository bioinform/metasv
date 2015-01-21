import sys
import os
import argparse


class Contig:
    def __init__(self, name, length):
        self.name = name
        self.length = length


def get_contigs(reference):
    contigs = []
    with open(reference + ".fai") as fai_file:
        for line in fai_file:
            line_items = line.strip().split("\t")
            name, length = line_items[0:2]
            name = name.split(" ")[0]
            contigs.append(Contig(name, int(length)))
    return contigs
