import os

class Contig:
    def __init__(self, name, length):
        self.name = name
        self.length = length


def get_contigs(reference):
    contigs = []
    reference_fai = str(reference) + ".fai"
    if reference is not None and os.path.isfile(reference_fai):
        with open(reference_fai) as fai_file:
            for line in fai_file:
                line_items = line.strip().split("\t")
                name, length = line_items[0:2]
                name = name.split(" ")[0]
                contigs.append(Contig(name, int(length)))
    return contigs
