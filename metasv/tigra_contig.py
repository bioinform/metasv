from svregion import *


class TigraContig:
    def __init__(self, name="", sequence=None):
        if sequence is not None:
            self.sequence = sequence
        self.parse_name(name)

    def parse_name(self, name):
        name_fields = name.strip().split()
        name_split = name_fields[0].split(".")
        self.sequence_len = int(name_fields[1])
        self.covs = float(name_fields[2])
        self.raw_name = name_fields[0]

        if len(name_fields) == 4:
            self.contig_kmerutil = int(name_fields[3])
        else:
            if name_fields[4][1:]:
                self.contig_i = map(lambda y: (y[0], int(y[1])),
                                    map(lambda x: x.split(":"), name_fields[4][1:].split(",")[:-1]))

            if name_fields[5][1:]:
                self.contig_o = map(lambda y: (y[0], int(y[1])),
                                    map(lambda x: x.split(":"), name_fields[5][1:].split(",")[:-1]))

            self.tags = int(name_fields[6])
            self.contig_kmerutil = int(name_fields[7])
            self.left_branches_num = int(name_fields[3][0:1])  # should equal len(self.contig_i) if contig_kmerutil > 0
            self.right_branches_num = int(name_fields[3][1:2])  # should equal len(self.contig_o) if contig_kmerutil > 0

        self.sv_region = SVRegion(name_split[0], int(name_split[1]), name_split[2], int(name_split[3]))
        self.sv_type = name_split[4]
        self.sv_len = int(name_split[5])

        if name_split[6] != "+-":
            self.plus_minus = name_split[6]
        self.contig_num = name_split[7]

        self.path = map(int, name_split[8:]) if name_split[8:] else []

    def update_sequence(self, sequence):
        self.sequence = sequence.upper()
        self.sequence_len = len(self.sequence)

    def __str__(self):
        return str(self.__dict__)

    def __cmp__(self, other):
        if self.covs != other.covs:
            return cmp(self.covs, other.covs)

        return cmp(self.raw_name, other.raw_name)


