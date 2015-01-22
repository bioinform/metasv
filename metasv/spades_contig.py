from svregion import *


class SpadesContig:
    def __init__(self, name="", sequence=None):
        if sequence is not None:
            self.sequence = sequence
        self.parse_name(name)

    def parse_name(self, name):
        name_fields = name.strip().split("_")

        self.sv_region = SVRegion(name_fields[0], int(name_fields[1]), name_fields[0], int(name_fields[2]))
        self.sv_type = name_fields[3]
        self.sv_len = 0
        self.sequence_len = int(name_fields[8])
        self.covs = float(name_fields[10])

        self.raw_name = name.strip()

    def update_sequence(self, sequence):
        self.sequence = sequence.upper()
        self.sequence_len = len(self.sequence)

    def __str__(self):
        return str(self.__dict__)

    def __repr__(self):
        return str(self)

    def __cmp__(self, other):
        if self.covs != other.covs:
            return cmp(self.covs, other.covs)

        return cmp(self.raw_name, other.raw_name)


