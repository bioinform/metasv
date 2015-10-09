from svregion import *
import re

pattern = re.compile(r'^(.+)_(\d+)_(\d+)_(INS|DEL|INV|DUP)_\d+_NODE_\d+_length_(\d+)_cov_(\d*\.\d+|\d+)')


class SpadesContig:
    def __init__(self, name="", sequence=None):
        if sequence is not None:
            self.sequence = sequence
        self.parse_name(name)

    def parse_name(self, name):
        name_match = pattern.search(name)

        self.sv_region = SVRegion(name_match.group(1), int(name_match.group(2)), name_match.group(1),
                                  int(name_match.group(3)))
        self.sv_type = name_match.group(4)
        self.sv_len = 0
        self.sequence_len = int(name_match.group(5))
        self.covs = float(name_match.group(6))

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
