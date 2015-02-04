import logging
import sys
import os

import vcf

from sv_interval import SVInterval


logger = logging.getLogger(__name__)

mydir = os.path.dirname(os.path.realpath(__file__))
valid_svs = set(["DEL", "INS"])
tool_name = "BreakSeq"
source = set([tool_name])
sv_name_dict = {"Deletion": "DEL", "Insertion": "INS", "Inversion": "INV"}

class BreakSeqRecord:
    def __init__(self, record_string):
        self.name = tool_name
        fields = record_string.split()
        self.chromosome = fields[0]
        self.start = int(fields[3]) - 1
        self.end = int(fields[4])
        self.score = int(fields[5])

        info_fields = {f.split(" ")[0]: f.split(" ")[1] for f in fields[8].split(";")}
        self.filter = info_fields["QUAL"]
        self.abc_counts = zip(["A", "B", "C"], info_fields["ABC"].split(" ")[1].split(","))
        self.pe_count = int(info_fields["PE"])
        self.sv_type = sv_name_dict(fields[2])
        self.sv_len = 0 if self.sv_type == "INS" else (self.end - self.start)
        self.info = {
            "BS_CHR": self.chromosome,
            "BS_START": self.start,
            "BS_END": self.end,
            "BS_SCORE": self.score,
            "BS_FILTER": self.filter,
            "BS_ABC_COUNTS": ",".join(map(str, self.abc_counts)),
            "BS_PE_COUNT": self.pe_count
        }

    def __str__(self):
        return str(self.__dict__)

    def __repr__(self):
        return "<" + self.__class__.__name__ + " " + str(self.__dict__) + ">"

    def to_sv_interval(self):
        if self.sv_type not in valid_svs:
            return None

        if self.sv_type == "DEL":
            return SVInterval(self.chromosome,
                              self.start,
                              self.end,  # fudge
                              name=self.name,
                              sv_type=self.sv_type,
                              length=self.sv_len,
                              sources=source,
                              cipos=[],
                              info=self.info,
                              native_sv=self)
        elif self.sv_type == "INS":
            return SVInterval(self.chromosome,
                              self.start,
                              self.end,  # fudge
                              name=self.name,
                              sv_type=self.sv_type,
                              length=self.sv_len,
                              sources=source,
                              cipos=[],
                              info=self.info,
                              native_sv=self)
        else:
            logger.error("Bad SV type: " + repr(self))
            return None

    def to_vcf_record(self, sample):
        alt = ["<%s>" % self.sv_type]
        sv_len = -self.sv_len if self.sv_type == "DEL" else self.sv_len
        info = {"SVLEN": sv_len,
                "SVTYPE": self.sv_type}
        info["END"] = self.end

        info.update(self.info)

        vcf_record = vcf.model._Record(self.chromosome,
                                       self.start,
                                       ".",
                                       "N",
                                       alt,
                                       ".",
                                       ".",
                                       info,
                                       "GT",
                                       [0],
                                       [vcf.model._Call(None, sample, vcf.model.make_calldata_tuple("GT")(GT="1/1"))])
        return vcf_record


class BreakSeqReader:
    def __init__(self, file_name, reference_handle=None):
        logger.info("File is " + str(file_name))
        self.file_fd = open(file_name) if file_name is not None else sys.stdin
        self.reference_handle = reference_handle

    def __iter__(self):
        return self

    def next(self):
        while True:
            line = self.file_fd.next()
            if line[0] != "#":
                return BreakSeqRecord(line.strip())
            else:
                continue