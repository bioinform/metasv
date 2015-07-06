import logging
import sys

import vcf

from sv_interval import SVInterval


logger = logging.getLogger(__name__)

tool_name = "BreakSeq"
source = set([tool_name])
sv_name_dict = {"Deletion": "DEL", "Insertion": "INS", "Inversion": "INV"}


class BreakSeqRecord:
    def __init__(self, record_string):
        self.name = tool_name
        fields = record_string.split("\t")
        self.chromosome = fields[0]
        self.start = int(fields[3]) - 1
        self.end = int(fields[4])
        self.score = int(fields[5])

        info_fields = dict([f.split(" ") for f in fields[8].split(";")])
        self.filter = info_fields["QUAL"]
        self.abc_counts = zip(["A", "B", "C"], map(int, info_fields["ABC"].split(",")))
        self.pe_count = int(info_fields["PE"])
        self.sv_type = sv_name_dict[fields[2]]
        self.sv_len = 0 if self.sv_type == "INS" else (self.end - self.start)
        self.info = {
            "BS_START": self.start,
            "BS_END": self.end,
            "BS_SCORE": self.score,
            "BS_FILTER": self.filter,
            "BS_ABC_COUNTS": ",".join(map(lambda x: str(x[1]), self.abc_counts)),
            "BS_PE_COUNT": self.pe_count
        }

    def __str__(self):
        return str(self.__dict__)

    def __repr__(self):
        return "<" + self.__class__.__name__ + " " + str(self.__dict__) + ">"

    def to_sv_interval(self):
        if self.sv_type not in BreakSeqReader.svs_supported:
            return None

        return SVInterval(self.chromosome,
                          self.start,
                          self.end,
                          name=self.name,
                          sv_type=self.sv_type,
                          length=self.sv_len,
                          sources=source,
                          cipos=[],
                          info=self.info,
                          native_sv=self)

    def to_vcf_record(self, sample):
        alt = ["<%s>" % self.sv_type]
        sv_len = -self.sv_len if self.sv_type == "DEL" else self.sv_len
        info = {"SVLEN": sv_len, "SVTYPE": self.sv_type, "END": self.end}

        info.update(self.info)

        return vcf.model._Record(self.chromosome,
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


class BreakSeqReader:
    svs_supported = set(["DEL", "INS"])

    def __init__(self, file_name, reference_handle=None, svs_to_report=None):
        logger.info("File is " + str(file_name))
        self.file_fd = open(file_name) if file_name is not None else sys.stdin
        self.reference_handle = reference_handle
        self.svs_supported = BreakSeqReader.svs_supported
        if svs_to_report is not None:
            self.svs_supported &= set(svs_to_report)

    def __iter__(self):
        return self

    def next(self):
        while True:
            line = self.file_fd.next().strip()
            if line:
                if line[0] != "#":
                    record = BreakSeqRecord(line)
                    if record.sv_type in self.svs_supported:
                        return record
                else:
                    continue