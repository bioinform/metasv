import logging
import pysam
from sv_interval import SVInterval

logger = logging.getLogger(__name__)

pindel_source = set(["Pindel"])
min_coverage = 10
het_cutoff = 0.2
hom_cutoff = 0.8

GT_REF = "0/0"
GT_HET = "0/1"
GT_HOM = "1/1"

PINDEL_TO_SV_TYPE = {"I": "INS", "D": "DEL", "LI": "INS", "TD": "DUP:TANDEM", "INV": "INV"}


class PindelRecord:
    def __init__(self, record_string, reference_handle=None):
        fields = record_string.split()
        self.sv_type = fields[1]

        if self.sv_type != "LI":
            self.sv_len = int(fields[2])
            self.num_nt_added = map(int, fields[4].split(":"))
            self.nt_added = map(lambda x: x.replace('"', ''), fields[5].split(":"))
            self.chromosome = fields[7]
            self.start_pos = int(fields[9])
            self.end_pos = int(fields[10])
            self.bp_range = (int(fields[12]), int(fields[13]))
            self.homlen = self.bp_range[1] - self.end_pos
            self.homseq = reference_handle.fetch(self.chromosome, self.end_pos - 1,
                                                 self.bp_range[1] - 1) if reference_handle else ""
            self.samples = [{"name": fields[i], "ref_support_at_start": int(fields[i + 1]),
                             "ref_support_at_end": int(fields[i + 2]), "plus_support": int(fields[i + 3]),
                             "minus_support": int(fields[i + 4])} for i in xrange(31, len(fields), 5)]
        else:
            self.sv_len = 0
            self.chromosome = fields[3]
            self.start_pos = int(fields[4])
            self.end_pos = int(fields[7])
            self.bp_range = (self.start_pos, self.end_pos)
            self.homlen = 0
            self.homseq = ""
            self.samples = [{"name": fields[i], "plus_support": int(fields[i + 2]), "minus_support": int(fields[i + 4])}
                            for i in xrange(10, len(fields), 5)]

        self.derive_genotype()

    def derive_genotype(self):
        if self.sv_type == "LI" or self.sv_type == "I":
            self.gt = GT_HET if sum(
                [sample["plus_support"] + sample["minus_support"] for sample in self.samples]) > 0 else GT_REF
            return

        total_event_reads = sum([sample["plus_support"] + sample["minus_support"] for sample in self.samples])
        total_ref_reads = sum(
            [sample["ref_support_at_start"] + sample["ref_support_at_end"] for sample in self.samples])
        if total_event_reads + total_ref_reads < min_coverage:
            self.gt = GT_REF
            return

        allele_fraction = float(total_event_reads) / (float(total_event_reads) + float(total_ref_reads))
        if allele_fraction < het_cutoff:
            self.gt = GT_REF
        elif allele_fraction < hom_cutoff:
            self.gt = GT_HET
        else:
            self.gt = GT_HOM

    def to_sv_interval(self):
        sv_type = PINDEL_TO_SV_TYPE[self.sv_type]
        if sv_type != "INS":
            return SVInterval(self.chromosome, self.start_pos, self.end_pos, "Pindel", sv_type=sv_type,
                              length=self.sv_len, sources=pindel_source, native_sv=self)
        else:
            return SVInterval(self.chromosome, self.start_pos, self.start_pos, "Pindel", sv_type=sv_type,
                              length=self.sv_len, sources=pindel_source, native_sv=self, wiggle=100, gt=self.gt)

    def __str__(self):
        return str(self.__dict__)


class PindelReader:
    def __init__(self, file_name, reference_handle=None):
        logger.info("File is " + file_name)
        self.file_fd = open(file_name)
        self.reference_handle = reference_handle

    def __iter__(self):
        return self

    def next(self):
        while True:
            line = self.file_fd.next()
            if line.find("ChrID") >= 1:
                return PindelRecord(line.strip(), self.reference_handle)
