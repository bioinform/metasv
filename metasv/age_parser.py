import logging
import itertools

logger = logging.getLogger(__name__)


def get_unique_age_records(age_records):
    unique_records = []

    for a1 in age_records:
        if a1.duplicate: continue
        unique_records.append(a1)

        for a2 in age_records:
            if a2.duplicate: continue
            if a1.contig.raw_name == a2.contig.raw_name: continue

            if a1.start1_end1s == a2.start1_end1s:
                a2.duplicate = True

    return unique_records


class AgeInput:
    def __init__(self, fname, length):
        self.fname = fname
        self.length = length

    def __str__(self):
        return "<%s(%s, %d)>" % (self.__class__.__name__, self.fname, self.length)

    def __repr__(self):
        return str(self)


class AgeRecord:
    def __init__(self, age_out_file=None):
        self.aligned_bases = 0
        self.breakpoint_identities = []
        self.inputs = []  # inputs[1] is expected to have an insertion w.r.t. inputs[0]
        self.percent = 100
        self.percents = [100, 100]
        self.start1_end1s = []
        self.start2_end2s = []
        self.excised_regions = []
        self.n_alt = 0
        self.alternate_regions = []

        self.flanking_regions = [0, 0]
        self.ref_flanking_regions = [0, 0]
        self.hom = 0
        self.flank_percent = 0
        self.nfrags = 1

        self.cs = 0
        self.ce = 0

        self.very_bad_assembly = False
        self.bad_assembly = False

        self.start = -1
        self.end = -1
        self.excised_lengths = [0, 0]

        self.assembly_contig = None
        self.used = False
        self.duplicate = False

        self.age_file = age_out_file
        self.score = 0

        if age_out_file is None:
            return

        with open(age_out_file) as age_fd:
            while True:
                line = age_fd.readline()
                if not line: break

                line = line.strip()

                if line.startswith("Aligned:"):
                    self.aligned_bases = int(line.split()[1])
                    continue

                if line.startswith("First seq"):
                    words = line.split()
                    file1 = words[-1]
                    len1 = int(words[-3])
                    self.inputs.append(AgeInput(file1, len1))
                    continue

                if line.startswith("Second seq"):
                    words = line.split()
                    file2 = words[-1]
                    len2 = int(words[-3])
                    self.inputs.append(AgeInput(file2, len2))
                    continue

                if line.startswith("Identic:"):
                    nums = map(int, line.split()[2::2])
                    self.percent = nums[0]
                    if len(nums) > 1:
                        self.percents = nums[1:3]
                    continue

                if line.startswith("Score:"):
                    self.score = int(line.split()[1])
                    continue

                if line.startswith("Alignment:"):
                    self.start1_end1s = map(lambda y: map(int, y),
                                            map(lambda x: x.split(","), line.split(":")[1].split()))
                    self.start2_end2s = map(lambda y: map(int, y),
                                            map(lambda x: x.split(","), line.split(":")[2].split()))

                    self.nfrags = len(self.start1_end1s)
                    self.flanking_regions[0] = abs(self.start2_end2s[0][1] - self.start2_end2s[0][0] + 1)
                    self.flanking_regions[1] = 0 if len(self.start2_end2s) == 1 else abs(
                        self.start2_end2s[1][1] - self.start2_end2s[1][0] + 1)

                    self.ref_flanking_regions[0] = abs(self.start1_end1s[0][1] - self.start1_end1s[0][0] + 1)
                    self.ref_flanking_regions[1] = 0 if len(self.start1_end1s) == 1 else abs(
                        self.start1_end1s[1][1] - self.start1_end1s[1][0] + 1)
                    if len(self.start2_end2s) > 1:
                        cs, ce = self.start1_end1s[0][1], self.start2_end2s[1][0]
                    continue

                if line.startswith("EXCISED REGION(S):"):
                    self.excised_regions = map(lambda y: map(int, y),
                                               map(lambda x: x.split(","), line.split(":")[1].split()))
                    continue

                if line.startswith("ALTERNATIVE REGION(S):"):
                    words = line.split(":")[1].split()
                    self.n_alt = int(words[0])
                    self.alternate_regions = map(lambda y: map(int, y), map(lambda x: x.split(","), words[1:]))
                    continue

                if line.startswith("Identity at b"):
                    self.breakpoint_identities = map(lambda y: map(int, y),
                                                     map(lambda x: x.split(","), line.split(":")[1].split()))
                    self.hom = max(0, self.breakpoint_identities[1][0])
                    continue

        if len(self.inputs) == 0:
            self.flank_percent = 0
            logger.warn("%s has problems" % (age_out_file))
        else:
            self.flank_percent = int(round(100.0 * sum(self.flanking_regions) / self.inputs[0].length))

        return

    def has_long_ref_flanks(self, min_len=50):
        return len(self.ref_flanking_regions) == 2 and min(self.ref_flanking_regions) >= min_len

    def has_ref_deletion(self, min_diff=20):
        return len(self.start1_end1s) == 2 and abs(self.start1_end1s[1][0] - self.start1_end1s[0][1] - 1) >= min_diff

    def has_insertion(self, min_diff=1, max_diff=49, window=20):
        return len(self.start2_end2s) == 2 and min_diff <= abs(
            self.start2_end2s[1][0] - self.start2_end2s[0][1] - 1) <= max_diff and abs(
            self.start1_end1s[1][0] - self.start1_end1s[0][1] - 1) <= window

    def insertion_length(self):
        return self.excised_regions[1][0]

    def breakpoint_match(self, breakpoint, window=20):
        return min(map(lambda x: abs(x - breakpoint), list(itertools.chain.from_iterable(self.start1_end1s)))) <= window

    def has_long_flanks(self, min_len):
        return self.flanking_regions[0] >= min_len and self.flanking_regions[1] >= min_len + self.hom and min(
            self.ref_flanking_regions[0], self.ref_flanking_regions[1]) >= min_len

    def has_only_long_left_flank(self, min_len):
        return min(self.ref_flanking_regions[0], self.flanking_regions[0]) >= min_len and self.flanking_regions[
                                                                                              1] <= min_len

    def has_only_long_right_flank(self, min_len):
        return self.flanking_regions[0] <= min_len and self.flanking_regions[1] >= min_len + self.hom and \
               self.ref_flanking_regions[1] >= min_len

    def flanks_cover_first_seq(self, min_r):
        r = int(round(100.0 * sum(self.flanking_regions) / self.inputs[0].length))
        return r >= min_r

    def has_enough_identical(self, min_percent_identical):
        return self.percent >= min_percent_identical  # and min(self.percents) >= min_percent_identical

    def is_reference(self, min_flank_length=100, max_length_diff=20):
        if len(self.start1_end1s) == 1:
            if abs(self.start1_end1s[0][1] - self.start1_end1s[0][
                0] + 1) >= self.contig.sv_region.length() - 100: return True
            return False

        if min(self.flanking_regions) >= min_flank_length and min(self.ref_flanking_regions) >= min_flank_length:
            l1 = abs(self.start1_end1s[1][0] - self.start1_end1s[0][1] - 1)
            l2 = abs(self.start2_end2s[1][0] - self.start2_end2s[0][1] - 1)

            if abs(l1 - l2) <= max_length_diff:
                return True
        return False

    def almost_all_bases_aligned(self, min_unaligned=10):
        return self.aligned_bases + min_unaligned >= self.contig.sequence_len

    def __str__(self):
        return repr(self.__dict__)

    def __repr__(self):
        return str(self)
