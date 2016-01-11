import logging
import itertools
import re

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


class AgeFormatError(Exception):
    def __init__(self, context, line):
        self.context = context
        self.line = line
    def __str__(self):
        return "Error while reading AGE output, L%d (section %s)." % (self.line, self.context)


class AgeRecord:
    def __init__(self, age_out_file=None,tr_region_1=[]):
        # Index 0 corresponds to the reference sequence
        self.aligned_bases = 0
        self.inputs = []  # inputs[1] is expected to have an insertion w.r.t. inputs[0]
        self.percent = 100  # Percent of matching bases in alignment
        self.percents = [100, 100]  # Percent of matching bases in alignment for each flank
        self.start1_end1s = []  # Alignment intervals for the first sequence
        self.start1_end1s_orig = []  # Alignment intervals for the first sequence (without considering truncations)
        self.start2_end2s = []  # Alignment intervals for the second sequence
        self.polarities1 = []  # Alignment polarities of intervals for the second sequence
        self.polarities2 = []  # Alignment polarities of intervals for the second sequence
        self.excised_regions = []  # List of excised regions as intervals
        self.n_alt = 0  # Number of alternate regions
        self.alternate_regions = []

        self.flanking_regions = [0, 0]  # Length of each flank on the donor sequence
        self.ref_flanking_regions = [0, 0]  # Length of each flank on the ref requence
        self.hom = 0  # Length of sequence homology
        self.flank_percent = 0  # Percent of total length in the flanks
        self.nfrags = 1  # Number of fragments in the alignment

        # The following fields are set later on when processing the alignments
        self.very_bad_assembly = False  # Alignment failed some basic checks
        self.bad_assembly = False  # Only one fragment in the alignment

        self.start = -1
        self.end = -1
        self.excised_lengths = [0, 0]  # Lengths of each excised region

        self.assembly_contig = None
        self.used = False
        self.duplicate = False

        self.age_file = age_out_file
        self.score = 0
        self.tr_region_1=tr_region_1 # Truncation region in first sequence
        self.invalid_tr_intervals=[] #Indicies of invalid intervals in first sequence due to truncation
        
        if age_out_file is not None:
            self.read_from_age_file(age_out_file)

    class LineReader:
        """Behaves like file objects do but also counts lines for tracing the AGE parser."""
        def __init__(self, fd):
            self.fd = fd
            self.line_num = 0
        def readline(self):
            line = self.fd.readline()
            if line: self.line_num += 1
            return line

    rx_rng = re.compile(r"\[\s*(\d+),\s*(\d+)\]")
    rx_perc = re.compile(r"\(\s*(\d+)%\)")
    rx_input = re.compile(r"=>\s+(\d+) nucs '(.*?)'")
    rx_ex = re.compile(r"seq =>\s+(\d+) nucs( \[(\d+),(\d+)\])?")
    rx_ident = re.compile(r"seq =>\s+(\d+) nucs( \[(\d+),(\d+)\] to \[(\d+),(\d+)\])?")

    def read_alignment_ranges(self, age_fd, name):
        line = age_fd.readline().strip()
        if not line.startswith(name):
            raise AgeFormatError("ALIGNMENT RANGES", age_fd.line_num)
        start_end = []
        for m in self.rx_rng.finditer(line):
            start_end.append([int(m.group(1)), int(m.group(2))])
        return start_end

    def parse_input_descriptor(self, line):
        """Returns a tuple of filename and sequence length."""
        m = self.rx_input.search(line)
        if m is None:
            raise AgeFormatError("INPUT DESCRIPTOR", age_fd.line_num)
        return (m.group(2), int(m.group(1)))

    def read_excluded_range(self, age_fd, name, context):
        line = age_fd.readline().strip()
        if len(line) == 0:
            return None
        if not line.startswith(name):
            raise AgeFormatError(context, age_fd.line_num)
        m = self.rx_ex.search(line)
        seq_len = int(m.group(1))
        # ATTN: Do you intend to ever use start or stop of excluded regions? So far those are never read anywhere. 
        return [0] if seq_len == 0 else [seq_len, int(m.group(3)), int(m.group(4))]

    def read_excluded_regions(self, age_fd, context):
        excluded_regions = [];
        while True:
            seq_rng = self.read_excluded_range(age_fd, "first", context)
            if seq_rng is None:
                break
            excluded_regions.append(seq_rng)
            seq_rng = self.read_excluded_range(age_fd, "second", context)
            if seq_rng is None:
                raise AgeFormatError(context, age_fd.line_num)
            excluded_regions.append(seq_rng)
            # ATTN: Is it intentional that while AGE allows multiple excised regions, only the first is ever looked at?
        return excluded_regions

    def read_identities(self, age_fd, name, context):
        line = age_fd.readline().strip()
        if not line.startswith(name):
            raise AgeFormatError(context, age_fd.line_num)
        m = self.rx_ident.search(line)
        seq_len = int(m.group(1))
        # ATTN: Do you intend to ever use positions of identities? So far those are never read anywhere. 
        return [0] if seq_len == 0 else [seq_len, int(m.group(3)), int(m.group(4)), int(m.group(5)), int(m.group(6))]

    def read_from_age_file(self, age_out_file):
        with open(age_out_file) as raw_age_fd:
            age_fd = self.LineReader(raw_age_fd)
            while True:
                line = age_fd.readline()
                if not line: break
                line = line.strip()

                if line.startswith("Aligned:"):
                    self.aligned_bases = int(line.split()[1])
                    continue

                if line.startswith("First  seq"):
                    file1, len1 = self.parse_input_descriptor(line)
                    if self.tr_region_1:
                        len1+=self.tr_region_1[1]                        
                    self.inputs.append(AgeInput(file1, len1))
                    continue

                if line.startswith("Second seq"):
                    file2, len2 = self.parse_input_descriptor(line)
                    self.inputs.append(AgeInput(file2, len2))
                    continue

                if line.startswith("Identic:"):
                    nums = map(int, self.rx_perc.findall(line))
                    self.percent = nums[0]
                    if len(nums) > 1:
                        self.percents = nums[1:3]
                    continue

                if line.startswith("Score:"):
                    self.score = int(line.split()[1])
                    continue

                if line == "Alignment:":
                    self.start1_end1s_orig = self.read_alignment_ranges(age_fd, "first")
                    self.start2_end2s = self.read_alignment_ranges(age_fd, "second")
                    if self.tr_region_1:
                        self.update_pos_tr()
                        self.invalid_tr_intervals=[i for i,interval in enumerate(self.start1_end1s_orig) if 
                                                     (interval[0]< self.tr_region_1[0]) ^ (interval[1]< self.tr_region_1[0])]
                        if self.invalid_tr_intervals:
                            logger.warn("These intervals has problems (spanned over truncation): %s" % self.invalid_tr_intervals)
                    else:
                        self.start1_end1s=self.start1_end1s_orig
                    self.polarities1 = map(lambda y: 1 if y[1]>y[0] else -1,self.start1_end1s)
                    self.polarities2 = map(lambda y: 1 if y[1]>y[0] else -1,self.start2_end2s)
                    self.nfrags = len(self.start1_end1s)
                    self.flanking_regions[0] = abs(self.start2_end2s[0][1] - self.start2_end2s[0][0] + 1)
                    self.flanking_regions[1] = 0 if len(self.start2_end2s) == 1 else abs(
                        self.start2_end2s[1][1] - self.start2_end2s[1][0] + 1)
                    self.ref_flanking_regions[0] = abs(self.start1_end1s[0][1] - self.start1_end1s[0][0] + 1)
                    self.ref_flanking_regions[1] = 0 if len(self.start1_end1s) == 1 else abs(
                        self.start1_end1s[1][1] - self.start1_end1s[1][0] + 1)
                    continue

                #TODO: May need to fix EXCISED REGIONS for truncated regions
                if line == "EXCISED REGION(S):":
                    self.excised_regions = self.read_excluded_regions(age_fd, line)
                    continue

                #TODO: May need to fix ALTERNATIVE REGION for truncated regions
                if line == "ALTERNATIVE REGION(S):":
                    self.alternate_regions = self.read_excluded_regions(age_fd, line)
                    continue

                #TODO: May need to fix breakpoint_identities for truncated regions
                if line == "Identity at breakpoints:":
                    self.read_identities(age_fd, "first", "IDENTITIES") # skip first sequence
                    breakpoint_identities = self.read_identities(age_fd, "second", "IDENTITIES")
                    # ATTN: Is it intentional to only examine the first occurrence?
                    self.hom = max(0, breakpoint_identities[0])
                    continue

        if len(self.inputs) == 0:
            self.flank_percent = 0
            logger.warn("%s has problems" % age_out_file)
        else:
            self.flank_percent = int(round(100.0 * sum(self.flanking_regions) / self.inputs[0].length))

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
        return min(self.ref_flanking_regions[0], self.flanking_regions[0]) >= min_len >= self.flanking_regions[1]

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

    def update_pos_tr(self):
        self.start1_end1s = []
        tr_p,tr_l = self.tr_region_1
        for i,ends in enumerate(self.start1_end1s_orig):
            dist_to_tr = map(lambda x: x-tr_p, ends)
            new_ends = [0,0]
            if (dist_to_tr[0] < 0 and  dist_to_tr[0] >=0) or (dist_to_tr[0] >= 0 and  dist_to_tr[1] <0) :
                j = filter(lambda x: ends[x]<tr_p,[0,1])[0]
                if abs(ends[j]-tr_p)>abs(ends[1-j]-tr_p):                    
                    new_ends[j] = ends[j]
                    new_ends[1-j] = tr_p-1
                    self.start2_end2s[i][1-j] -= (ends[1-j] - tr_p+1)
                else:
                    new_ends[1-j] = ends[1-j] + tr_l
                    new_ends[j] = tr_p + tr_l 
                    self.start2_end2s[i][j] += (tr_p - ends[j])
            elif dist_to_tr[0] >= 0:
                new_ends = map(lambda x: x+tr_l, ends)
            else:
                new_ends = ends                
            self.start1_end1s.append(new_ends)

    def __str__(self):
        return repr(self.__dict__)

    def __repr__(self):
        return str(self)
