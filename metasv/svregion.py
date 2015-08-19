class SVRegion:
    def __init__(self, chrom1, pos1, chrom2, pos2):
        self.chrom1 = chrom1
        self.pos1 = pos1
        self.chrom2 = chrom2
        self.pos2 = pos2

    def __str__(self):
        return "(%s, %d, %s, %d)" % (self.chrom1, self.pos1, self.chrom2, self.pos2)

    def __repr__(self):
        return str(self)

    def to_tuple(self):
        return self.chrom1, self.pos1, self.chrom2, self.pos2

    def length(self):
        return self.pos2 - self.pos1
