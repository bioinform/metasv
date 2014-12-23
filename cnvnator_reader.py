import logging
import re

logger = logging.getLogger(__name__)

'''
From the CNVnator README
The output is as follows:

CNV_type coordinates CNV_size normalized_RD e-val1 e-val2 e-val3 e-val4 q0

normalized_RD -- normalized to 1.
e-val1        -- is calculated using t-test statistics.
e-val2        -- is from probability of RD values within the region to be in
the tails of gaussian distribution describing frequencies of RD values in bins.
e-val3        -- same as e-val1 but for the middle of CNV
e-val4        -- same as e-val2 but for the middle of CNV
q0            -- fraction of reads mapped with q0 quality
'''

pattern = re.compile("(:|-)")

sv_type_dict = {"deletion": "DEL", "duplication": "DUP"}
class CNVnatorRecord:
  def __init__(self, record_string):
    fields = record_string.split()
    self.sv_type = sv_type_dict[fields[0]]

    coordinates = pattern.split(fields[1])
    self.chromosome = coordinates[0]
    self.start = int(coordinates[2])
    self.end = int(coordinates[4])
    self.sv_len = self.end - self.start
    self.normalized_rd = float(fields[3])
    self.e_val1 = float(fields[4])
    self.e_val2 = float(fields[5])
    self.e_val3 = float(fields[6])
    self.e_val4 = float(fields[7])
    self.q0 = float(fields[8])

  def __str__(self):
    return str(self.__dict__)

class CNVnatorReader:
  def __init__(self, file_name):
    logger.info("File is " + file_name)
    self.file_fd = open(file_name)

  def __iter__(self):
    return self

  def next(self):
    while True:
      line = self.file_fd.next()
      if line[0] != "#":
        return CNVnatorRecord(line.strip())
