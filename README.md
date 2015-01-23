MetaSV
===========

MetaSV: An accurate and integrative structural-variant caller for next generation sequencing

# Software Dependencies
MetaSV uses the following Python packages which must be available on the `PYTHONPATH` environment variable:
* [pysam](http://pysam.readthedocs.org/en/latest/): Only version 0.7.7 has been tested to work with MetaSV
* [pybedtools](http://pythonhosted.org/pybedtools/)
* [pyvcf](https://github.com/jamescasbon/PyVCF)

In addition, it requires the following software to be installed:
* [SPAdes](http://bioinf.spbau.ru/spades): Use for assembly around the SV breakpoints
* [AGE](https://github.com/abyzovlab/AGE): Used for determining the SV breakpoints using assembled contigs
The paths to the above tools are specified on the command-line options for MetaSV.

# Tools
* `metasv.py`: Top-level MetaSV script. Use this when running all of MetaSV.
* `breakdancer_reader.py`: Reader for BreakDancer output files. Can also convert to a VCF file.

# Installing
`pip install https://github.com/bioinform/metasv/archive/<version>-alpha.tar.gz`. The current version of MetaSV is 0.1-alpha.

# Running
Type `metasv.py -h` for help on the command-line options.

