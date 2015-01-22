#!/bin/bash


METASV="run_metasv.py"
$METASV --pindel_native pindel/22_* --breakdancer_native breakdancer/22.out --cnvnator_native cnvnator/22.out --breakseq_vcf breakseq/breakseq.vcf.gz --reference reference/human_g1k_v37_decoy.fasta --outdir out --sample NA12878 --disable_assembly --filter_gaps --keep_standard_contigs

