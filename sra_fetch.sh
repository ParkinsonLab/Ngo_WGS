#!/bin/bash

#Load software/programs
module load CCEnv StdEnv/2020 gcc/9.3.0 sra-toolkit/3.0.0

#Prefetch then fasterq-dump to get SRA runs to fastq
parallel -a final_ref_acc.txt prefetch -O SRA_prefetch/
parallel -a final_ref_acc.txt fasterq-dump -O raw_fastq/ -f -p SRA_prefetch/{}