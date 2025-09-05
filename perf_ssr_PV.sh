#!/bin/bash

module load CCEnv StdEnv/2020
source ~/.virtualenvs/perf_ssr/bin/activate

parallel -a final_continent_acc.txt PERF -i phase/genomes/{}.fna -o phase/perf_out/{}_perf.tsv -g phase/gff/{}.gff --gene-key locus_tag