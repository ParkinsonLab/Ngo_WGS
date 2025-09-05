#!/bin/bash

#Start by cleaning input file
echo -e "1\naligned_gene_sequences/$1_gene_0.fas\n5\n$1_clean.fasta" | hyphy rmv

#Make ML tree using iqtree

iqtree2 -s $1_clean.fasta -B 1000 --prefix $1 --quiet

#Start with BUSTED for per gene selection

hyphy busted --alignment $1_clean.fasta --tree $1.contree --multiple-hits Double

#Feed consensus ML tree into hyphy fubar and redirect to txt file

hyphy fubar --alignment $1_clean.fasta --tree $1.contree 

