#!/bin/bash

#import software/programs
module load CCEnv StdEnv/2020 gcc/9.3.0 python/3.11.5 blast+/2.13.0 spades/3.15.4 prokka/1.14.5 prank mafft clustal-omega mash cd-hit

#Run unicycler for assemblies
parallel -j 8 -S $HOSTS --env PATH,LD_LIBRARY_PATH --wd $PWD -a continent_acc.txt ~/bin/Unicycler/unicycler-runner.py -1 trim/{}_1.fastq -2 trim/{}_2.fastq -o assemblies/{} -t 10

#Run prokka for gene annotation with ngo training and ngo ref proteins
prokka \
  --outdir prokka/$1 \
  --prefix $1 \
  --addgenes \
  --locustag $1 \
  --genus Neisseria \
  --kingdom Bacteria \
  --proteins $HOME/prokka_databases/ngo_prokka_ref.gb \
  --force \
  --cpus 20 \
  --centre X \
  --compliant \
  --prodigaltf $HOME/prokka_databases/ngo_training.trn \
  assembled_genomes/$1.fasta

#Run panaroo to construct pangenome
panaroo -i continent_gff/*.gff -o panaroo_results/ --clean-mode strict -a core --aligner mafft -t 80
