#!/bin/bash

#bgzip, index and list filtered vcf files
parallel bgzip {} ::: gatk_out/filtered_vcf/*.vcf
parallel bcftools index {} ::: gatk_out/filtered_vcf/*.vcf.gz

#merge files with bcftools only keeping PASS on filtering
bcftools merge -O z -o Calbicans_vcf.vcf.gz -m all --threads 80 -l Calbicans_vcfs.txt

#Index vcf for gatk
gatk IndexFeatureFile -I popnet_continent_vcf.vcf.gz

#make tsv table of variants
gatk VariantsToTable \
  -V popnet_continent_vcf.vcf.gz \
  -F CHROM -F POS -F REF -GF GT \
  -O continent_variant_table.tsv

#Run PopNet
python $HOME/bin/PopNet/PopNet.py /gpfs/fs0/scratch/j/jparkin/dcarru/ngo_WGS/IPNC/data/PopNet/popnet_config.txt