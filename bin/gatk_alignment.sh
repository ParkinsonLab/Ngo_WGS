#!/bin/bash

sample=$1

#fastp trimming
  fastp \
		--in1 ./raw_fastq/${sample}_R1.fastq \
    --in2 ./raw_fastq/${sample}_R2.fastq \
    --out1 ./trim/${sample}_1.fastq \
    --out2 ./trim/${sample}_2.fastq \
    --qualified_quality_phred 20 \
    --unqualified_percent_limit 40 \
    --cut_mean_quality 20 \
    --cut_window_size 5 \
    --cut_front \
    --cut_tail \
    --correction \
    --thread 4 \
    --html ./report/html/${sample}.html \
    --json ./report/json/${sample}.json 

#fastq to SAM
java -jar $EBROOTPICARD/picard.jar FastqToSam \
			    FASTQ=trim/${sample}_1.fastq \
			    FASTQ2=trim/${sample}_2.fastq \
			    OUTPUT=gatk_out/unmapped/${sample}.unmap.sam \
			    READ_GROUP_NAME=$sample \
			    SAMPLE_NAME=$sample \
			    LIBRARY_NAME=$sample \
			    PLATFORM="illumina" \
			    VERBOSITY='ERROR' \
			    SORT_ORDER="queryname"

#map to reference
bwa mem -t 4 reference/ngo_reference.fasta \
			trim/${sample}_1.fastq \
			trim/${sample}_2.fastq \
			-o gatk_out/aligned/${sample}.align.sam

#merge with unmapped
java -jar $EBROOTPICARD/picard.jar MergeBamAlignment \
			ALIGNED=gatk_out/aligned/${sample}.align.sam \
			UNMAPPED=gatk_out/unmapped/${sample}.unmap.sam \
			O=gatk_out/merged/${sample}.merged.bam \
			R=reference/ngo_reference.fasta \
			VERBOSITY='ERROR' \
			CLIP_ADAPTERS=false \
			CLIP_OVERLAPPING_READS=false

#De-duplicate  
gatk MarkDuplicatesSpark \
			-I gatk_out/merged/${sample}.merged.bam \
			-O gatk_out/markedDups/${sample}.mark.bam \
			--verbosity 'ERROR'  \
			--conf "spark.executor.cores=2"

#bcftools for calling variants and generating gvcf
bcftools mpileup -Ou -d 5200 -C 50 -q 20 -a AD,ADF,ADR,DP,SP -f reference/ngo_reference.fasta gatk_out/markedDups/${sample}.mark.bam | \
bcftools call --ploidy 1 -mv -Ob | \
bcftools filter -g 3 -G 10 -e '%QUAL < 50 && FMT/AD < 8 && FMT/ADF < 3 && FMT/ADR < 3 && MQ < 20' -o gatk_out/filtered_vcf/${sample}.filter.vcf 
bgzip gatk_out/filtered_vcf/${sample}.filter.vcf 
bcftools index gatk_out/filtered_vcf/${sample}.filter.vcf.gz
      
#remove some temporary files
rm -f gatk_out/unmapped/${sample}.unmap.sam
rm -f gatk_out/aligned/${sample}.align.sam
rm -f gatk_out/merged/${sample}.merged.bam
