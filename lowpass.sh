#!/bin/bash

# example execution: ./lowpass.sh A15 500

# variable definitions

sample=$1
binSize=$2
genome=Mus_musculus/UCSC/mm10/Sequence/BWAIndex/genome.fa

# directory definitions

outDir=/user/home/working_dir/results
fastqDir=/user/home/working_dir/fastq

# constructing read file names from sample name

r1=$fastqDir/${sample}_R1.fastq.gz
r2=$fastqDir/${sample}_R2.fastq.gz
r3=$fastqDir/${sample}_unpaired.fastq.gz

# read mappings with BWA

echo "running read mapping with bwa and postprocessing with samtools"

bwa mem -t 4 -B 2 -O 4,4 -L 3,3 -U 5 $genome $r1 $r2 -R "@RG\tID:${sample}\tLB:${sample}\tSM:${sample}\tPL:"ILLUMINA"" | samtools view -b - | samtools sort - -o $outDir/${sample}.bam
samtools index $outDir/${sample}.bam

bwa mem -t 4 -B 2 -O 4,4 -L 3,3 -U 5 $genome $r3 -R "@RG\tID:${sample}\tLB:${sample}\tSM:${sample}\tPL:"ILLUMINA"" | samtools view -b - | samtools sort - -o $outDir/${sample}_unpaired.bam
samtools index $outDir/${sample}_unpaired.bam

# merging the BAM files from paired- and single-ended mappings

samtools merge -r - $outDir/${sample}.bam $outDir/${sample}_unpaired.bam | samtools sort - -o $outDir/${sample}_merged.bam
samtools index $outDir/${sample}_merged.bam

# starting the analysis in R

echo "running Rscript 'qdnaseq.r' with bin size $binSize for sample $sample"

Rscript qdnaseq.r $sample $binSize $outDir
