#!/bin/bash

cd /camp/lab/turnerj/working/shared_projects/LP_WGS_karyo/round2

LIB=$(sed -n "${SLURM_ARRAY_TASK_ID}p" trimmed)
VAR=$(echo $LIB | rev | cut -d "/" -f 1 | rev)
SAMPLE=$(echo $VAR | cut -d "_" -f 1)

READ1="$LIB"_R1_001_trimmed.fq.gz
READ2="$LIB"_R2_001_val_2.fq.gz
READ3="$LIB"_unpaired.fq.gz

OUT=/camp/lab/turnerj/working/shared_projects/LP_WGS_karyo/round2/bams
cd $OUT
echo "cd $OUT"

GENOME=/camp/svc/reference/Genomics/iGenomes/Mus_musculus/UCSC/mm10/Sequence/BWAIndex/genome.fa

echo "bwa mem -t 4 -B 2 -O 4,4 -L 3,3 -U 5 $GENOME $READ1 -R "@RG\tID:${SAMPLE}\tLB:${SAMPLE}\tSM:${SAMPLE}\tPL:"ILLUMINA"" > "$SAMPLE".sam"
bwa mem -t 4 -B 2 -O 4,4 -L 3,3 -U 5 $GENOME $READ1 -R "@RG\tID:${SAMPLE}\tLB:${SAMPLE}\tSM:${SAMPLE}\tPL:"ILLUMINA"" > "$SAMPLE".sam
bwa mem -t 4 -B 2 -O 4,4 -L 3,3 -U 5 $GENOME $READ3 -R "@RG\tID:${SAMPLE}\tLB:${SAMPLE}\tSM:${SAMPLE}\tPL:"ILLUMINA"" > "$SAMPLE"_unpaired.sam


echo "samtools view -b "$SAMPLE".sam -o "$SAMPLE".bam"
samtools view -b "$SAMPLE".sam -o "$SAMPLE".bam
echo "samtools sort "$SAMPLE".bam -o "$SAMPLE".sorted.bam"
samtools sort "$SAMPLE".bam -o "$SAMPLE".sorted.bam
echo "samtools index "$SAMPLE".sorted.bam"
samtools index "$SAMPLE".sorted.bam
rm "$SAMPLE".sam "$SAMPLE".bam

samtools view -b "$SAMPLE"_unpaired.sam -o "$SAMPLE"_unpaired.bam
samtools sort "$SAMPLE"_unpaired.bam -o "$SAMPLE"_unpaired.sorted.bam
samtools index "$SAMPLE"_unpaired.sorted.bam
rm "$SAMPLE"_unpaired.sam "$SAMPLE"_unpaired.bam

samtools merge -r "$SAMPLE".merged.bam "$SAMPLE".sorted.bam "$SAMPLE"_unpaired.sorted.bam
samtools sort "$SAMPLE".merged.bam -o "$SAMPLE".merged.sorted.bam
samtools index "$SAMPLE".merged.sorted.bam
rm "$SAMPLE".merged.bam



BIN=750
DIR=/camp/lab/turnerj/working/shared_projects/LP_WGS_karyo/round2
FILES=$DIR/bamfiles
ALLCHR="T"
OUT=$DIR/plots/

echo "running Rscript 'qdnaseq.r' with bin size $BIN for file number $SLURM_ARRAY_TASK_ID from list $FILES in $DIR"
echo "plotting all chromosomes is set to: $ALLCHR"

Rscript qdnaseq.r $SLURM_ARRAY_TASK_ID $BIN $FILES $DIR $ALLCHR $OUT
