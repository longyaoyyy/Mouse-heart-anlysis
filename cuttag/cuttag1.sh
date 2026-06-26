#!/bin/bash

##== linux 命令 ==##
cores=16
projPath="/mnt/g/cuttag_3"
mkdir -p "${projPath}/fastqFileQC"
mkdir -p "${projPath}/alignment/sam/bowtie2_summary"
mkdir -p "${projPath}/alignment/bam"
mkdir -p "${projPath}/alignment/bed"
mkdir -p "${projPath}/alignment/bedgraph"

ref="/root/Database/bowtie2Index/hg38"
spikeInRef="/root/Database/bowtie2Index/E.colin"
chromSize="/root/Database/hg38/hg38.chrom.sizes"
rawdata_dir="${projPath}/Rawdata"
fastp_path=/usr/bin/fastp
output_dir="${projPath}/fastqFileQC"
log_file="${projPath}/alignment/sam/bowtie2_summary/script.log"

for sample_dir in "$rawdata_dir"/*/; do
    sample_name=$(basename "$sample_dir" "/")

    input_R1="$sample_dir/$sample_name"_R1.fq.gz
    input_R2="$sample_dir/$sample_name"_R2.fq.gz

    output_R1="${output_dir}/${sample_name}_trimmed_R1.fq.gz"
    output_R2="${output_dir}/${sample_name}_trimmed_R2.fq.gz"

    $fastp_path -i "$input_R1" -I "$input_R2" -o "$output_R1" -O "$output_R2" --thread ${cores} --json "${output_dir}/${sample_name}_fastp.json" --html "${output_dir}/${sample_name}_fastp.html"

    bowtie2 --local --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700 -p ${cores} -x ${ref} -1 "$output_R1" -2 "$output_R2" -S "${projPath}/alignment/sam/${sample_name}_bowtie2.sam" &> "${projPath}/alignment/sam/bowtie2_summary/${sample_name}_bowtie2.txt"

    bowtie2 --local --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700 -p ${cores} -x ${spikeInRef} -1 "$output_R1" -2 "$output_R2" -S "${projPath}/alignment/sam/${sample_name}_bowtie2_spikeIn.sam" &> "${projPath}/alignment/sam/bowtie2_summary/${sample_name}_bowtie2_spikeIn.txt"

    seqDepthDouble=$(samtools view -F 0x04 "${projPath}/alignment/sam/${sample_name}_bowtie2_spikeIn.sam" | wc -l)
    seqDepth=$((seqDepthDouble / 2))

    echo "$seqDepth" > "${projPath}/alignment/sam/bowtie2_summary/${sample_name}_bowtie2_spikeIn.seqDepth"

    echo "$(date): Finished processing sample $sample_name" >> "$log_file"
done
