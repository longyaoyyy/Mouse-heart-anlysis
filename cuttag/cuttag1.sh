#!/bin/bash
###
 # @Author: error: error: git config user.name & please set dead value or install git && error: git config user.email & please set dead value or install git & please set dead value or install git
 # @Date: 2024-09-25 22:37:28
 # @LastEditors: error: error: git config user.name & please set dead value or install git && error: git config user.email & please set dead value or install git & please set dead value or install git
 # @LastEditTime: 2024-12-09 19:56:33
 # @FilePath: /root/code/upper/cuttag1.sh
 # @Description: 这是默认设置,请设置`customMade`, 打开koroFileHeader查看配置 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
### 
###
 # @Author: error: error: git config user.name & please set dead value or install git && error: git config user.email & please set dead value or install git & please set dead value or install git
 # @Date: 2024-09-25 22:37:28
 # @LastEditors: error: error: git config user.name & please set dead value or install git && error: git config user.email & please set dead value or install git & please set dead value or install git
 # @LastEditTime: 2024-12-08 14:50:59
 # @FilePath: /root/code/upper/cuttag1.sh
 # @Description: 这是默认设置,请设置`customMade`, 打开koroFileHeader查看配置 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
### 

##== linux 命令 ==##
# 设置线程数
cores=16

# 初始化项目路径变量
projPath="/mnt/g/cuttag_3"

# 创建用于存放结果的目录
mkdir -p "${projPath}/fastqFileQC"
mkdir -p "${projPath}/alignment/sam/bowtie2_summary"
mkdir -p "${projPath}/alignment/bam"
mkdir -p "${projPath}/alignment/bed"
mkdir -p "${projPath}/alignment/bedgraph"

# 设置参考基因组索引路径
ref="/root/Database/bowtie2Index/hg38"

# 设置spike-in参考基因组路径
spikeInRef="/root/Database/bowtie2Index/E.colin"

# 设置染色体大小文件路径
chromSize="/root/Database/hg38/hg38.chrom.sizes"

# 指定原始数据目录
rawdata_dir="${projPath}/Rawdata"

# 接头去除工具路径
fastp_path=/usr/bin/fastp
# 指定输出目录
output_dir="${projPath}/fastqFileQC"

# 日志文件
log_file="${projPath}/alignment/sam/bowtie2_summary/script.log"

# 遍历Rawdata目录下的所有子文件夹
for sample_dir in "$rawdata_dir"/*/; do
    # 获取样本名称
    sample_name=$(basename "$sample_dir" "/")

    # 定义输入文件路径
    input_R1="$sample_dir/$sample_name"_R1.fq.gz
    input_R2="$sample_dir/$sample_name"_R2.fq.gz

    # 定义输出文件路径
    output_R1="${output_dir}/${sample_name}_trimmed_R1.fq.gz"
    output_R2="${output_dir}/${sample_name}_trimmed_R2.fq.gz"

    # 使用fastp进行接头去除
    $fastp_path -i "$input_R1" -I "$input_R2" -o "$output_R1" -O "$output_R2" --thread ${cores} --json "${output_dir}/${sample_name}_fastp.json" --html "${output_dir}/${sample_name}_fastp.html"

    # 使用bowtie2工具进行双端测序数据的比对
    bowtie2 --local --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700 -p ${cores} -x ${ref} -1 "$output_R1" -2 "$output_R2" -S "${projPath}/alignment/sam/${sample_name}_bowtie2.sam" &> "${projPath}/alignment/sam/bowtie2_summary/${sample_name}_bowtie2.txt"

    # 使用bowtie2工具对spike-in序列进行比对
    bowtie2 --local --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700 -p ${cores} -x ${spikeInRef} -1 "$output_R1" -2 "$output_R2" -S "${projPath}/alignment/sam/${sample_name}_bowtie2_spikeIn.sam" &> "${projPath}/alignment/sam/bowtie2_summary/${sample_name}_bowtie2_spikeIn.txt"

    # 计算序列深度
    seqDepthDouble=$(samtools view -F 0x04 "${projPath}/alignment/sam/${sample_name}_bowtie2_spikeIn.sam" | wc -l)
    seqDepth=$((seqDepthDouble / 2))

    # 输出序列深度到文件
    echo "$seqDepth" > "${projPath}/alignment/sam/bowtie2_summary/${sample_name}_bowtie2_spikeIn.seqDepth"

    # 记录日志
    echo "$(date): Finished processing sample $sample_name" >> "$log_file"
done