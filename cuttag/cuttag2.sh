#!/bin/bash
###
 # @Author: error: error: git config user.name & please set dead value or install git && error: git config user.email & please set dead value or install git & please set dead value or install git
 # @Date: 2024-09-30 18:42:55
 # @LastEditors: error: error: git config user.name & please set dead value or install git && error: git config user.email & please set dead value or install git & please set dead value or install git
 # @LastEditTime: 2024-12-09 22:01:33
 # @FilePath: /root/code/upper/cuttag2.sh
 # @Description: 这是默认设置,请设置`customMade`, 打开koroFileHeader查看配置 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
### 

# 设置线程数
cores=16
# 初始化项目路径变量
projPath="/mnt/g/cuttag_3"
# 指定原始数据目录
rawdata_dir="${projPath}/raw2"
# 指定sam文件路径
sam_dir="${projPath}/alignment/sam"
# 指定bam目录
bam_dir="${projPath}/alignment/bam"

# 遍历Rawdata目录下的所有子文件夹
for sample_dir in "$rawdata_dir"/*/; do
    # 获取样本名称
    sample_name=$(basename "$sample_dir" "/")

    # 转换SAM到BAM格式
    samtools view -@ 15 -S ${sam_dir}/${sample_name}_bowtie2.sam -b > ${bam_dir}/${sample_name}.bam
    
    # 排序BAM文件
    samtools sort -@ 15 -l 5 -o ${bam_dir}/${sample_name}.bam.sort ${bam_dir}/${sample_name}.bam
    
    # 创建BAM索引
    samtools index -@ 25 ${bam_dir}/${sample_name}.bam.sort
    
    # 生成BAM统计信
    samtools flagstat ${bam_dir}/${sample_name}.bam.sort > ${bam_dir}/${sample_name}.stat
    
    # 生成BigWig文件
    bamCoverage -p 20 -v --binSize 50 --normalizeUsing RPKM --extendReads -b ${bam_dir}/${sample_name}.bam.sort -o ${bam_dir}/${sample_name}.bam.sort.bw
done
# 指定输出目录
output_dir="${projPath}/peaks"
# 指定阴性对照BAM文件路径
neg_bam_sort="${bam_dir}/CM-neg.bam.sort"
# 遍历Rawdata目录下的所有子文件夹
for sample_dir in "$rawdata_dir"/*/; do
    # 获取样本名称
    sample_name=$(basename "$sample_dir" "/")

    # 对样本进行callpeak处理，同时使用阴性对照
    macs3 callpeak -t "${bam_dir}/${sample_name}.bam.sort" -c "${neg_bam_sort}" -q 0.01 -g hs -B --nomodel --extsize 263 -n "${sample_name}" --outdir "${output_dir}"

    # 提取特定的统计信息
    egrep "tags after filtering in treatment" "${output_dir}/${sample_name}_peaks.xls" > "${output_dir}/${sample_name}_tags.txt"
done
