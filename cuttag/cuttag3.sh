#!/bin/bash

# 设置线程数
cores=16
# 初始化项目路径变量
projPath="/mnt/g/cuttag_3"
# 指定原始数据目录
rawdata_dir="${projPath}/Rawdata"
# 指定sam文件路径
sam_dir="${projPath}/alignment/sam"
# 指定bam目录
bam_dir="${projPath}/alignment/bam"
# 指定输出目录
output_dir="${projPath}/peaks"
# 指定阴性对照BAM文件路径
neg_bam_sort="${bam_dir}/AC-neg.bam.sort"
# 遍历Rawdata目录下的所有子文件夹
for sample_dir in "$rawdata_dir"/*/; do
    # 获取样本名称
    sample_name=$(basename "$sample_dir" "/")

    # 对样本进行callpeak处理，同时使用阴性对照
    macs3 callpeak -t "${bam_dir}/${sample_name}.bam.sort" -c "${neg_bam_sort}" -q 0.01 -g hs -B --nomodel --extsize 263 -n "${sample_name}" --outdir "${output_dir}"

    # 提取特定的统计信息
    egrep "tags after filtering in treatment" "${output_dir}/${sample_name}_peaks.xls" > "${output_dir}/${sample_name}_tags.txt"
done

