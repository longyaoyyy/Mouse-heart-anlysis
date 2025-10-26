#!/bin/bash

# 设置基础路径
CLEANDATA_DIR="/mnt/g/scRNA/Cleandata"
TRANSCRIPTOME="/root/Database/10x/GRCh38"
# 使用本地目录作为输出路径，避免跨文件系统问题
OUTPUT_DIR="/root/cellranger_results"
FINAL_OUTPUT_DIR="/mnt/g/scRNA/cellranger_processed"
LOG_DIR="/root/cellranger_logs"

# 创建必要的目录
mkdir -p $OUTPUT_DIR
mkdir -p $FINAL_OUTPUT_DIR
mkdir -p $LOG_DIR

# 遍历每个样本目录
for sample_dir in $CLEANDATA_DIR/*/
do
    # 提取样本名称
    sample_name=$(basename $sample_dir)
    
    echo "Processing sample: $sample_name"
    
    # 检查最终输出目录是否已存在
    if [ -d "$FINAL_OUTPUT_DIR/$sample_name" ]; then
        echo "Output directory for $sample_name already exists. Skipping..."
        continue
    fi
    
    # 运行 cellranger count，使用本地目录作为输出路径
    cellranger count --id=$sample_name \
       --fastqs=$sample_dir \
       --sample=$sample_name \
       --transcriptome=$TRANSCRIPTOME \
       --create-bam=true \
       --localcores=20 \
       --localmem=120 \
       --output-dir=$OUTPUT_DIR/$sample_name \
       2>&1 | tee $LOG_DIR/${sample_name}_log.txt
    
    # 检查命令执行结果
    if [ ${PIPESTATUS[0]} -eq 0 ]; then
        echo "Successfully processed sample: $sample_name"
        # 处理成功后移动结果到最终目标位置
        mkdir -p $FINAL_OUTPUT_DIR/$sample_name
        mv $OUTPUT_DIR/$sample_name/outs $FINAL_OUTPUT_DIR/$sample_name
    else
        echo "Failed to process sample: $sample_name. Check log file: $LOG_DIR/${sample_name}_log.txt"
    fi
done

echo "Batch processing completed. Results in $FINAL_OUTPUT_DIR. Logs in $LOG_DIR."