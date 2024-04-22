#!/bin/bash
# 定义文件名模式
file_pattern="*_1.fq.gz"

if [ ! -d "$BigWig" ]; then
    # 如果文件夹不存在，则创建它
    mkdir "BigWig"
    echo "文件夹 'BigWig' 已创建"
else
    echo "文件夹 'BigWig' 已经存在"
fi

if [ ! -d "$Bam" ]; then
    # 如果文件夹不存在，则创建它
    mkdir "Bam"
    echo "文件夹 'Bam' 已创建"
else
    echo "文件夹 'Bam' 已经存在"
fi

if [ ! -d "$TXT" ]; then
    # 如果文件夹不存在，则创建它
    mkdir "TXT"
    echo "文件夹 'TXT' 已创建"
else
    echo "文件夹 'TXT' 已经存在"
fi

# 使用 for 循环遍历匹配的文件
for file in $file_pattern
do
    echo "Step1: 解压文件: $file和${file/_1.fq.gz/_2.fq.gz}"
	
    gunzip "$file" & gunzip "${file/_1.fq.gz/_2.fq.gz}"

    echo "
	Step2: 使用 trim_galore 处理文件
	"
    trim_galore --cores 4 --paired "${file/_1.fq.gz/_1.fq}" "${file/_1.fq.gz/_2.fq}"

    echo "
	删除fastq文件
	"
     rm "${file/_1.fq.gz/_1.fq}" "${file/_1.fq.gz/_2.fq}"

    echo "
	Step3: 使用 bowtie2 对进行序列比对
	"
    bowtie2 -p 20 --reorder --very-sensitive-local -x mm10 --local --no-mixed --no-discordant -X 1000 --phred33 -1 "${file/_1.fq.gz/_1_val_1.fq}" -2 "${file/_1.fq.gz/_2_val_2.fq}" -S "${file/_1.fq.gz/}.sam" 
    
	echo "
	删除trim文件
	"
	rm "${file/_1.fq.gz/_1_val_1.fq}" "${file/_1.fq.gz/_2_val_2.fq}"

    echo "
	Step4: 将 SAM 文件转换为 BAM 文件
	"
    samtools view -@ 20 -bS "${file/_1.fq.gz/}.sam" > "${file/_1.fq.gz/}.bam"

    echo "
	删除中间生成的 SAM 文件
	"
    rm "${file/_1.fq.gz/}.sam"

    echo "
	Step5: bam文件排序
	"
    samtools sort  -@ 20 "${file/_1.fq.gz/}.bam" -o "${file/_1.fq.gz/_sorted}.bam" 

    echo "
	删除中间生成的 bam 文件
	"
    rm "${file/_1.fq.gz/}.bam"

    echo "
	Step6: picard标记重复reads
	"
    picard MarkDuplicates  --showHidden true --INPUT "${file/_1.fq.gz/_sorted}.bam"   --OUTPUT "${file/_1.fq.gz/_rmdup}.bam"    --METRICS_FILE "${file/_1.fq.gz/}".picard_info.txt --VALIDATION_STRINGENCY LENIENT
    rm "${file/_1.fq.gz/_sorted}.bam" 
      
	echo "
	Step7: 建立index文件
	"
    samtools index  -@ 20 "${file/_1.fq.gz/_rmdup}.bam" 
      
    echo "
	Step8: 利用bamcoverage生成可视化文件，参数： --normalizeUsing RPGC --binSize 20 --effectiveGenomeSize 2494787188（mm10） -p max --smoothLength 40  --ignoreDuplicates  --extendReads
	"
    bamCoverage --normalizeUsing RPGC --binSize 20 --effectiveGenomeSize 2494787188 -p max --smoothLength 40  --ignoreDuplicates  --extendReads  -b "${file/_1.fq.gz/_rmdup}.bam"  -o ./BigWig/"${file/_1.fq.gz/}.bw" 

	echo "
	Step9: 开始callpeak
	"
	
    macs2 callpeak -t "${file/_1.fq.gz/_rmdup}.bam"  -f BAMPE -g mm --outdir ./MACS2/"${file/_1.fq.gz/}"/ -n "${file/_1.fq.gz/-SPMR}" -q 0.01 --cutoff-analysis
    
    # echo "开始转格式BDGtoBigWig"
    # bedGraphToBigWig ./MACS2/"${file/_1.fq.gz/}"/"${file/_1.fq.gz/-SPMR}"_treat_pileup.bdg  ./BWA-MEM2_index/hg19.chrom.sizes ./BigWig/"${file/_1.fq.gz/_convert-SPMR}".bw
    find . -type f -name '*.bdg' -exec rm -f {} +
	mv *.bam *.bai ./Bam
    mv *.txt ./TXT
done
