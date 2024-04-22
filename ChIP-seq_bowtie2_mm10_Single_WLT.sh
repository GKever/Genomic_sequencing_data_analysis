#!/bin/bash
# 定义文件名模式
file_pattern="*.man"

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
    echo "解压文件: $file"
	
    fastq-dump  --split-files "$file" 

    echo "
	使用 trim_galore 处理文件
	"
    trim_galore --cores 4 "${file/.man/.man_1.fastq}" 

    echo "
	删除fastq文件
	"
    rm  "${file/.man/.man_1.fastq}" "$file" 

    echo "
	使用 bowtie2 对处理后的文件进行比对
	"
    bowtie2 -p 20 --reorder --sensitive-local -x mm10 --local --no-mixed --no-discordant -X 700 --phred33  "${file/.man/.man_1_trimmed.fq}" -S "${file/.man/}.sam" 
    
	echo "
	删除trim文件
	"
	rm "${file/.man/.man_1_trimmed.fq}" 

    echo "
	将 SAM 文件转换为 BAM 文件
	"
    samtools view -@ 20 -bS "${file/.man/}.sam" > "${file/.man/}.bam"

    echo "
	删除中间生成的 SAM 文件
	"
    rm "${file/.man/}.sam"

    echo "
	bam文件排序
	"
    samtools sort  -@ 20 "${file/.man/}.bam" -o "${file/.man/_sorted}.bam" 

    echo "
	删除中间生成的 bam 文件
	"
    rm "${file/.man/}.bam"

    echo "
	picard标记重复reads
	"
    picard MarkDuplicates  --showHidden true --INPUT "${file/.man/_sorted}.bam"   --OUTPUT "${file/.man/_rmdup}.bam"    --METRICS_FILE "${file/.man/}".picard_info.txt --VALIDATION_STRINGENCY LENIENT
    rm "${file/.man/_sorted}.bam" 
      
	echo "
	建立index文件
	"
    samtools index  -@ 20 "${file/.man/_rmdup}.bam" 
      
    echo "
	利用bamcoverage生成可视化文件，参数： --normalizeUsing RPGC --binSize 20 --effectiveGenomeSize 2864785220（hg19） -p max --smoothLength 40  --ignoreDuplicates  --extendReads
	"
    bamCoverage --normalizeUsing RPGC --binSize 20 --effectiveGenomeSize 2494787188 -p max --smoothLength 40  --ignoreDuplicates  -b "${file/.man/_rmdup}.bam"  -o ./BigWig/"${file/.man/}.bw" 

	echo "
	开始callpeak
	"
	
    macs2 callpeak -t "${file/.man/_rmdup}.bam"  -f BAMPE -g mm --SPMR --outdir ./MACS2/"${file/.man/}"/ -n "${file/.man/-SPMR}" -B -q 0.01 --cutoff-analysis
    
    # echo "开始转格式BDGtoBigWig"
    # bedGraphToBigWig ./MACS2/"${file/.man/}"/"${file/.man/-SPMR}"_treat_pileup.bdg  ./BWA-MEM2_index/hg19.chrom.sizes ./BigWig/"${file/.man/_convert-SPMR}".bw
    # find . -type f -name '*.bdg' -exec rm -f {} +
	mv *.bam *.bai ./Bam
    mv *.txt ./TXT
done
