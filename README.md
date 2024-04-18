这里主要记录一些我用于基因组测序数据分析的一些笔记和代码，有空的话我会不断地进行一些更新。

![20240325-1](https://github.com/GKever/Genomic_sequencing_data_analysis/assets/111635048/6bbf124c-aca9-4e18-b781-235fa99e85d9)

**一、搭建分析平台**

测序分析过程基本都时基于Linux平台，可以在服务器上部署也可以在windows/Mac电脑上部署，强烈建议在Windows中使用WSL功能，安装Linux子系统，本文后续都是基于WSL平台（ubuntu 18.04)进行。具体方法为：

1. windows启用wsl功能

   win11中在设置中，系统->可选功能->更多windows功能->勾选“适用于linux的windows子系统”,“虚拟机平台”，“Windows虚拟机监控程序平台”，确定后重启电脑。

2. 安装Linux分发版本

   在Windows商店中搜索Linux，安装自己喜欢的Linux版本，我个人常用Ubuntu系统。

3. Linux系统中搭建分析流程

   启动Linux系统，安装conda软件（https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html），

```
mkdir -p ~/miniconda3
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
rm -rf ~/miniconda3/miniconda.sh
```
   安装完成后，重启下Linux系统以启用conda。
   
4. 搭建环境

   在Linux中新建一个文件夹“seq”，下载名为"Genomic_sequence_data_analysis_workflow_GK.yaml"的conda虚拟环境配置文件，该虚拟环境基于python3.7，包含众多常见分析软件：  
  - bedtools=2.30.0=hc088bd4_0
  - bowtie2=2.3.5.1=py37he513fc3_0
  - bwa-mem2=2.2.1=h9a82719_1
  - bzip2=1.0.8=hd590300_5
  - circos=0.69.8=hdfd78af_1
  - cutadapt=1.18=py37h14c3975_1
  - deeptools=3.5.1=py_0
  - deeptoolsintervals=0.1.9=py37h8902056_4
  - fastp=0.12.4=0
  - fastqc=0.11.9=hdfd78af_1
  - intervene=0.6.4=pyh864c0ab_3
  - macs2=2.2.6=py37h516909a_0
  - picard=2.25.1=hdfd78af_1
  - samtools=1.7=1
  - seqkit=2.5.1=h9ee0642_0
  - trim-galore=0.6.10=hdfd78af_0

```
conda env create Genomic_sequence_data_analysis_workflow_GK.yaml
```
   等待安装完成，随后
```
conda activate seq
```
**二、数据分析**

写在开头，对于测序过程中的一些基本概念，比如sequence，fragments，reads，peak，bin等，需要有一些基本理解，

fragments：这个一般指片段化之后的DNA片段，和sequence的含义基本一致；

reads：读长，指测序后的碱基序列，不等同于fragments;

alignment: 比对，把reads比对至基因组对应的位置上，等同于mapping;

genomic build：指基因组版本号，序列比对时需要参考基因组，不同物种有各自的参考序列，同一物种的基因组序列也有不同的组装版本，如GRCh38（hg19）/GRCh39，GRCm38（mm10）/GRCm39。

concordantly mapping：对于双端测序，一段fragment产生的两个reads应为相向分布，此时称为concordantly。但由于比对错误或者某些特殊序列，导致paired reads分布不符合相向特征，如同向，不在同一位置等，称为disconcordantly，这样的序列会在后续分析中剔除。

pile up：堆积，一些reads同时集中于某一区域，产生堆积效应，容易形成峰状图形；

另外，一些常见文件格式也要熟悉，如fasta（fa），fastq（fq），SAM，BAM，bw（bigwig/wiggle），bed，bdg（bedgraph），gff，narrowpeak，broadpeak等


**1. 质控Quality control**

   一般使用fastQC工具来了解测序结果的质量，一般看看就行了，没啥太大用，结果不好也改变不了什么，也可省略这一步

```fastqc xxx.fq```

   得到xxx.html的网页文件，打开就可看到对reads评价的各项结果。
   
**2. 数据过滤，去除barcode/adaptor**

   原始的fastq文件往往是包含adaptor以及barcode的。如果不去除adaptor将会对下游的mapping产生影响。如果序列很短而3'端又存在adaptor的话，影响会较为严重。这里主要使用trim_galore工具进行数据过滤和接头序列修剪，注意，对于包含接头的序列
   
   对于双端测序（paired-end）使用以下指令：
   
```trim_galore --paired AAA_1.fq AAA_2.fq```
   
   Trim之后会得到以“AAA_1.val_1.fq”和“AAA_2.val_2.fq”和结尾的文件。
   
   对于单端测序（single-end）使用以下指令：
   
```trim_galore AAA.fq ```
   
   Trim之后会得到以“AAA.trimmed.fq”结尾的文件。

**3. 序列比对Mapping**

   对于基因组测序结果分析，常用的mapping软件有BWA-MEM2以及bowtie2。这里以bowtie2为例，演示比对过程。需要提前准备index文件，常见物种的index文件可以在bowtie2网站中下载（https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml）

 ```
 ##对于单端测序使用一下指令：
   
 bowtie2 -p 20 --reorder --very-sensitive-local -x hg19 -U AAA. val_1.fq -S AAA.sam

##对于双端测序：
   
bowtie2 -p 20 --reorder --very-sensitive-local -x hg19 -X 700 -1 AAA_1_val_1.fq -2 AAA_2_val_2.fq -S AAA.sam

##对于CUT&Tag测序结果可以使用以下优化分析策略：

bowtie2 -p 20 –reorder --very-sensitive-local -x hg19 –local --no-mixed --no-discordant --phred33 -I 10 -X 700 -1 AAA_1_val_1.fq -2 AAA_2_val_2.fq -S XXXX.sam

##对于ATAC-seq测序结果可考虑以下分析策略：

bowtie2 -p 20 --reorder --very-sensitive-local -X 2000 -x hg19 -1 AAA_1_val_1.fq -2 AAA_2_val_2.fq -S AAA.sam
```
   
   Bowtie2和BWA-MEM2都是比较优秀的比对软件，BWA-MEM2的速度更快，比对率也相对更优秀些，而Bowtie2优势在于可以设置不同参数对结果进行优化，比如-N可以设置允许错配碱基数（最多允许错配一个碱基），对于某些特定存在碱基突变细胞，该参数尤其有效。
   比对策略中，--very-sensitive-local是灵敏度最高的模式，但比对时间会大大增加，可优先考虑--sensitive-local模式，在该模式下通常一个3000w reads的测序文件能够在10分钟左右比对完成（i9 10900K），而前者大概需要2-3h。
   另外有一点，对于双端测序文件，bowtie2默认允许的最大比对sequence长度为500，可以通过设置-X参数将这一条件放宽值1000（根据文库的片段分布情况），特别是对于ATAC-seq文库，片段大小分布比较宽，短片段（低于200bp）和长片段（大于500bp）都很多，为了保证数据不丢失，在trim-galore（清除接头）过程中要保留短reads（低于150bp），而比对过程中要保留长reads。

**4. SAM/BAM文件处理**

   比对结束后，会得到一个扩展名为sam的文件。SAM的全称是sequence alignment/map format。而BAM就是SAM的二进制文件(B取自binary)。
```
##1) SAM，BAM格式转换，缩减文件体积，方便后续处理

samtools view -@ 20 -bS AAA.sam > AAA.bam

##2) sorting，将测序结果按照其在基因组上的坐标位置进行排序

samtools sort -@ 20 AAA.bam -o AAA.sort.bam

##3）mark duplicates，标记重复reads，由于PCR会产生PCR duplicates，测序仪器可能产生optical duplicates，带来冗余信息，需要对其标记

picard MarkDuplicates  --showHidden true --INPUT AAA.sort.bam --OUTPUT AAA.rmdup.bam --METRICS_FILE AAA.picard_info.txt --VALIDATION_STRINGENCY LENIENT

##4) 创建index文件，一些基于bam文件的处理需要先建立索引文件，生成AAA.rmdup.bai文件

samtools index -@ 20 AAA.rmdup.bam 
```

**5. Visualizing data，结果可视化**

   在基因组浏览器UCSC genome browser或IGV中可以查看mapping后的结果，可以用于在genome browser或者viewer中显示的文件格式有：bed, wiggle, bedGraph, bam, gff/gtf等等。最常见的就是bigwig文件，该文件可以通过bam文件或者bgd文件生成。bam文件是比对后最接近原始结果，建议用bam文件生成bigwig。
```
bamCoverage --normalizeUsing RPGC --effectiveGenomeSize 2864785220 --binSize 10 -p max –smoothLength 40 --ignoreDuplicates --centerReads -b AAA.rmdup.bam -o ./BigWig/AAA.bw
```
   对于基因组测序结果建议使用RPGC进行归一化，需要设置--effectiveGenomeSize参数，不同基因组版本的参数值不一样，具体参考网站（https://deeptools.readthedocs.io/en/latest/content/feature/effectiveGenomeSize.html）。建议选择弃掉duplicates，添加blacklist。
   
   对于双端测序结果，**一定要使用--extendReads参数进行优化**（不用设置参数值），而单端测序结果需要再根据平均sequence长度设置参数值。**双端测序不可同时设置--centerReads和--extendReads**。

   在IGV中观察可视化结果时，对于ChIP-seq，应观察1）峰是否完整，清楚。还是连成一片；2）它是否出现的正确的位置；3）如果是Transcript factor，峰是否是在它应该在的位置上，比如TSS的上游？如果是H3K4me3，峰是否富集的TSS附近？如果是H3K4me1/2,H3/H4ac, DNase，峰是否在TSS附近或者调控单元的末端？如果是H3K36me3，峰是否全基因都富集？如果是H3K27me3，是否富集在inactive基因的CpG island附近？如果H3K9me3，峰是否出现在broad domains和repeat elements上；4）背景是否足够低；5）检查有无明显的尖刺样的信号（spikes）。它们可能是由污染等原因造成的。

**6. Peak calling，找峰**

   Call peak有几种不同的策略，可根据目的需要选择：

```
##转录因子的尖峰策略，参数含义用代码macs2 callpeak -h查看

macs2 callpeak -t AAA.sort.rmdup.bam -c input.sort.rmdup.bam -f BAMPE -g hs --SPMR  -B -q 0.01 --outdir ./ -n AAA --cutoff-analysis

##组蛋白的宽峰策略

macs2 callpeak -t AAA.sort.rmdup.bam -c input.sort.rmdup.bam -f BAMPE  -g hs --SPMR --broad --broad-cutoff  0.01 --outdir ./ -n AAA --cutoff-analysis

```
   **对于双端测序一定要选择BAMPE模式**，很过旧的分析方法里都没提到过这个问题，如果是CUT&Tag或者ATAC-seq，参考文件“-c”可以不设置。--SPMR 是MACS2自带的均一化方法，对峰弱化很多，如果只是peakcalling，可以不设置，而且可以添加```--keep-dup all```参数，添加blacklist。
   
   Call peak之后会得到生成“AAA_summit.bed”、“AAA_peaks.narrowPeak”、“AAA_control_lambda.bdg” 、“AAA_treat_pileup.bdg”和“AAA_peaks.xls”五个文件。



$ 

