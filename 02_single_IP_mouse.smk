configfile: "env/config.yaml"
import os
import re

print("单端测序只有IP样本的命名格式为:Sample_IP.fq.gz")

samples = []
groups  = []
for item in os.listdir(os.path.abspath(".")+"/"+config["dir"]["raw_data_dir"]+"/"):
     #输出指定后缀类型的文件
    if(item.endswith('gz')):
          m=re.findall(r'(.*)_(.*).f.*q.gz',item)               #findall匹配返回一个list
          sample_name_list = [x[0] for x in m]               #提取样本
          sample_name_list="".join(sample_name_list)         #list 转 str
          samples.append(sample_name_list)                   #样本加入list
          samples=list(set(samples))                         #去除重复样本名字
          samples = [i for i in samples if(len(str(i))!=0)]  #去除空

          #提取Group信息列表
          # groups_name_list = [x[1] for x in m]              #提取read名
          # groups_name_list="".join(groups_name_list)        #list 转 str
          # groups.append(groups_name_list)                   #样本加入list
          # groups=list(set(groups))                          #去除重复样本名字
          # groups = [i for i in groups if(len(str(i))!=0)]   #去除空

print(samples)

rule all:
     input:
          #qc
          raw_data_fastqc = expand(os.path.join( config["dir"]["raw_data_qc_result_dir"],"{sample}_IP_fastqc.html" ),sample=samples),
          clean_data = expand(os.path.join( config["dir"]["clean_data_dir"], "{sample}_IP_trimmed.fq.gz"),sample=samples),
          clean_data_fastqc=expand(os.path.join( config["dir"]["clean_data_qc_result_dir"],"{sample}_IP_trimmed_fastqc.html"),sample=samples),
          cut_base_data = expand(os.path.join(config["dir"]["cut_base_dir"],"{sample}_IP_cut.fq.gz"),sample=samples),
          cut_data_fastqc = expand(os.path.join(config["dir"]["cut_data_qc_result_dir"],"{sample}_IP_cut_fastqc.html"),sample=samples),
        
          # align
          sam_file=expand(os.path.join(config["dir"]["sam_dir"], "{sample}_IP.sam"),sample=samples),
          sort_bam=expand(os.path.join(config["dir"]["sort_bam_dir"],"{sample}_IP.sort.bam"),sample=samples),
          sort_bam_dup=expand(os.path.join(config["dir"]["dup_dir"],"{sample}_IP_dup.sort.bam"),sample=samples),
          sort_bam_bai = expand(os.path.join(config["dir"]["dup_dir"],"{sample}_IP_dup.sort.bam.bai"),sample=samples),
        
          # peaks calling
          peaks = expand(os.path.join(config["dir"]["macs2_call_peak"],"{sample}_peaks_narrowPeak"),sample=samples),
          
          #bam to bw
          bw = expand(os.path.join(config["dir"]["bw_file"],"{sample}_IP.bw"),sample=samples)


# 01 raw data quality control--------------------------------
rule raw_fastqc:
     input:
          raw_data = os.path.join(config["dir"]["raw_data_dir"],"{sample}_IP.fq.gz")
     output:
          raw_data_fastqc = os.path.join(config["dir"]["raw_data_qc_result_dir"],"{sample}_IP_fastqc.html")
     params:
          raw_data_fastqc_dir = os.path.join(config["dir"]["raw_data_qc_result_dir"])
     message: """----------------- Quality check of raw data with Fastqc.--------------"""
     shell: 
          "fastqc {input.raw_data} -o {params.raw_data_fastqc_dir} "

# 02 filter low-quality sequences--------------------------------
rule trim:
     input:
          raw_data = os.path.join(config["dir"]["raw_data_dir"],"{sample}_IP.fq.gz")
     output:
          clean_data = os.path.join(config["dir"]["clean_data_dir"],"{sample}_IP_trimmed.fq.gz")
     message: """-------------------------- Trim. ------------------------------------"""
     params:
          quality = config["para"]["trim_quality"],   #碱基最小质量分数
          stringency = config["para"]["trim_stringency"], #设定可以忍受的前后adapter重叠的碱基数,默认为1
          length = config["para"]["trim_length"],    #设定输出reads长度阈值,小于设定值会被抛弃
          output_dir = config["dir"]["clean_data_dir"]
     shell:
          "trim_galore --gzip --phred33 -q {params.quality} -j 12 --stringency {params.stringency} --length {params.length} {input[0]} -o {params.output_dir}"


# 03 clean data quality control--------------------------------
rule trim_fastqc:
     input:
          clean_data = os.path.join(config["dir"]["clean_data_dir"], "{sample}_IP_trimmed.fq.gz")
     output:
          clean_data_fastqc=os.path.join(config["dir"]["clean_data_qc_result_dir"],"{sample}_IP_trimmed_fastqc.html")         
     params:
          os.path.join(config["dir"]["clean_data_qc_result_dir"])
     message: """------------ Quality check of clean data with Fastqc.---------------"""
     shell:
          "fastqc {input.clean_data} -o {params[0]}"


# 04 cut base --------------------------------
rule cut_base:
     input:
          clean_data = os.path.join(config["dir"]["clean_data_dir"],"{sample}_IP_trimmed.fq.gz")
     output:
          cut_base_data = os.path.join(config["dir"]["cut_base_dir"],"{sample}_IP_cut.fq.gz")
     message: """----------------------- Cut base ---------------------------------------."""
     params:
          base = config["para"]["cut_left_base"]
     shell:
          "cutadapt -u {params.base} -o {output[0]} {input[0]}"

# 05 cut fastqc--------------------------------
rule cut_fastqc:
     input:
          cut_base_data = os.path.join(config["dir"]["cut_base_dir"],"{sample}_IP_cut.fq.gz")
     output:
          cut_data_fastqc=os.path.join(config["dir"]["cut_data_qc_result_dir"],"{sample}_IP_cut_fastqc.html")
     params:
          os.path.join(config["dir"]["cut_data_qc_result_dir"])
     message: """------------ Cut base data with Fastqc.---------------"""
     shell:
          "fastqc {input.cut_base_data} -o {params[0]}"


#--------------------------------------------------------------------------------------------------------------------------------
# 06 align to reference and sort bam file
rule align_sam:
     input:
          cut_base_data = os.path.join(config["dir"]["cut_base_dir"], "{sample}_IP_cut.fq.gz")
     output:
          sam_file = os.path.join(config["dir"]["sam_dir"], "{sample}_IP.sam")
     log:
          align_log = os.path.join(config["dir"]["align_qc_dir"], "{sample}_IP.align.log")
     message: """---------------------- Alignment ---------------------------------------."""
     params:
          reference = config["reference"]["mouse_reference"]
     shell:
          "bowtie2 -p 12 -x {params.reference} -U {input.cut_base_data} -S {output.sam_file} > {log.align_log} 2>&1"
        
# 07 rm low quality reads and rm chrM
rule rm_reads:
     input:
          sam_file = os.path.join(config["dir"]["sam_dir"], "{sample}_IP.sam")
     output:
          sort_bam = os.path.join(config["dir"]["sort_bam_dir"],"{sample}_IP.sort.bam")
     params:
          bam_qc = config["para"]["bam_quality"]
     shell:
          "samtools view -h -f 2 -@ 12 {input.sam_file} |grep -v 'chrM' |samtools sort -@ 12 - > {output.sort_bam}"

# 08 rm PCR duplication
rule rm_duplication:
     input:
          sort_bam = os.path.join(config["dir"]["sort_bam_dir"],"{sample}_IP.sort.bam")
     output:
          sort_bam_dup = os.path.join(config["dir"]["dup_dir"],"{sample}_IP_dup.sort.bam"),
          sort_bam_dup_csv = os.path.join(config["dir"]["dup_dir"],"{sample}_IP_dup.csv")
     shell:
          "picard MarkDuplicates \
          I={input.sort_bam} \
          O={output.sort_bam_dup} \
          M={output.sort_bam_dup_csv} \
          REMOVE_DUPLICATES=true"

# 09 construct sort bam file index
rule sort_bam_index:
     input:
          sort_bam_dup = os.path.join(config["dir"]["dup_dir"],"{sample}_IP_dup.sort.bam")
     output:
          sort_bam_bai = os.path.join(config["dir"]["dup_dir"],"{sample}_IP_dup.sort.bam.bai")
     shell:
          "samtools index {input.sort_bam_dup} {output.sort_bam_bai}"

# 10 peaks calling 
rule peaks_calling:
     input:
          treat_bam_dup = os.path.join(config["dir"]["dup_dir"],"{sample}_IP_dup.sort.bam")
     output:
          peaks = os.path.join(config["dir"]["macs2_call_peak"],"{sample}_peaks_narrowPeak")
     params:
          names = "{sample}",
          species = "mm"
     shell:
          "macs2 callpeak -t {input.treat_bam_dup} -n {params.names} --nomodel -f BAMPE --keep-dup all -q 0.05 -B --SPMR -g {params.species} --outdir {output.peaks}"

# 11 bam transfor bw
rule bam_transfor_bw:
     input:
          sort_bam_dup = os.path.join(config["dir"]["dup_dir"],"{sample}_IP_dup.sort.bam")
     output:
          bw = os.path.join(config["dir"]["bw_file"],"{sample}_IP.bw")
     shell:
          "bamCoverage -p 16 \
            -b {input.sort_bam_dup} \
            --binSize 10 \
            --normalizeUsing RPKM \
            -o {output.bw}"