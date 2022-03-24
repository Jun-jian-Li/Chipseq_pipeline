configfile: "env/config.yaml"
import os
import re

print("双端测序有Input对照样本的命名格式为:Sample_Input_R1.fq.gz/Sample_IP_R1.fq.gz")
print("双端测序有Input对照样本的命名格式为:Sample_Input_R2.fq.gz/Sample_IP_R2.fq.gz")

samples = []
reads   = []
groups  = []

for item in os.listdir(os.path.abspath(".")+"/"+config["dir"]["raw_data_dir"]+"/"):
     #输出指定后缀类型的文件
     if(item.endswith('gz')):
         m=re.findall(r'(.*)_(.*)_(.*\d).*q.gz',item)       #findall匹配返回一个list
         sample_name_list = [x[0] for x in m]          #提取样本
         sample_name_list="".join(sample_name_list)    #list 转 str
         samples.append(sample_name_list)              #样本加入list
         samples=list(set(samples))                    #去除重复样本名字
         samples = [i for i in samples if(len(str(i))!=0)]#去除空

         #提取Group信息列表
         groups_name_list = [x[1] for x in m]            #提取read名
         groups_name_list="".join(groups_name_list)        #list 转 str
         groups.append(groups_name_list)                  #样本加入list
         groups=list(set(groups))                        #去除重复样本名字
         groups = [i for i in groups if(len(str(i))!=0)] #去除空

         #提取Read信息列表
         read_name_list = [x[2] for x in m]            #提取read名
         read_name_list="".join(read_name_list)        #list 转 str
         reads.append(read_name_list)                  #样本加入list
         reads=list(set(reads))                        #去除重复样本名字
         reads = [i for i in reads if(len(str(i))!=0)] #去除空

print (samples)  #获取的样本的名字
print (groups)   #获取的样本的分组
print (reads)    #获取的样本的测序read方向

rule all:
     input:
          #qc
          raw_data_fastqc = expand(os.path.join( config["dir"]["raw_data_qc_result_dir"],"{sample}_{group}_{R}_fastqc.html" ),sample=samples,group=groups,R=reads),
          clean = expand(os.path.join( config["dir"]["clean_data_dir"], "{sample}_{group}_{trim_read}.fq.gz"),sample=samples,group=groups,trim_read=['R1_val_1','R2_val_2']),
          clean_data_fastqc1=expand(os.path.join( config["dir"]["clean_data_qc_result_dir"],"{sample}_{group}_R1_val_1_fastqc.html"),sample=samples,group=groups),
          clean_data_fastqc2=expand(os.path.join( config["dir"]["clean_data_qc_result_dir"],"{sample}_{group}_R2_val_2_fastqc.html"),sample=samples,group=groups),
          cut_base_data1 = expand(os.path.join(config["dir"]["cut_base_dir"],"{sample}_{group}_cut_R1.fq.gz"),sample=samples,group=groups),
          cut_base_data2 = expand(os.path.join(config["dir"]["cut_base_dir"],"{sample}_{group}_cut_R2.fq.gz"),sample=samples,group=groups),
          cut_data_fastqc1=expand(os.path.join(config["dir"]["cut_data_qc_result_dir"],"{sample}_{group}_cut_R1_fastqc.html"),sample=samples,group=groups),
          cut_data_fastqc2=expand(os.path.join(config["dir"]["cut_data_qc_result_dir"],"{sample}_{group}_cut_R2_fastqc.html"),sample=samples,group=groups),
        
          # align
          sam_file=expand(os.path.join(config["dir"]["sam_dir"], "{sample}_{group}.sam"),sample=samples,group=groups),
          sort_bam=expand(os.path.join(config["dir"]["sort_bam_dir"],"{sample}_{group}.sort.bam"),sample=samples,group=groups),
          sort_bam_dup=expand(os.path.join(config["dir"]["dup_dir"],"{sample}_{group}_dup.sort.bam"),sample=samples,group=groups),
          sort_bam_bai = expand(os.path.join(config["dir"]["dup_dir"],"{sample}_{group}_dup.sort.bam.bai"),sample=samples,group=groups),
          #peaks calling
          peaks = expand(os.path.join(config["dir"]["macs2_call_peak"],"{sample}_peaks.narrowPeak"),sample=samples),
          #bam to bw
          bw = expand(os.path.join(config["dir"]["bw_file"],"{sample}_{group}.bw"),sample=samples,group=groups)



# 01 raw data quality control
rule raw_fastqc:
    input:
         raw_data = os.path.join(config["dir"]["raw_data_dir"],"{sample}_{group}_{R}.fq.gz")
    output:
         raw_data_fastqc = os.path.join(config["dir"]["raw_data_qc_result_dir"],"{sample}_{group}_{R}_fastqc.html")
    params:
         raw_data_fastqc_dir = os.path.join(config["dir"]["raw_data_qc_result_dir"])
    message: """----------------- Quality check of raw data with Fastqc.--------------"""
    shell:
         "fastqc {input.raw_data} -o {params.raw_data_fastqc_dir} "

# 02 filter low-quality sequences
rule trim:
    input:
         raw_data1 = os.path.join(config["dir"]["raw_data_dir"],"{sample}_{group}_R1.fq.gz"),
         raw_data2 =os.path.join(config["dir"]["raw_data_dir"],"{sample}_{group}_R2.fq.gz")
    output:
         clean1 = os.path.join(config["dir"]["clean_data_dir"], "{sample}_{group}_R1_val_1.fq.gz"),
         clean2 = os.path.join(config["dir"]["clean_data_dir"], "{sample}_{group}_R2_val_2.fq.gz")
    message: """-------------------------- Trim. ------------------------------------"""
    params:
         quality = config["para"]["trim_quality"],   #碱基最小质量分数
         stringency = config["para"]["trim_stringency"], #设定可以忍受的前后adapter重叠的碱基数,默认为1
         length = config["para"]["trim_length"],    #设定输出reads长度阈值,小于设定值会被抛弃
         output_dir = config["dir"]["clean_data_dir"]
    shell:
         "trim_galore --gzip --phred33 -q {params.quality} -j 16\
         --stringency {params.stringency} --length {params.length} --paired {input[0]} {input[1]} -o {params.output_dir}"

# 03 clean data quality control
rule trim_fastqc:
     input:
          clean_data1 = os.path.join(config["dir"]["clean_data_dir"], "{sample}_{group}_R1_val_1.fq.gz"),
          clean_data2 = os.path.join(config["dir"]["clean_data_dir"], "{sample}_{group}_R2_val_2.fq.gz")
     output:
          clean_data_fastqc1=os.path.join(config["dir"]["clean_data_qc_result_dir"],"{sample}_{group}_R1_val_1_fastqc.html"),
          clean_data_fastqc2=os.path.join(config["dir"]["clean_data_qc_result_dir"],"{sample}_{group}_R2_val_2_fastqc.html")
     params:
          os.path.join(config["dir"]["clean_data_qc_result_dir"])
     message: """------------ Quality check of clean data with Fastqc.---------------"""
     shell:
          "fastqc {input.clean_data1} -o {params[0]} && fastqc {input.clean_data2} -o {params[0]}"

# 04 cut base 
rule cut_base:
     input:
          clean_data1 = os.path.join(config["dir"]["clean_data_dir"], "{sample}_{group}_R1_val_1.fq.gz"),
          clean_data2 = os.path.join(config["dir"]["clean_data_dir"], "{sample}_{group}_R2_val_2.fq.gz")
     output:
          cut_base_data1 = os.path.join(config["dir"]["cut_base_dir"],"{sample}_{group}_cut_R1.fq.gz"),
          cut_base_data2 = os.path.join(config["dir"]["cut_base_dir"],"{sample}_{group}_cut_R2.fq.gz")
     message: """----------------------- Cut base ---------------------------------------."""
     params:
         left_base = config["para"]["cut_left_base"],
         right_base = config["para"]["cut_right_base"]
     shell:
          "cutadapt -u {params.left_base} -U {params.right_base} -o {output[0]} -p {output[1]} {input[0]} {input[1]}"

# 05 cut fastqc
rule cut_fastqc:
     input:
          cut_base_data1 = os.path.join(config["dir"]["cut_base_dir"],"{sample}_{group}_cut_R1.fq.gz"),
          cut_base_data2 = os.path.join(config["dir"]["cut_base_dir"],"{sample}_{group}_cut_R2.fq.gz")
     output:
          cut_data_fastqc1=os.path.join(config["dir"]["cut_data_qc_result_dir"],"{sample}_{group}_cut_R1_fastqc.html"),
          cut_data_fastqc2=os.path.join(config["dir"]["cut_data_qc_result_dir"],"{sample}_{group}_cut_R2_fastqc.html")
     params:
          os.path.join(config["dir"]["cut_data_qc_result_dir"])
     message: """------------ Cut base data with Fastqc.---------------"""
     shell:
          "fastqc {input.cut_base_data1} -o {params[0]} && fastqc {input.cut_base_data2} -o {params[0]} "


#--------------------------------------------------------------------------------------------------------------------------------
# 06 align to reference and sort bam file
rule align_sam:
     input:
         cut_base_data1 = os.path.join(config["dir"]["cut_base_dir"], "{sample}_{group}_cut_R1.fq.gz"),
         cut_base_data2 = os.path.join(config["dir"]["cut_base_dir"], "{sample}_{group}_cut_R2.fq.gz")
     output:
         sam_file = os.path.join(config["dir"]["sam_dir"], "{sample}_{group}.sam")
     log:
         align_log = os.path.join(config["dir"]["align_qc_dir"], "{sample}_{group}.align.log")
     message: """---------------------- Alignment ---------------------------------------."""
     params:
         reference = config["reference"]["huamn_reference"]
     shell:
         "bowtie2 -p 12 -x {params.reference} -1 {input.cut_base_data1} -2 {input.cut_base_data2} -S {output.sam_file} > {log.align_log} 2>&1"
        
# 07 rm low quality reads and rm chrM
rule rm_reads:
     input:
          sam_file = os.path.join(config["dir"]["sam_dir"], "{sample}_{group}.sam")
     output:
          sort_bam = os.path.join(config["dir"]["sort_bam_dir"],"{sample}_{group}.sort.bam")
     params:
          bam_qc = config["para"]["bam_quality"]
     shell:
          "samtools view -h -f 2 -@ 12 {input.sam_file} |grep -v 'chrM' |samtools sort -@ 12 - > {output.sort_bam}"

# 08 rm PCR duplication
rule rm_duplication:
     input:
          sort_bam = os.path.join(config["dir"]["sort_bam_dir"],"{sample}_{group}.sort.bam")
     output:
          sort_bam_dup = os.path.join(config["dir"]["dup_dir"],"{sample}_{group}_dup.sort.bam"),
          sort_bam_dup_csv = os.path.join(config["dir"]["dup_dir"],"{sample}_{group}_dup.csv")
     shell:
          "picard MarkDuplicates \
          I={input.sort_bam} \
          O={output.sort_bam_dup} \
          M={output.sort_bam_dup_csv} \
          REMOVE_DUPLICATES=true"

# 09 construct sort bam file index
rule sort_bam_index:
     input:
          sort_bam_dup = os.path.join(config["dir"]["dup_dir"],"{sample}_{group}_dup.sort.bam")
     output:
          sort_bam_bai = os.path.join(config["dir"]["dup_dir"],"{sample}_{group}_dup.sort.bam.bai")
     shell:
          "samtools index {input.sort_bam_dup} {output.sort_bam_bai}"

# 10 peaks calling 
rule peaks_calling:
     input:
          treat_bam_dup = os.path.join(config["dir"]["dup_dir"],"{sample}_IP_dup.sort.bam"),
          control_bam_dup = os.path.join(config["dir"]["dup_dir"],"{sample}_Input_dup.sort.bam")
     output:
          peaks = os.path.join(config["dir"]["macs2_call_peak"],"{sample}_peaks.narrowPeak")
     params:
          names = "{sample}",
          species = "hs"
     shell:
          "macs2 callpeak -t {input.treat_bam_dup} -c {input.control_bam_dup} -n {params.names} --nomodel -f BAMPE --keep-dup all -q 0.05 -B --SPMR -g {params.species} --outdir {output.peaks}"

# 11 bam transfor bw
rule bam_transfor_bw:
     input:
          sort_bam_dup = os.path.join(config["dir"]["dup_dir"],"{sample}_{group}_dup.sort.bam")
     output:
          bw = os.path.join(config["dir"]["bw_file"],"{sample}_{group}.bw")
     shell:
          "bamCoverage -p 16 \
            -b {input.sort_bam_dup} \
            --binSize 10 \
            --normalizeUsing RPKM \
            -o {output.bw}"













