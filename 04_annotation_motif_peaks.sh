#----------------------------------------------------------------
#对每个样本peaks区域的信号进行基因组注释
#----------------------------------------------------------------
#创建文件夹以及传递参数

mkdir 09_anno_motif

sample=$1

narrowPeak_dir=07_macs2_peaks_calling

mkdir 09_anno_motif/01_anno
mkdir 09_anno_motif/02_motif

anno_outdir=09_anno_motif/01_anno
motif_outdir=09_anno_motif/02_motif

annotatePeaks.pl ${narrowPeak_dir}/${sample}_peaks.narrowPeak mm10 > ${anno_outdir}/${sample}_peaks_anno.txt


# 提取peaks的位置信息文件
awk '{print $4"\t"$1"\t"$2"\t"$3"\t+"}' ${narrowPeak_dir}/${sample}_peaks.narrowPeak > ${motif_outdir}/${sample}_homer_peaks.tmp
findMotifsGenome.pl ${motif_outdir}/${sample}_homer_peaks.tmp mm10 ${motif_outdir} -len 8,10,12

