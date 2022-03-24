#!/bin/bash

#删除运行完的文件节约计算机储存资源
rm -rf 01_raw_data/*
rm -rf 02_trim_data/*
rm -rf 03_cut_base/*
rm -rf 04_sam/*
rm -rf 05_sort_bam/*
rm -rf 06_rm_dup_bam/*
rm -rf 07_macs2_peaks_calling/*
rm -rf 08_bw/*
rm -rf 09_anno_motif/*
rm -rf result/01_raw_data_qc/*
rm -rf result/02_clean_data_qc/*
rm -rf result/03_cut_data_qc/*
rm -rf result/04_bam_qc/*


echo "文件移除完成,记得进行下游分析,希望会有好的结果!!!"