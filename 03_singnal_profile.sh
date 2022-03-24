#----------------------------------------------------------------
#绘制每个样本peaks区域的信号分布图,热图(TSS 上下3kb 以及gene body区域分布)
#----------------------------------------------------------------
#创建文件夹以及传递参数


sample=$1
color="#f03e3e"  # #1c7ed6
bw_dir=08_bw

cd ${bw_dir}
mkdir heatmap_profiles_plot
cd ..


narrowPeak_dir=07_macs2_peaks_calling
outdir=08_bw/heatmap_profiles_plot
#----------------------------------------------------------------


#绘制单独每个样本的 peaks的 profiles 
computeMatrix reference-point \
    --referencePoint TSS \
    -R ${narrowPeak_dir}/${sample}_peaks.narrowPeak \
    -S ${bw_dir}/${sample}.bw \
    -b 3000 -a 3000 \
    --skipZeros \
    -p 16 \
    -o $outdir/matrix_${sample}_point.gz


plotProfile -m $outdir/matrix_${sample}_point.gz \
            -o $outdir/${sample}_TSS_profile.pdf \
            --dpi 750 \
            --plotWidth 10 \
            --plotHeight 8 \
            --startLabel TSS \
            --colors ${color}

plotHeatmap -m $outdir/matrix_${sample}_point.gz \
     -out $outdir/${sample}_TSS_heatmap.pdf \
     --colorList "white,${color} " \
     --whatToShow 'heatmap and colorbar' \
     --dpi 750 \
     --boxAroundHeatmaps yes \
     --averageTypeSummaryPlot mean \
     --missingDataColor white
#--heatmapHeight 30 \
#--zMin 20 20 \
#--zMax 70 70 \
#--alpha 0.8 \
#--colorMap "Blues" 
#绘制单独每个样本的 peaks的 heatmap

#----------------------------------------------------------------

#绘制gene Body区域的分布图
computeMatrix scale-regions \
    -R ${narrowPeak_dir}/${sample}_peaks.narrowPeak \
    -S ${bw_dir}/${sample}.bw \
    -b 3000 -a 3000 \
    --skipZeros \
    -p 16 \
    -o $outdir/matrix_${sample}_scaled.gz

plotProfile -m $outdir/matrix_${sample}_scaled.gz \
            -o $outdir/${sample}_Body_profile.pdf \
            --dpi 750 \
            --plotWidth 10 \
            --plotHeight 8 \
            --startLabel TSS \
            --colors ${color}

plotHeatmap -m $outdir/matrix_${sample}_scaled.gz \
            -o $outdir/${sample}_Body_heatmap.pdf \
            --colorList "white,${color} " \
            --whatToShow 'heatmap and colorbar' \
            --dpi 750 \
            --boxAroundHeatmaps yes \
            --heatmapHeight 30 \
            --averageTypeSummaryPlot mean \
            --missingDataColor white
#--zMin 20 20 \
#--zMax 70 70 \