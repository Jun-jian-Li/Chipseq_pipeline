library(data.table)
library(tidyverse)
library(data.table)
library(stringr)
library(Hmisc)
library(ggplot2)
library(scales)
library(ggmap)   #为了引用主题theme_nothing，用来消除原始ggplot绘图自带的一切标签
library(gridExtra)


args <- commandArgs(T)

sample<-args[1]

path<-getwd()

peak<-fread(paste0(path,"/",sample,"_peaks_anno.txt"),header = T,sep="\t")

anno.df<-as.data.frame(table(str_split(peak$Annotation,"\\(",simplify = T)[,1]))
anno.df$Var1<-capitalize(as.character(anno.df$Var1))

colnames(anno.df)<-c("Group","value")
anno.df<-anno.df[order(anno.df$value,decreasing = T),]
anno.df$Group

df1<-anno.df %>% 
     as.data.frame() %>% 
     mutate(Group = factor(Group,levels = c("Intron ","Promoter-TSS ","Exon ","5' UTR ","Intergenic","Non-coding ","TTS ","3' UTR ")),
            cumulative = cumsum(value),
            midpoint = cumulative - value / 2,
            label = paste0(round(value / sum(value) * 100, 2), "%"))
                                             

#group是分类情况，value是所占的比例：百分之50的话值就是50

p1<-ggplot(df1, aes(x=1,y=value, fill=Group))+
      geom_bar(width = 1, stat = "identity",position = "stack")+#identity是直接引用数据集中变量的值（表示不要计数，而是直接使用数据本身作为频数。）
      coord_polar("y", start=0)+   #转换为极坐标
      # scale_fill_manual(values=c("#63cdda", "#f8a5c2", "#f5cd79"))+ #修改饼图的颜色（自定义颜色）
      scale_fill_brewer(palette="RdPu",direction = -1)+ #使用来自RColorBrewer的调色板
      #添加文字
      geom_text(aes(x = 2, y = midpoint, label = label),size =4, colour="Black")+ ## 加上百分比标签的位置和数值
      #添加主题    
      theme_nothing(legend = TRUE)+
      theme(legend.title=element_blank(),
            legend.text=element_text(size=12),
            legend.position = "right")

ggsave(p1,filename = paste0(path,"/",sample,"_peaks_anno.pdf"),height = 8,width = 8) 
