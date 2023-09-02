library(tidyverse)
library(edgeR)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(readxl)
library(openxlsx)
library(ggpubr)
library(patchwork)
library(ggsci)
library(rtracklayer)
library(peakimport) #devtools::install_github('CharlesJB/peakimport')

#Read integration site----
integration<-read_xlsx('TSS.xlsx')%>%
  mutate(Orientation=ifelse(Orientation=='Same', 1, -1))%>%
  arrange(PID)
integration<-add_column(integration, Entry=c(1:dim(integration)[1]))%>%
  relocate(Entry)%>%
  mutate(Class=factor(Class, levels=c('Producer', 'Non-Producer', 'Defective')))

integration%>%mutate(Entry2=paste0(PID, Chr, Loc, Seq))%>%
  filter(duplicated(Entry2)) #Make sure there is no duplications

#TSS using nearest TSS from edgeR package----
#Scanning up and downstream TSS with same and opposite orientation
nafunction<-function(x, s){
  if(dim(x)[1]==0){
    x<-tibble(symbol=0, strand=s, tss=1E10)
  } else x<-x
}
i=1
tss<-tibble()
tmp_final<-NULL
tmp<-NULL
tmp_pos_pos<-NULL
tmp_pos_neg<-NULL
tmp_neg_pos<-NULL
tmp_neg_neg<-NULL
for (i in 1:dim(integration)[1]){
  tmp<-
    nearestTSS(chr=rep(integration$Chr[i], 2001), 
               locus = seq(integration$Loc[i]-1E6, integration$Loc[i]+1E6, by=1E3), species = 'Hs')%>%
    as_tibble()%>%dplyr::select(-distance)%>%unique()%>%
    mutate(Distance=integration$Loc[i]-tss)%>%
    arrange(abs(Distance))%>%
    mutate(strand=as.numeric(paste0(strand, '1')))
  tmp_pos_pos<-tmp%>%filter(strand==1, Distance>0)%>%
    dplyr::select(c(symbol, strand, tss))%>%
    dplyr::slice(1)%>%nafunction(., s=1)
  tmp_pos_neg<-tmp%>%filter(strand==1, Distance<0)%>%
    dplyr::select(c(symbol, strand, tss))%>%
    dplyr::slice(1)%>%nafunction(., s=1)
  tmp_neg_pos<-tmp%>%filter(strand==-1, Distance>0)%>%
    dplyr::select(c(symbol, strand, tss))%>%
    dplyr::slice(1)%>%nafunction(., s=-1)
  tmp_neg_neg<-tmp%>%filter(strand==-1, Distance<0)%>%
    dplyr::select(c(symbol, strand, tss))%>%
    dplyr::slice(1)%>%nafunction(., s=-1)
  tmp_final<-rbind(tmp_pos_pos, tmp_pos_neg, tmp_neg_pos, tmp_neg_neg)%>%
    add_column(PID=integration$PID[i], .before = 'symbol')%>%
    add_column(Entry=i, .before='PID')%>%
    mutate(stream=c(1, -1, 1, -1))%>%relocate(stream, .before = 'tss')
  tss<-rbind(tss, tmp_final)
}
View(tss)

rm(tmp_final, tmp, tmp_pos_neg, tmp_pos_pos, tmp_neg_neg, tmp_neg_pos)

data_tss<-left_join(tss, integration, by=c('Entry', 'PID'))%>%
  mutate(strand=strand*Orientation, 
         stream=stream*Orientation)%>%
  mutate(Distance=abs(Loc-tss))%>%
  mutate(Class=factor(Class, levels=c('Producer', 'Non-Producer', 'Defective')))%>%
  mutate(Distance=ifelse(Distance>1E8, 1E8, Distance))

cairo_pdf('1. TSS_nearestTSS.pdf', width=15, height=10)
ggarrange(
  data_tss%>%group_by(Entry)%>%summarise(min=min(Distance), Class=Class)%>%
    dplyr::select(Class, min)%>%unique()%>%
    ggplot(aes(x=Class, y=min, color=Class))+
    geom_jitter(width = 0.2)+
    scale_y_log10()+
    geom_boxplot(alpha=0.2, outlier.shape = NA, width=0.4)+
    stat_compare_means()+
    scale_color_nejm()+
    theme_classic()+guides(color='none')+xlab('')+
    ylab('Distance')+
    ggtitle('Nearest TSS'), 
  
  data_tss%>%filter(strand==1, stream==1)%>%
  ggplot(aes(x=Class, y=Distance, color=Class))+
    geom_jitter(width = 0.2)+
    scale_y_log10()+
    geom_boxplot(alpha=0.2, outlier.shape = NA, width=0.4)+
    stat_compare_means()+
    scale_color_nejm()+
    theme_classic()+guides(color='none')+xlab('')+
  ggtitle('Upstream, Same direction'),
  
  data_tss%>%filter(strand==1, stream==-1)%>%
  ggplot(aes(x=Class, y=Distance, color=Class))+
  geom_jitter()+
  scale_y_log10()+
  geom_boxplot(alpha=0.2, outlier.shape = NA, width=0.4)+
  stat_compare_means()+scale_color_nejm()+
  theme_classic()+guides(color='none')+xlab('')+
  ggtitle('Downstream, Same direction'),

  data_tss%>%filter(strand==-1, stream==1)%>%
  ggplot(aes(x=Class, y=Distance, color=Class))+
  geom_jitter()+
  scale_y_log10()+
  geom_boxplot(alpha=0.2, outlier.shape = NA, width=0.4)+
  stat_compare_means()+scale_color_nejm()+
  theme_classic()+guides(color='none')+xlab('')+
  ggtitle('Upstream, Opposite direction'),

  data_tss%>%filter(strand==-1, stream==-1)%>%
  ggplot(aes(x=Class, y=Distance, color=Class))+
  geom_jitter()+
  scale_y_log10()+
  geom_boxplot(alpha=0.2, outlier.shape = NA, width=0.4)+
  stat_compare_means()+scale_color_nejm()+
  theme_classic()+guides(color='none')+xlab('')+
  ggtitle('Downstream, Opposite direction'),

common.legend=TRUE, legend = 'none')
dev.off()

#TSS by refTSS, another database----
reftss<-import.bed('refTSS_v3.1_human_coordinate.hg38.bed')
reftss<-
  reftss%>%as_tibble()%>%dplyr::select(seqnames, start, end, strand, name)
reftss_anno<-read_delim('refTSS_v3.2_human_annotation.txt')%>%
  dplyr::select(`#CAGE_Peak_ID`, Gene_symbol)%>%
  dplyr::rename('name'='#CAGE_Peak_ID')
reftss<-left_join(reftss, reftss_anno, by='name')
reftss$Gene_symbol[is.na(reftss$Gene_symbol)]<-
  reftss$name[is.na(reftss$Gene_symbol)]

nafunction<-function(x, s){
  if(dim(x)[1]==0){
    x<-tibble(Gene_symbol=0, strand=s, start=1E12)
  } else x<-x
}

i=1
tss<-tibble()
tmp_final<-NULL
tmp<-NULL
tmp_pos_pos<-NULL
tmp_pos_neg<-NULL
tmp_neg_pos<-NULL
tmp_neg_neg<-NULL

for (i in 1:dim(integration)[1]){
  tmp<-
    reftss%>%filter(seqnames==paste0('chr', integration$Chr[i]))%>%
    filter(end>=integration$Loc[i]-1E6)%>%
    filter(start<=integration$Loc[i]+1E6)%>%
    mutate(strand=as.numeric(paste0(strand, '1')))%>%
    mutate(Distance=integration$Loc[i]-0.5*(start+end))%>%
    arrange(abs(Distance))
  
  if(dim(tmp)[1]==0) {
    tmp<-reftss%>%filter(seqnames==paste0('chr', integration$Chr[i]))%>%
      filter(end>=integration$Loc[i]-1E8)%>%
      filter(start<=integration$Loc[i]+1E8)%>%
      mutate(strand=as.numeric(paste0(strand, '1')))%>%
      mutate(Distance=integration$Loc[i]-0.5*(start+end))%>%
      arrange(abs(Distance))
  }
  
  tmp_pos_pos<-tmp%>%filter(strand==1, Distance>0)%>%
    dplyr::select(c(Gene_symbol, strand, start))%>%
    dplyr::slice(1)%>%nafunction(., s=1)
  tmp_pos_neg<-tmp%>%filter(strand==1, Distance<0)%>%
    dplyr::select(c(Gene_symbol, strand, start))%>%
    dplyr::slice(1)%>%nafunction(., s=1)
  tmp_neg_pos<-tmp%>%filter(strand==-1, Distance>0)%>%
    dplyr::select(c(Gene_symbol, strand, start))%>%
    dplyr::slice(1)%>%nafunction(., s=-1)
  tmp_neg_neg<-tmp%>%filter(strand==-1, Distance<0)%>%
    dplyr::select(c(Gene_symbol, strand, start))%>%
    dplyr::slice(1)%>%nafunction(., s=-1)
  tmp_final<-rbind(tmp_pos_pos, tmp_pos_neg, tmp_neg_pos, tmp_neg_neg)%>%
    add_column(PID=integration$PID[i], .before = 'Gene_symbol')%>%
    add_column(Entry=i, .before='PID')%>%
    mutate(stream=c(1, -1, 1, -1))%>%relocate(stream, .before = 'start')
  tss<-rbind(tss, tmp_final)
}
View(tss)

rm(tmp_final, tmp, tmp_pos_neg, tmp_pos_pos, tmp_neg_neg, tmp_neg_pos)

data_reftss<-left_join(tss, integration, by=c('Entry', 'PID'))%>%
  mutate(strand=strand*Orientation, 
         stream=stream*Orientation)%>%
  mutate(Distance=abs(Loc-start))%>%
  mutate(Class=factor(Class, levels=c('Producer', 'Non-Producer', 'Defective')))%>%
  mutate(Distance=ifelse(Distance>1E8, 1E8, Distance))

cairo_pdf('1. TSS_RefTSS.pdf', width=15, height=10)
ggarrange(
  data_reftss%>%group_by(Entry)%>%summarise(min=min(Distance), Class=Class)%>%
    dplyr::select(Class, min)%>%unique()%>%
    ggplot(aes(x=Class, y=min, color=Class))+
    geom_jitter(width = 0.2)+
    scale_y_log10()+
    geom_boxplot(alpha=0.2, outlier.shape = NA, width=0.4)+
    stat_compare_means()+
    scale_color_nejm()+
    theme_classic()+guides(color='none')+xlab('')+
    ylab('Distance')+
    ggtitle('Nearest TSS'), 
  
  data_reftss%>%filter(strand==1, stream==1)%>%
    ggplot(aes(x=Class, y=Distance, color=Class))+
    geom_jitter(width = 0.2)+
    scale_y_log10()+
    geom_boxplot(alpha=0.2, outlier.shape = NA, width=0.4)+
    stat_compare_means()+
    scale_color_nejm()+
    theme_classic()+guides(color='none')+xlab('')+
    ggtitle('Upstream, Same direction'),
  
  data_reftss%>%filter(strand==1, stream==-1)%>%
    ggplot(aes(x=Class, y=Distance, color=Class))+
    geom_jitter()+
    scale_y_log10()+
    geom_boxplot(alpha=0.2, outlier.shape = NA, width=0.4)+
    stat_compare_means()+scale_color_nejm()+
    theme_classic()+guides(color='none')+xlab('')+
    ggtitle('Downstream, Same direction'),
  
  data_reftss%>%filter(strand==-1, stream==1)%>%
    ggplot(aes(x=Class, y=Distance, color=Class))+
    geom_jitter()+
    scale_y_log10()+
    geom_boxplot(alpha=0.2, outlier.shape = NA, width=0.4)+
    stat_compare_means()+scale_color_nejm()+
    theme_classic()+guides(color='none')+xlab('')+
    ggtitle('Upstream, Opposite direction'),
  
  data_reftss%>%filter(strand==-1, stream==-1)%>%
    ggplot(aes(x=Class, y=Distance, color=Class))+
    geom_jitter()+
    scale_y_log10()+
    geom_boxplot(alpha=0.2, outlier.shape = NA, width=0.4)+
    stat_compare_means()+scale_color_nejm()+
    theme_classic()+guides(color='none')+xlab('')+
    ggtitle('Downstream, Opposite direction'),
  
  common.legend=TRUE, legend = 'none')
dev.off()

#Peaks----
chainObject <- import.chain("hg19ToHg38.over.chain") #Because Roadmap used hg19 to map

peaksearch<-function(x){
  peak<-import_peaks(paste0('Chip/E043-', x, '.narrowPeak')) #This comes from https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/ E043-series
  peak<-liftOver(peak, chainObject)  #Import peak file and lift them to hg38
  
  peak_tidy<-peak%>%as_tibble()%>%
    filter(qValue>-log10(0.05)) #change them to tidy format
  
  i=1
  peak_peaks<-tibble()
  tmp<-NULL
  tmp1<-NULL
  
  for (i in 1:dim(integration)[1]){
    tmp<-integration[i,]
    tmp1<-peak_tidy%>%filter(seqnames==paste0('chr', tmp$Chr))%>%
      filter(end>=(tmp$Loc-1E4), 
             start<=(tmp$Loc+1E4))%>% #Find the overlap
      mutate(Width=end-start, 
             Distance1=abs(0.5*(end+start)-tmp$Loc))%>%
      summarise(Signal=sum(Width*signalValue),  #Signal means width of peaks * mean signal values
                Peaks=length(signalValue),   #Count number of peaks
                Distance=min(Distance1))%>% #Count nearest peak distance to integration site
      mutate(Entry=i, PID=tmp$PID)%>%relocate(Entry, PID) 
    peak_peaks<-rbind(peak_peaks, tmp1)
  }
  peak_peaks$Signal<-replace_na(peak_peaks$Signal, 0)
  colnames(peak_peaks)[3]<-paste(x, 'Signal', sep='_')
  colnames(peak_peaks)[4]<-paste(x, 'Number', sep='_')
  colnames(peak_peaks)[5]<-paste(x, 'Distance', sep='_')
  return(peak_peaks)
}

data_histone<-integration%>%
  left_join(., peaksearch('H3K4me1'), by=c('Entry', 'PID'))%>%
  left_join(., peaksearch('H3K4me3'), by=c('Entry', 'PID'))%>%
  left_join(., peaksearch('H3K9me3'), by=c('Entry', 'PID'))%>%
  left_join(., peaksearch('H3K27ac'), by=c('Entry', 'PID'))%>%
  left_join(., peaksearch('H3K27me3'), by=c('Entry', 'PID'))%>%
  left_join(., peaksearch('H3K36me3'), by=c('Entry', 'PID')) 
data_histone$Class<-factor(data_histone$Class, levels=c('Producer', 
                                                        'Non-Producer', 
                                                        'Defective')) #Generate peak files of different histone markers
data_histone[data_histone==Inf]<-1E6
data_histone[str_detect(colnames(data_histone), 'Distance')]<-
  data_histone[str_detect(colnames(data_histone), 'Distance')]%>%
  lapply(., replace_na, replace=1E6)

i=1
list=NULL
tmp=NULL
tmplab=NULL
data_histone_1<-NULL
for (i in 1:(dim(data_histone)[2]-8)) {
  tmp<-colnames(data_histone)[i+8]
  tmplab=tmp
  data_histone_1<-data_histone
  if(max(data_histone_1[i+8])>=20) {
    data_histone_1[i+8]<-log10(data_histone[i+8]+1)
    tmplab=paste0(tmp, '(log10)')
  }
  list[[i]]<-
    data_histone_1%>%
    #filter(Class!='Defective')%>%
    ggplot(aes_string(x='Class', y=tmp, color='Class'))+
    geom_jitter(width = 0.2, size=0.5)+
    geom_boxplot(outlier.shape = NA, alpha=0.2, width=0.2, lwd=0.4)+
    theme_classic()+
    stat_compare_means()+
    ylab(tmplab)+
    scale_color_nejm(alpha=0.8)+
    guides(color='none')+xlab('') + ylim(c(-1, max(data_histone_1[i+8]+2)))
  
}

pdf('2. CHIP With defective.pdf', height=6, width=10)
i=1 #Signal
(list[[i]]|list[[i+3]]|list[[i+6]]) /
  (list[[i+9]]|list[[i+12]]|list[[i+15]]) 
i=2 #Number
(list[[i]]|list[[i+3]]|list[[i+6]]) /
  (list[[i+9]]|list[[i+12]]|list[[i+15]])
i=3 #Distance
(list[[i]]|list[[i+3]]|list[[i+6]]) /
  (list[[i+9]]|list[[i+12]]|list[[i+15]])
dev.off()

for (i in 1:(dim(data_histone)[2]-8)) {
  tmp<-colnames(data_histone)[i+8]
  tmplab=tmp
  data_histone_1<-data_histone
  if(max(data_histone_1[i+8])>=20) {
    data_histone_1[i+8]<-log10(data_histone[i+8]+1)
    tmplab=paste0(tmp, '(log10)')
  }
  list[[i]]<-
    data_histone_1%>%
    filter(Class!='Defective')%>%
    ggplot(aes_string(x='Class', y=tmp, color='Class'))+
    geom_jitter(width = 0.2, size=0.5)+
    geom_boxplot(outlier.shape = NA, alpha=0.2, width=0.2, lwd=0.4)+
    theme_classic()+
    stat_compare_means()+
    ylab(tmplab)+
    scale_color_nejm(alpha=0.8)+
    guides(color='none')+xlab('') + ylim(c(-1, max(data_histone_1[i+8]+2)))
  
}
pdf('2. CHIP Without defective.pdf', height=6, width=10)
i=1 #Signal
(list[[i]]|list[[i+3]]|list[[i+6]]) /
  (list[[i+9]]|list[[i+12]]|list[[i+15]]) 
i=2 #Number
(list[[i]]|list[[i+3]]|list[[i+6]]) /
  (list[[i+9]]|list[[i+12]]|list[[i+15]])
i=3 #Distance
(list[[i]]|list[[i+3]]|list[[i+6]]) /
  (list[[i+9]]|list[[i+12]]|list[[i+15]])
dev.off()

#Save files----
save.image('TSS.rdata')
