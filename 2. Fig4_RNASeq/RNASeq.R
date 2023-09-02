library(tidyverse)
library(pheatmap)
library(readxl)
library(openxlsx)
library(RColorBrewer)
library(factoextra)
library(DESeq2)
library(fgsea)
library(ggsci)
library(MetBrewer)
library(ggpubr)
library(M3C)
library(patchwork)

control<-read_delim('Control_matrix.txt')%>%as_tibble()

llv<-read_delim('NSV_matrix.txt')%>%as_tibble()

group<-read_xlsx('Group.xlsx')%>%na.exclude()

counts<-inner_join(llv, control, by='Gene_ID')%>%
  dplyr::select(Gene_ID, group$Sample)%>%
  dplyr::filter(str_starts(Gene_ID, '(RNA[:digit:]|RNR\\d)', negate = T)) %>%
  column_to_rownames('Gene_ID')%>%
  dplyr::filter(rowSums(across(.cols=everything()))>13)

#DEG----
dds<-DESeqDataSetFromMatrix(counts, DataFrame(group), ~ Group)
dds <- DESeq(dds, test = 'Wald') 
res<- results(dds, tidy = T)%>%as_tibble()
results(dds)%>%summary()
write.xlsx(results(dds, format = 'DataFrame', tidy = TRUE), 'DEG.xlsx')

#Volcano----
#Volcano plot
res%>%
  na.exclude()%>%
  mutate(color=ifelse(padj<0.05 & abs(log2FoldChange)>0.5, 'Sig', 'NS'))%>%
  mutate(shape=ifelse((str_starts(row, '(MIR|CERNA|VT|SNOR|XIST|LOC\\d|SCARNA|TALAM1|FAM239B|RNU|RPL13AP6|SUMO1P3)') | 
                         str_detect(row, '-AS|-IT|\\dP\\d|\\dAP\\d|PTENP1|TMEM198B|NBR2|CROCCP2|SCART1|LINC\\d|FAAHP1|CYP4F35P|TALAM1|KCNQ1OT1')), 'NC', 'C'))%>%
  View()

genelabel<-res%>%
  filter(padj<0.1, abs(log2FoldChange)>0.5)%>%
  filter(str_detect(string = row, pattern = '(CXC|ADAMTS1|MLKL|MTRNR2L2|FUT|BTN3A1|LTB|CIITA|SEMA4A|FCGR|NFKBID|BAG5|RCAN|FGF|EGF|TYK|TGF|BNIP|LTBP|PTGE|HDAC|PIK3|FRS2|WNT|RIPK3|CTS|CX\\dC|CCL|CCR|FOXP|ZAP|OAS1|IFI|IRF|IKB|IFN|TNF)') |
           str_starts(string=row, pattern = '(IL\\d|CD\\d|TNF|FCG|PARP|CTLA4|FAS|TRIM|SMAD|ADAM|FC|IFI|IRF|ITG|OAS|FAST|CC|FOX|PI|BCL|KLF|PFK)'))%>%
  filter(str_detect(row, '-AS')==FALSE)%>%
  select(row)

volcano<-
  res%>%
  na.exclude()%>%
  mutate(color=ifelse(padj<0.1 & pvalue<0.05 & padj>0.05 & abs(log2FoldChange)>0.5, 
                      'Padj 0.05-0.1', 
                      ifelse(padj<0.05 & abs(log2FoldChange)>0.5, 'Padj<0.05', 'NS')))%>%
  mutate(shape=ifelse((str_starts(row, 'MIR|FAM66E|SNOR|TALAM|VTRNA|XIST|LOC\\d|SCARNA|TALAM1|FAM239B|RNU|RPL13AP6|SUMO1P3') | 
                         str_detect(row, '-AS|-IT|RASA4DP|WHAMMP3|FAM27E3|FAM30C|FAM293B|\\dP\\d|\\dAP\\d|PTENP1|TMEM198B|NBR2|CROCCP2|SCART1|LINC\\d|FAAHP1|CYP4F35P|TALAM1|KCNQ1OT1')), 'NC', 'C'))%>%
  mutate(label=ifelse(row %in% genelabel$row, row, NA))%>%
  ggplot(aes(x=log2FoldChange , y=-log10(pvalue), color=color, shape=shape))+
  geom_point(size=1)+
  scale_color_manual(values = c('grey', 'dodgerblue2', pal_locuszoom()(4)[1]))+
  theme_bw()+
  geom_vline(xintercept =c(-0.5, 0, 0.5), linetype='dashed')+
  geom_hline(yintercept=-log10(0.05))+ylim(0, 10)+
  ggrepel::geom_label_repel(aes(label=label), size=2, color='black', max.overlaps = 20)

cairo_pdf('1. Volcano.pdf', width =6, height = 4)
volcano
dev.off()

#GSEA----
ranks<- res %>%  
  dplyr::select(row, stat)%>%deframe()

pathways.hallmark<-c(gmtPathways('GSEA/h.all.v2023.1.Hs.symbols.gmt'),
                     gmtPathways('GSEA/c2.cp.reactome.v2023.1.Hs.symbols.gmt'))
set.seed(111)
fgseaRes <- 
  fgseaMultilevel(pathways=pathways.hallmark, stats=ranks, nPermSimple = 10000)%>%
  as_tibble()%>%
  arrange(desc(NES))
write.xlsx(fgseaRes, 'GSEA.xlsx')

#Plotting
#Plot the normalized enrichment scores. Color the bar indicating whether or not the pathway was significant:
pos_NES<-
  fgseaRes %>% 
  dplyr::filter(padj<0.05, NES>1, 
               str_detect(pathway, 'INFECT|RUNX1|VIRUS|HIV_|INTERFERON|TNF|TLR|IL[:digit:]|IFN|OXID|MHC|APOPTOSIS|UNIQUIT|IMMUN|CYTOTOX|INFLAMM|PHOS')==TRUE)%>%
  dplyr::filter(str_detect(pathway, 'INNATE|SARS|INFLUENZA|KNOWN|HSCS')==FALSE)%>%
  mutate(pathway=str_to_title(str_replace_all(pathway, '_', ' ')))%>%
  mutate(pathway=str_replace_all(pathway, 'Hiv|hiv', 'HIV'))%>%
  mutate(pathway=str_replace_all(pathway, 'Mhc', 'MHC'))%>%
  mutate(pathway=str_replace_all(pathway, ' Of', ' of'))%>%
  ggplot(aes(reorder(pathway, NES), NES, fill=-log10(padj), size=size)) +
  geom_point(color='grey30', pch=21) +
  coord_flip() +
  labs(x="", y="Normalized Enrichment Score",
       title="Enriched in LLV") + 
  theme_minimal()+
  scale_size_continuous(limits = c(5,1000), breaks = c(10, 100, 1000))+
  scale_fill_gradient(low = pal_locuszoom(alpha=0.5)(5)[2], 
                       high= pal_locuszoom(alpha=1)(5)[1], limits=c(1.34, 10), 
                       na.value = 'firebrick4')  +
  ylim(1.5, 2.5)+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 50)) 
pos_NES
neg_NES<-
  fgseaRes %>% 
  dplyr::filter(padj<0.05,NES< -1,  
                str_detect(pathway, 'INFECT|HIV_|INFLAMM|VIRUS|INTERFERON|TNF|TLR|IL[:digit:]|IFN|OXID|MHC|APOPTOSIS|UNIQUIT|IMMUN|CYTOTOX')==TRUE)%>%
  dplyr::filter(str_detect(pathway, 'INNATE|SARS|INFLUENZA|LEISHM|IMMUNE_SYSTEM')==FALSE)%>%
  mutate(pathway=str_to_title(str_replace_all(pathway, '_', ' ')))%>%
  mutate(pathway=str_replace_all(pathway, 'Nfkb', 'NFKB'))%>%
  mutate(pathway=str_replace_all(pathway, 'Tnfa', 'TNF-alpha'))%>%
  mutate(pathway=str_replace_all(pathway, ' Of', ' of'))%>%
  mutate(pathway=str_replace_all(pathway, 'Ifns', ' IFNs'))%>%
  mutate(pathway=str_replace_all(pathway, 'Zbp1 Dai', ' ZBP1-DAI'))%>%
  ggplot(aes(reorder(pathway, NES), NES, fill=-log10(padj), size=size)) +
  geom_point(color='grey30',pch=21) +
  coord_flip() +
  labs(x="", y="Normalized Enrichment Score",
       title="Enriched in Control") + 
  theme_minimal()+
  scale_size_continuous(limits = c(5,1000), breaks = c(10, 100, 1000))+
  scale_fill_gradient(low = pal_locuszoom(alpha=0.4)(5)[3], 
                       high= pal_locuszoom(alpha=1)(5)[5], 
                       limits=c(1.34, 5))  + xlim(-2.5, 1.5)+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 50)) +ylim(-2.5, -1.5)

neg_NES

cairo_pdf('2. GSEA_LLV_vs_Control_2.pdf', width =6, height = 6)
pos_NES/neg_NES + plot_layout(heights = c(1.5, 1))
dev.off()

#Plot GSEA enrichment----
plotEnrichment_mod<-
  function (pathway, stats, gseaParam = 1, 
            note, yaxislim, Segment_size=0.1) {
    rnk <- rank(-stats)
    ord <- order(rnk)
    statsAdj <- stats[ord]
    statsAdj <- sign(statsAdj) * (abs(statsAdj)^gseaParam)
    statsAdj <- statsAdj/max(abs(statsAdj))
    pathway <- unname(as.vector(na.omit(match(pathway, names(statsAdj)))))
    pathway <- sort(pathway)
    gseaRes <- calcGseaStat(statsAdj, selectedStats = pathway, 
                            returnAllExtremes = TRUE)
    bottoms <- gseaRes$bottoms
    tops <- gseaRes$tops
    n <- length(statsAdj)
    xs <- as.vector(rbind(pathway - 1, pathway))
    ys <- as.vector(rbind(bottoms, tops))
    toPlot <- data.frame(x = c(0, xs, n + 1), y = c(0, ys, 0))
    diff <- (max(tops) - min(bottoms))/8
    x = y = NULL
    g <- ggplot(toPlot, aes(x = x, y = y)) + 
      geom_line(color = "grey80", linetype='dashed') + theme_pubclean() + 
      geom_point(color = "#357EBDCC", size = 0.8, alpha=0.4) + 
      geom_segment(data = data.frame(x = pathway), 
                   mapping = aes(x = x, y = -0.02, 
                                 xend = x, yend = 0.02), 
                   size=0.1, color='grey50', alpha=0.75) + 
      theme(panel.border = element_blank(), panel.grid.minor = element_blank()) + 
      labs(x = "Rank", y = "Enrichment score")+ 
      guides(color='none')+
      theme(plot.title = element_text(size = 10))+
      ggtitle(note)+
      ylim(yaxislim)
    g
  }

fgseaRes%>%
  filter(pathway %in% c('REACTOME_HIV_INFECTION', 
                        'REACTOME_APOPTOSIS', 
                        'REACTOME_INTERFERON_ALPHA_BETA_SIGNALING', 
                        'HALLMARK_INTERFERON_GAMMA_RESPONSE'))
pdf('3. Enrichment IFN.pdf', height = 4, width=12)
plotEnrichment_mod(pathway = pathways.hallmark$REACTOME_HIV_INFECTION, 
                   ranks, Segment_size = 0.001, yaxislim = c(-0.6, 0.6),
                   note = 'Reactome HIV Infection, NES=2.0, Padj=3.6E-7') |
  plotEnrichment_mod(pathway = pathways.hallmark$REACTOME_APOPTOSIS, 
                     ranks, Segment_size = 0.001,yaxislim = c(-0.6, 0.6),
                     note = 'Reactome Apoptosis, NES=1.7, Padj=0.002') |
  plotEnrichment_mod(pathway = pathways.hallmark$REACTOME_INTERFERON_ALPHA_BETA_SIGNALING, 
                     ranks, Segment_size = 0.001, yaxislim = c(-0.6, 0.6),
                     note = 'Reactome Interferon Alpha/Beta Signaling, NES=-1.8, Padj=0.01') |
  plotEnrichment_mod(pathway = pathways.hallmark$HALLMARK_INTERFERON_GAMMA_RESPONSE, 
                     ranks, Segment_size = 0.001,yaxislim = c(-0.6, 0.6),
                     note = 'Hallmark Interferon Gamma Response, NES=-1.6, Padj=0.006') 
  
dev.off()

#Save----
save.image('RNAseq.rdata')
