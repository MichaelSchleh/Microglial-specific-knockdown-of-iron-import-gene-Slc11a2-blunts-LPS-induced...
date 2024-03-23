#Manuscript Title: "Microglial-specific knockdown of iron import gene, Slc11a2, blunts LPS-induced neuroinflammatory responses in a sex-specific manner
#Purpose: To determine the affect of LPS on microglia-specific SLC11a2 knockdown compared to WT floxed Control mice

#Citation: Volk Robertson K, Schleh MW, Harrison FE, Hasty AH. 
#          Microglial-specific knockdown of iron import gene, Slc11a2, blunts LPS-induced neuroinflammatory responses in a sex-specific manner. 
#          Brain Behav Immun. 2024 Feb;116:370-384. doi: 10.1016/j.bbi.2023.12.020. Epub 2023 Dec 22. 
#          PMID: 38141840; PMCID: PMC10874246.

#Creator: Michael Schleh (michael.w.schleh@vanderbilt.edu)

#Groups;
# 1) ControlSaline
# 2) ControlLPS
# 3) KDSaline
# 4) KDLPS

#Analysis Plan
  #1) ControlSaline vs. ControlLPS
  #2) KDSaline vs. KDLPS
  #2) Control LPS vs. KDLPS

#Packages used in lab pipelines
library(gplots)
library(dplyr)
library(tidyverse)
library(stringr)
library(calibrate)
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(vsn)
library(RColorBrewer)
library(edgeR)
library(textshape)
library(clusterProfiler)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(DOSE)
library(enrichplot)
library(EnhancedVolcano)
library(limma)
library(dplyr)
library(ggVennDiagram)
library(apeglm)
library(ggrepel)
library(AnnotationDbi)
library(ggpubr)
library(reshape2)

#Setwd and Figure Path
setwd("INPUT WORKING DIRECTORY")
path = "INPUT PATH TO SAVE FIGURES"

#Create folders in working directory
 dir.create("DEGs")
 dir.create("DEGs/UP")
 dir.create("DEGs/DOWN")
 dir.create("DEGs/DEG_Results")
 dir.create("Volcano_Plots")
 dir.create("MA_Plots")
 dir.create("VennDiagrams")
 dir.create("VennDiagrams/upregulated")
 dir.create("VennDiagrams/downregulated")
 dir.create("QCPlots")


## Load in Coldata files: Filter males & filter out KV5,KV10, KV11 from comparisons ----
Coldata <- rbind.data.frame(read.csv("MetaData", header = TRUE, sep = ","))
Coldata <- Coldata %>% filter(Sex=="Male" & !Sample %in% c("KV.5","KV.10","KV.11")) #FILTER OUT MALE MICE WITH <50% SLC11A2 KNOCKDOWN
Samples <- Coldata$Sample

## Load in Countdata files:

#Join together, 1) ControlSaline_vs_ControlLPS & 2) KDSaline_vs_KDLPS datasets. 
Countdata <- read.csv("RNA-Seq Countdata", header = TRUE).
Countdata <- Countdata[-2] #Removes the "Gene name" that comes with this upload

#Cleanup Countdata columns by Removes "X8233." from all columns, leaving only "KV" identifiers
Countdata <- data.frame(Countdata[1],
                        Countdata[-1] %>% 
                          rename_with(~ str_remove(.,"X8233.")))
Countdata <- Countdata %>% column_to_rownames("X")

#Filter the Countdata dataset by the Coldata samples
Countdata <- Countdata %>% dplyr::select(all_of(Samples))
Countdata <- as.matrix(Countdata)

saline <- factor(c(rep("Saline",4), rep("LPS",5), rep("Saline",7),rep("LPS",4)))
genotype <- factor(c(rep("Control",9),rep("KD",11)))
Coldata <- data.frame(row.names=colnames(Countdata),
                      #Coldata[3],
                      saline,
                      genotype)
#Run DESeq -----
dds <- DESeqDataSetFromMatrix(countData = Countdata, colData = Coldata, design=~ saline + genotype)
colnames(dds) <- Samples
dds$group <- factor(paste0(dds$saline, dds$genotype))
design(dds) <- ~ group
#Use For salinecontrol vs. LPScontrol comparison
dds$group <- relevel(dds$group, ref="SalineControl") 
#dds
dds <- DESeq(dds)

#Use resultsNames to select either 
  #1) "group_LPScontrol_vs_salinecontrol"
  #2) "group_LPSKD_vs_salineKD"
resultsNames(dds) # all comparison identifiers should be "_vs_salinecontrol"

#Set thresholds for differentially upregulated genes and make directories
padj_cutoff <- 0.05
fc_UP <- 1
fc_DOWN <- -1

#Get Differential Expression Results against the Saline reference control.
#1) "group_LPScontrol_vs_salinecontrol----

res_group_LPSControl_vs_SalineControl <- lfcShrink(dds,coef="group_LPSControl_vs_SalineControl",type="apeglm")
table(res_group_LPSControl_vs_SalineControl$padj<0.05)
res_group_LPSControl_vs_SalineControl <- res_group_LPSControl_vs_SalineControl[order(res_group_LPSControl_vs_SalineControl$padj), ]

#Annotate results object with gene symbols
res_group_LPSControl_vs_SalineControl$symbols  <- mapIds(org.Mm.eg.db, keys= rownames(res_group_LPSControl_vs_SalineControl), column= "SYMBOL", keytype= "ENSEMBL", multiVals="first")
res_group_LPSControl_vs_SalineControl <- data.frame(Ensemble=rownames(res_group_LPSControl_vs_SalineControl),res_group_LPSControl_vs_SalineControl)  

#Write DEG list to excel file
write.csv(res_group_LPSControl_vs_SalineControl, file="DEGs/DEresults_LPScontrol_vs_SalineControl_MALES.csv")

#Create files containing only upregulated DEGs
res_group_LPSControl_vs_SalineControl_DEGs <- res_group_LPSControl_vs_SalineControl %>% data.frame() %>% as_tibble() 
res_group_LPSControl_vs_SalineControl_DEGs_UP <- res_group_LPSControl_vs_SalineControl_DEGs %>% filter(padj < padj_cutoff) %>% filter(log2FoldChange > fc_UP)
write.csv(res_group_LPSControl_vs_SalineControl_DEGs_UP, file="DEGs/UP/DEGs_UP_LPScontrol_vs_SalineControl.csv")

#Create files containing only downregulated DEGs
res_group_LPSControl_vs_SalineControl_DEGs_DOWN <- res_group_LPSControl_vs_SalineControl_DEGs %>% filter(padj < padj_cutoff) %>% filter(log2FoldChange < fc_DOWN)
write.csv(res_group_LPSControl_vs_SalineControl_DEGs_DOWN, file="DEGs/DOWN/DEGs_DOWN_LPScontrol_vs_SalineControl.csv")

#Re-set dds to set comparisons against salineKD ----
#2) group_LPSKD_vs_salineKD 
dds$group <- relevel(dds$group, ref="SalineKD")
dds <- DESeq(dds)
resultsNames(dds) # all comparison identifiers should be "_vs_salinecontrol"


res_group_LPSKD_vs_SalineKD <- lfcShrink(dds,coef="group_LPSKD_vs_SalineKD",type="apeglm")
table(res_group_LPSKD_vs_SalineKD$padj<0.05)
res_group_LPSKD_vs_SalineKD <- res_group_LPSKD_vs_SalineKD[order(res_group_LPSKD_vs_SalineKD$padj), ]

#Annotate results object with gene symbols
res_group_LPSKD_vs_SalineKD$symbols  <- mapIds(org.Mm.eg.db, keys= rownames(res_group_LPSKD_vs_SalineKD), column= "SYMBOL", keytype= "ENSEMBL", multiVals="first")
res_group_LPSKD_vs_SalineKD <- data.frame(Ensemble=rownames(res_group_LPSKD_vs_SalineKD),res_group_LPSKD_vs_SalineKD)  

#Write DEG list to excel file
write.csv(res_group_LPSKD_vs_SalineKD, file="DEGs/DEresults_LPSKD_vs_SalineKD_MALES.csv")

#Create files containing only upregulated DEGs
res_group_LPSKD_vs_SalineKD_DEGs <- res_group_LPSKD_vs_SalineKD %>% data.frame() %>% as_tibble() 
res_group_LPSKD_vs_SalineKD_DEGs_UP <- res_group_LPSKD_vs_SalineKD_DEGs %>% filter(padj < padj_cutoff) %>% filter(log2FoldChange > fc_UP)
write.csv(res_group_LPSKD_vs_SalineKD_DEGs_UP, file="DEGs/UP/DEGs_UP_LPSKD_vs_SalineKD.csv")

#Create files containing only downregulated DEGs
res_group_LPSKD_vs_SalineKD_DEGs_DOWN <- res_group_LPSKD_vs_SalineKD_DEGs %>% filter(padj < padj_cutoff) %>% filter(log2FoldChange < fc_DOWN)
write.csv(res_group_LPSKD_vs_SalineKD_DEGs_DOWN, file="DEGs/DOWN/DEGs_DOWN_LPSKD_vs_SalineKD.csv")

#Re-set dds to set comparisons against salineKD ----
#3) group_KDLPS_vs_LPSKD
dds$group <- relevel(dds$group, ref="LPSKD")
dds <- DESeq(dds)
resultsNames(dds) # all comparison identifiers should be "_vs_salinecontrol"


res_group_LPSControl_vs_LPSKD <- lfcShrink(dds,coef="group_LPSControl_vs_LPSKD",type="apeglm")
table(res_group_LPSControl_vs_LPSKD$padj<0.05)
res_group_LPSControl_vs_LPSKD <- res_group_LPSControl_vs_LPSKD[order(res_group_LPSControl_vs_LPSKD$padj), ]

#Annotate results object with gene symbols
res_group_LPSControl_vs_LPSKD$symbols  <- mapIds(org.Mm.eg.db, keys= rownames(res_group_LPSControl_vs_LPSKD), column= "SYMBOL", keytype= "ENSEMBL", multiVals="first")
res_group_LPSControl_vs_LPSKD <- data.frame(Ensemble=rownames(res_group_LPSControl_vs_LPSKD),res_group_LPSControl_vs_LPSKD)  

#Write DEG list to excel file
write.csv(res_group_LPSControl_vs_LPSKD, file="DEGs/DEresults_LPSControl_vs_LPSKD_MALES.csv")

#Create files containing only upregulated DEGs
res_group_LPSControl_vs_LPSKD_DEGs <- res_group_LPSControl_vs_LPSKD %>% data.frame() %>% as_tibble() 
res_group_LPSControl_vs_LPSKD_DEGs_UP <- res_group_LPSControl_vs_LPSKD_DEGs %>% filter(padj < padj_cutoff) %>% filter(log2FoldChange > fc_UP)
write.csv(res_group_LPSControl_vs_LPSKD_DEGs_UP, file="DEGs/UP/DEGs_UP_LPSControl_vs_LPSKD.csv")

#Create files containing only downregulated DEGs
res_group_LPSControl_vs_LPSKD_DEGs_DOWN <- res_group_LPSControl_vs_LPSKD_DEGs %>% filter(padj < padj_cutoff) %>% filter(log2FoldChange < fc_DOWN)
write.csv(res_group_LPSControl_vs_LPSKD_DEGs_DOWN, file="DEGs/DOWN/DEGs_DOWN_LPSControl_vs_LPSKD.csv")

#Figure Generation

#Figure Panels for Male Data ----
#6B Principal components analysis----
vst <- vst(dds)
PCAdata1 <- plotPCA(vst, intgroup="group", returnData = TRUE)
percentVar <- round(100 * attr(PCAdata1, "percentVar")) 
PCAdata1$group <- factor(PCAdata1$group, levels=c("SalineControl","LPSControl","SalineKD","LPSKD"))


ggplot(PCAdata1, aes(x = PC1, y = PC2, color = group,label=rownames(PCAdata1))) + 
  stat_ellipse(geom = "polygon", alpha = 3/4, level=0.9, aes(fill = group,color=group),size=0.75) +
  geom_point(aes(fill=(group)), size =3,shape=21) + 
  scale_fill_manual(values=c("white","white","#ece6ca","#0f99b2")) + 
  scale_color_manual(values=c("#ece6ca","#0f99b2","black","black")) + 
  # scale_x_continuous(limits=c(-45,15)) +
  # scale_y_continuous(limits=c(-6,22)) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) + 
  ylab(paste0("PC2: ", percentVar[2], "% variance")) + 
  #geom_text_repel() +
  ggtitle("PCA") +
  guides(fill=guide_legend(ncol=4)) +
  theme(plot.title = element_text(hjust=0.5,size=10,face="bold"),
        axis.text = element_text(size=10,color="black"),
        axis.title = element_text(size=10,color="black"),
        axis.line = element_line(size=1),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size=8),
        panel.background = element_blank())
ggsave(filename="PCA_Males.png", height=2.5,width=2.5,units="in",dpi=300,
       path = "INSERT PATH")



#6C Venn Diagrams ----
res_group_LPSControl_vs_SalineControl_DEGs_UP <- res_group_LPSControl_vs_SalineControl_DEGs_UP %>% filter(symbols!="NA")
res_group_LPSKD_vs_SalineKD_DEGs_UP <- res_group_LPSKD_vs_SalineKD_DEGs_UP %>% filter(symbols!="NA")

x <- list(A=sample(res_group_LPSControl_vs_SalineControl_DEGs_UP$pvalue),
          B=sample(res_group_LPSKD_vs_SalineKD_DEGs_UP$Ensemble))
ggVennDiagram(x,
              category.names = c("Control","KD"),
              label = "count") +
  labs(title = "Upreguated
by LPS") +
  scale_fill_gradient(low="#0f99b2",high = "white") +
  scale_color_manual(values = c("#0f99b2","black")) +
  theme(plot.title = element_text(hjust=0.5,size=12,face="bold"),
        legend.position = "none")

ggsave(filename="VennDiagram_DEGs_UP_Males.png",height=2.5,width=2.5,units="in",dpi=300,
       path = "PATH")


res_group_LPSControl_vs_SalineControl_DEGs_DOWN <- res_group_LPSControl_vs_SalineControl_DEGs_DOWN %>% filter(symbols!="NA")
res_group_LPSKD_vs_SalineKD_DEGs_DOWN <- res_group_LPSKD_vs_SalineKD_DEGs_DOWN %>% filter(symbols!="NA")

y <- list(A=sample(res_group_LPSControl_vs_SalineControl_DEGs_DOWN$symbols),
          B=sample(res_group_LPSKD_vs_SalineKD_DEGs_DOWN$symbols))
ggVennDiagram(y,
              category.names = c("Control","KD"),
              label = "count") +
  labs(title = "Downregulated
by LPS") +
  scale_fill_gradient(low="#ece6ca",high = "white") +
  scale_color_manual(values = c("#ece6ca","black")) +
  theme(plot.title = element_text(hjust=0.5,size=12,face="bold"),
        legend.position = "none")

ggsave(filename="VennDiagram_DEGs_DOWN_Males.png",height=2.5,width=2.5,units="in",dpi=300,
       path = "PATH")

#6D Heatmap ----- 

#Clustering Coldata key
  #1 = UP in Control ONLY
  #2 = UP in Knockdown Only
  #3 = UP in BOTH Control and KD
  #4 = Down in Control ONLY
  #5 Down in LPS Only
  #6 Down in LPS and CONTROL

Clustering.Coldata <- read.csv("Clustering Coldata.csv", header = TRUE, sep = ",")
Clustering.Coldata <- data.frame(Clustering.Coldata,
                                 Cluster2 = ifelse(Clustering.Coldata$Cluster==1,"UP: Control",
                                                   ifelse(Clustering.Coldata$Cluster==2,"UP: KD",
                                                          ifelse(Clustering.Coldata$Cluster==3,"UP: Control & KD",
                                                                 ifelse(Clustering.Coldata$Cluster==4,"DOWN: Control",
                                                                        ifelse(Clustering.Coldata$Cluster==5,"DOWN: KD",
                                                                               ifelse(Clustering.Coldata$Cluster==6,"DOWN: Control & KD","")))))))

Clustering.Coldata$Cluster2 <- factor(Clustering.Coldata$Cluster2, levels = c("UP: Control & KD","UP: Control","UP: KD", "DOWN: Control & KD","DOWN: Control","DOWN: KD")) 
                                                                               
cluster.ensembles <- Clustering.Coldata$Ensemble

a <- data.frame(Countdata)
a <- rownames_to_column(a, var="Ensemble")
a <- rbind.data.frame(data.frame((a %>% filter(Ensemble %in% res_group_LPSControl_vs_SalineControl_DEGs_UP$Ensemble))),
                       data.frame((a %>% filter(Ensemble %in% res_group_LPSKD_vs_SalineKD_DEGs_UP$Ensemble))),
                       data.frame((a %>% filter(Ensemble %in% res_group_LPSControl_vs_SalineControl_DEGs_DOWN$Ensemble))),
                       data.frame((a %>% filter(Ensemble %in% res_group_LPSKD_vs_SalineKD_DEGs_DOWN$Ensemble))))
a <- data.frame(a[1],
                t(scale(t(a[2:21]))))
                #t(scale(t(a[11:21]))))

b <- a %>% filter(Ensemble %in% cluster.ensembles) 

b <- right_join(Clustering.Coldata,
                               a,
                               by="Ensemble")

b <- b %>% filter(Cluster!="NA")
b <- arrange(b, match(b$Cluster2, levels(b$Cluster2)))
b <- data.frame(b,
                placeholder= 1:5395)

countdata.heatmap <- melt(b,id.vars = c("Cluster","Cluster2","Ensemble","placeholder"))
countdata.heatmap$Cluster2 <- factor(countdata.heatmap$Cluster2)

C1.heatmap <- ggplot(countdata.heatmap, aes(x = rev(placeholder), y = variable, color = value)) +
  geom_tile() +
  coord_flip() +
  #scale_x_discrete(limits=rev) +
  labs(color="Z-score") +
  scale_colour_gradientn(colours=c('#075AFF','#FFFFCC','#FF0000'),
                         limits=c(-3,3)) +
  scale_y_discrete(labels=c('ControlSaline01','ControlSaline02','ControlSaline03','ControlSaline04',
                            'ControlLPS01','ControlLPS02','ControlLPS03','ControlLPS04','ControlLPS05',
                            'KDSaline01','KDSaline02','KDSaline03','KDSaline04','KDSaline05','KDSaline06','KDSaline07',
                            'KDLPS01','KDLPS02','KDLPS03','KDLPS04')) +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(angle=90,vjust=0.5),
        legend.position = "left",
        axis.ticks.y = element_blank(),
        panel.background = element_blank(),
        axis.title=element_blank())

C1.heatmap

Cluster <- ggplot(countdata.heatmap, aes(x = rev(placeholder), y = 1, fill = Cluster2)) +
  geom_tile() +
  coord_flip() +
  # #scale_fill_discrete(name = "Cluster Profile",
  #                     labels=c('UP: Control', 'UP: KD','UP: Control & KD',
  #                              'DOWN: Control','DOWN: KD','DOWN: Control & KD')) +
  theme(legend.position="none",
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.background = element_blank())
Cluster

cluster.legend <- ggplot(countdata.heatmap, aes(x = placeholder, y = 1, fill = Cluster2)) +
  geom_tile() +
  coord_flip() +
  lims(x = c(0,0), y = c(0,0))+
  theme_void()+
  scale_fill_discrete(name = "Cluster Profile",
                      labels=c("UP: Control & KD",
                               "UP: Control",
                               "UP: KD",
                               "DOWN: Control & KD",
                               "DOWN: Control","
                               DOWN: KD")) +
  theme(legend.position="top",
        #legend.margin = c(0.5,0.5),
        
        legend.key.size = unit(1, "cm"),
        legend.text = element_text(size =  12),
        legend.title = element_text(size = 14, face = "bold"))+
  guides(colour = guide_legend(override.aes = list(size=8)))

cluster.legend
ggsave(filename="Cluster.Legend.png",height=1,width=6.85,units="in",dpi=300,
       path="PATH")

final.heatmap <- ggarrange(C1.heatmap,Cluster,
                           nrow=1,
                           widths = c(1, 0.07),
                           #common.legend = T,
                           align="h")
final.heatmap
ggsave(filename="Heatmap_Males.png",height=6,width=5,units="in",dpi=300,
       path = "PATH")

#6E - Create Volcano Plots for SalineLPS vs. KDLPS comparison. ----

#Saline Control vs. LPS Control

keyvals.shape <- ifelse(
  rownames(res) %in% celltype1, 17,
  ifelse(rownames(res) %in% celltype2, 64,
         3))
keyvals.shape[is.na(keyvals.shape)] <- 3
names(keyvals.shape)[keyvals.shape == 3] <- 'PBMC'
names(keyvals.shape)[keyvals.shape == 17] <- 'Cell-type 1'
names(keyvals.shape)[keyvals.shape == 64] <- 'Cell-type 2'

#Modify Shape
keyvals.shape <- ifelse(
  res_group_LPSControl_vs_LPSKD$log2FoldChange < -1, 19,
  ifelse(res_group_LPSControl_vs_LPSKD$log2FoldChange > 1, 21,
         21))
keyvals.shape[is.na(keyvals.shape)] <- 21
names(keyvals.shape)[keyvals.shape == 21] <- 'high'
names(keyvals.shape)[keyvals.shape == 19] <- 'mid'
names(keyvals.shape)[keyvals.shape == 21] <- 'low'

#Modify Color
keyvals.colour <- ifelse(
  res_group_LPSControl_vs_LPSKD$log2FoldChange < -1, "#0f99b2",
  ifelse(res_group_LPSControl_vs_LPSKD$log2FoldChange > 1, "#0f99b2",
         ''))
keyvals.colour[is.na(keyvals.colour)] <- 'white'
  names(keyvals.colour)[keyvals.colour == "#0f99b2"] <- 'high'
  names(keyvals.colour)[keyvals.colour == 'white'] <- 'mid'
  names(keyvals.colour)[keyvals.colour == "#0f99b2"] <- 'low'
  
  #Switch FC values to make DOWN in SLC11a2KD on the right            
  res_group_LPSControl_vs_LPSKD <- res_group_LPSControl_vs_LPSKD %>% dplyr::mutate(LFC=log2FoldChange*-1) 
  
  Volcano <- EnhancedVolcano(res_group_LPSControl_vs_LPSKD,
                             lab = res_group_LPSControl_vs_LPSKD$symbols,
                             x='LFC',
                             y='padj',
                             selectLab = c('Slc40a1','Slc11a2','Il1b','Il6','Tnfa', 'Fth1','Nos2','Gpx4','Cx3cr1',
                                           'Trem2','Cx3cr1','Slc39a14','Il19','Arg1',
                                           'Eef2k','Slc46a3','Gpr155','Slamf1','Cxcl11','Il2ra',
                                           'Cdk19','Nfkb2','Hif1a','Il4ra','Arg1','Tnfrsf9','Cxcl10','Ccl5'),
                             xlab = bquote(~Log[2]~ 'fold change'),
                             ylab = bquote(~-Log[10]~adjusted~italic(P)),
                             
                             title='LPSControl vs. KDLPS',
                             pCutoff = 0.05, 
                             FCcutoff = 1,
                             labSize = 2,
                             shapeCustom = keyvals.shape,
                             colCustom = keyvals.colour,
                             subtitle = NULL,
                             caption = NULL,
                             drawConnectors = T,
                             typeConnectors = "closed",
                             maxoverlapsConnectors =Inf,
                             min.segment.length = -0.6,
                             widthConnectors = 0.5,
                             legendPosition = 'none',
                             gridlines.major = FALSE,
                             lengthConnectors = unit(0, "npc"),
                             boxedLabels = T,
                             shape = 21,
                             titleLabSize = 12,
                             axisLabSize = 12,
                             legendIconSize = 0.5)
  Volcano
  ggsave(filename="Volcano.png", height=4.3,width=5,units="in",dpi=300,
         path = "PATH")

  
#Figure 6F 
#Gene Ontology Variable String_Response to LPS in Control vs. KD
#Data for DEGs up and DOWN were input into the DAVID respository and figures generated by ggplot 
# https://david.ncifcrf.gov/tools.jsp

DAVID_UP <- read.delim("DAVID_UP_ ControlLPS vs. KDLPS.txt")
DAVID_DOWN <- read.delim("DAVID_DOWN_ ControlLPS vs. KDLPS.txt")
DAVID <- rbind.data.frame(
  data.frame(DAVID_UP %>% filter(FDR<0.05),
             Direction="Up"),
  data.frame(DAVID_DOWN[1:15,] %>% filter(FDR<0.05),
             Direction="Down"))
DAVID <- data.frame(DAVID,
                    FDR.NEW=ifelse(DAVID$Direction=="Down",(-log(DAVID$FDR)*-1),-log(DAVID$FDR)))
DAVID[c('GO', 'Term')] <- str_split_fixed(DAVID$Term, '~', 2)
DAVID <- data.frame(DAVID,
                    Fold.Enrichment2 = ifelse(DAVID$Direction=="Up",
                                              DAVID$Fold.Enrichment, DAVID$Fold.Enrichment*-1))



ggplot(DAVID, aes(x=FDR.NEW,y=reorder(Term,FDR.NEW),color=Fold.Enrichment2)) +
  geom_segment(aes(yend = reorder(Term,FDR.NEW), xend = 0),color ="black",size=0.7,alpha=0.5) +
  geom_point(aes(size=Count)) +
  labs(x="-log(FDR)",
       title="GO: Biological process
KDLPS vs. ControlLPS") +
  scale_colour_gradient2(low='#075AFF',
                         mid='#FFFFCC',
                         high='#FF0000',
                         midpoint = 0) +
  geom_vline(xintercept = 0,color="black",size=1) +
  theme(plot.title = element_text(hjust=0.5,size=10,face="bold"),
        axis.text = element_text(size=10,color="black"),
        axis.title = element_text(size=12,color="black",face='bold'),
        axis.title.y = element_blank(),
        axis.line = element_line(size=1),
        #panel.grid.major.y = element_line(colour = "grey60", linetype = "dashed"),
        legend.position = "none",
        panel.background = element_blank(),
        strip.text = element_text(face="bold",size=10)) +
  facet_grid(factor(Direction,levels=c("Up","Down")) ~ ., scales = "free_y", space = "free_y")
ggsave(filename="GO_KDLPS vs. ControlLPS.png", height=4,width=8,units="in",dpi=300,
       path = "PATH")

#Upregulated genes in control 
DAVID_UP_Control <- data.frame(read.delim("2_GO_DAVID_Up_Control.Only.txt"),Group="Control")
DAVID_UP_Control <- data.frame(DAVID_UP_Control)[c(1:10),]
DAVID_UP_Control[c('GO', 'Term')] <- str_split_fixed(DAVID_UP_Control$Term, '~', 2)


#DAVID_UP_Control 
DAVID_Up_Control_Plot <- ggplot(DAVID_UP_Control, aes(x=-log(FDR),y=reorder(Term,-log(FDR))),color=Count) +
  geom_segment(aes(yend = reorder(Term,-log(FDR)), xend = 0),color ="black",size=0.7,alpha=0.5) +
  geom_point(aes(size=Count,fill=Count),shape=21) +
  labs(x="-log(FDR)",
       title = "Upregulated by LPS") +
  scale_fill_gradientn(colours=c("white","#0f99b2"),
                        limits = c(0,100)) +
  scale_x_continuous(limits = c(0,50)) +
  scale_size_continuous(limits = c(0,100)) +
  geom_vline(xintercept = 0,color="black",size=1) +
  theme(plot.title = element_text(hjust=0.5,size=14,face="bold"),
        axis.text.y = element_text(size=10,color="black"),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.line.y = element_line(size=1),
        axis.line.x=element_line(color="white"),
        #panel.grid.major.y = element_line(colour = "grey60", linetype = "dashed"),
        #legend.position = "bottom",
        panel.background = element_blank()) +
  facet_grid(factor(Group,levels=c("Control","KD")) ~ ., scales = "free_y", space = "free_y")
DAVID_Up_Control_Plot

#Upregulated genes in CONTROL & KD 
DAVID_UP_Control.KD <- data.frame(read.delim("2_GO_DAVID_Up_Control+KD.txt"),Group="Control & KD")
DAVID_UP_Control.KD <- data.frame(DAVID_UP_Control.KD)[c(1:10),]
DAVID_UP_Control.KD[c('GO', 'Term')] <- str_split_fixed(DAVID_UP_Control.KD$Term, '~', 2)


#DAVID_UP_Control & KD
DAVID_UP_Control.KD_Plot <- ggplot(DAVID_UP_Control.KD, aes(x=-log(FDR),y=reorder(Term,-log(FDR))),color=Count) +
  geom_segment(aes(yend = reorder(Term,-log(FDR)), xend = 0),color ="black",size=0.7,alpha=0.5) +
  geom_point(aes(size=Count,fill=Count),shape=21) +
  labs(x="-log(FDR)") +
  scale_fill_gradientn(colours=c("white","#C3B294"),
                       limits = c(0,100)) +
  scale_x_continuous(limits = c(0,50)) +
  scale_size_continuous(limits = c(0,100)) +
  geom_vline(xintercept = 0,color="black",size=1) +
  theme(plot.title = element_text(hjust=0.5,size=10,face="bold"),
        axis.text.y = element_text(size=10,color="black"),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.line.y = element_line(size=1),
        axis.line.x=element_line(color="white"),
        #panel.grid.major.y = element_line(colour = "grey60", linetype = "dashed"),
        #legend.position = "bottom",
        panel.background = element_blank()) +
  facet_grid(factor(Group,levels=c("Control & KD")) ~ ., scales = "free_y", space = "free_y")
DAVID_UP_Control.KD_Plot

#Saline Control vs. LPS Control
keyvals <- ifelse(
  res_group_LPSControl_vs_SalineControl$log2FoldChange < -1, "#ece6ca",
  ifelse(res_group_LPSControl_vs_SalineControl$log2FoldChange > 1, "#0f99b2",
         'black'))
         keyvals[is.na(keyvals)] <- 'black'
           names(keyvals)[keyvals == "#0f99b2"] <- 'high'
           names(keyvals)[keyvals == 'black'] <- 'mid'
           names(keyvals)[keyvals == "#ece6ca"] <- 'low'
           
#END