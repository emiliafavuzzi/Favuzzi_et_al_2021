library(dplyr)
library(stringr)
library(SingleCellExperiment)
library(scran)
library(scater)
library(igraph)
library(ggplot2)
library(umap)
library(viridis)
library(Seurat)
#this R script starts with the output of the matlab script, for each brain there is a folder with subfolders in it organized by sliceXsideX (the original folders have been moved to a different directory but the organization is preserved)
#MERFISH count files are named as sliceX_sideXfastCountsFixedZCounts.csv
#smFISH count files are named as slice1_side2counted4genesNUCLEITHRESHHIGH.csv
setwd("/Volumes/GoogleDrive/My Drive/MERFISH_analysis")
#add info of brain, slices, side to the smFISH file and save a new file
files_to_process <- list.files("/Volumes/GoogleDrive/My Drive/MERFISH_analysis", pattern="HHIGH", recursive=T, full.names=F)
for(i in files_to_process)
{
  print(i)
  tmp <- read.csv(i, header = T, stringsAsFactors = F)
  name <- gsub("brain.*?/", "", i)
  name <- gsub("/.*", "", name)
  tmp$name <- name
  write.csv(tmp, file=paste0(i,'_smFISH.csv'))
}
#merge all files
files_to_process <- list.files("/Volumes/GoogleDrive/My Drive/MERFISH_analysis", pattern="Fixed", recursive=T, full.names=F)
out <- do.call(rbind,lapply(files_to_process,function(x){
  read.csv(x, header = T, stringsAsFactors = F)
}))
write.csv(file = 'MERFISH.csv',out)
files_to_process <- list.files("/Volumes/GoogleDrive/My Drive/MERFISH_analysis", pattern="smFISH", recursive=T, full.names=F)
out <- do.call(rbind,lapply(files_to_process,function(x){
  read.csv(x, header = T, stringsAsFactors = F)
}))
write.csv(file = 'smFISH.csv',out)
#merge MERFISH and smFISH files and add correct gene names/info
MERFISH <- read.csv("MERFISH.csv", row.names=1, sep=",", header = TRUE)
smFISH <- read.csv("smFISH.csv", row.names=1, sep=",", header = TRUE)
data <- cbind(MERFISH, smFISH[,-(1:5)]) #remove Tspan and Gabbr1 smFISH probes that we did not use
genenames <- c("Tmem119", "Gabbr1", "Gabbr2", "Cd164", "Clec4a3", "Ecscr", "Fcrls", "Gpr34", "Laptm4a", "Laptm5", "P2ry12", "P2ry13", "Ppib", "Selenok", "Selenop", "Tmem14c", "Trem2", "Tspan3", "Tspan4", "Tspan7", "cell.size", "zdepth", "C1qc", "Tmsb4x", "Rps29", "Ftl1", "mt_Nd2", "mt_Co3", "info")  
colnames(data) <- genenames
#add genotype label
data <- dplyr::mutate(data, genotype = ifelse(str_detect(data$info, "control.*"), "ctr", "ko"))
write.csv(file = 'ctrko_NN.csv',data)

#normalize by total counts/cell for each slice_side using a table derived from a separate count of all spots and DAPI cells
slice_factor <- read.table('scale_emc.tsv',header = T, row.names = 1, check.names = F, stringsAsFactors = F)
data[,'avg_scaling'] <- sapply(data$info,function(x){slice_factor[x,3]})
data1 <- data[,1:28] #selection of cols with counts
for(i in 1:ncol(data1)) {data1[ , i] <- data1[ , i] /(data$avg_scaling)*100} #apply slice scaling factor
data[,1:20] <- data1[,1:20] #we select the cols we want to replace so we do not change cellsize and zdepth
data[,23:28] <- data1[,23:28] #same as above
write.csv(file = 'merfish_batchcorr.csv',data) #save batch corrected file

# we will work with the batch corrected file but potentially the code below could also be applied to the raw data if you load ctrko_NN.csv file
data <- read.csv('merfish_batchcorr.csv', row.names=1, sep=",", header = TRUE)
#remove cells that have Tmem119+Fcrls_P2ry12 < 7 in any z plane (=0 z planes)
data <- dplyr::filter(data, zdepth >= 1)
#obtain RNA density (i.e. divide counts by cell volume)
volume <- data$cell.size*data$zdepth
p <- matrix(volume, nrow=nrow(data), ncol=ncol(data))
norm.data <- (data[, 1:28])/p*100000
data[, 1:28] <- norm.data
data <- data[,-c(21, 22, 29, 31)]#remove info, scaling, cell.size and z_depth cols as we do not need them anymore
data <- data %>% mutate(mt_Nd2 = mt_Nd2/100000, mt_Co3 = mt_Co3/100000) #scale mt genes because we had to use integrated density for these (note, we will remove them later as the signal was not clearly marking any cell groups)
data <- data %>% select(genotype, everything()) #move the genotype as first col
#remove cells that have half inside the mask 
data <- data[!(data$C1qc==0 & data$Tmsb4x== 0 & data$Rps29== 0 & data$Ftl1== 0),]
write.csv(file = 'merfish_norm.csv',data) #save normalized RNA density data that we will use for downstream analysis

#downstream analysis (clustering and DEGs)
data <- read.csv("merfish_norm.csv", sep=",", header = TRUE)
#import as single cell experiment
rownames(data) <- stringr::str_c('cell', seq(1:nrow(data)))
exp <- as.matrix(t(data[,3:28]))
cell.labels <- rownames(data)
genes <- rownames(exp)
sce <- SingleCellExperiment(assays = list(counts = exp),
                            colData=DataFrame(label=cell.labels),
                            rowData=DataFrame(label=genes),
                            metadata=list(study="MERFISH"))
colLabels(sce) <- data$genotype
genotype <- data$genotype
#normalize cell-specific biases (deconvolution strategy)
set.seed(100)
sce <- computeSumFactors(sce)
sce <- scater::logNormCounts(sce)
#export assay as we will need it for Seurat later
data <- assay(sce, "logcounts")
data2 <- as.data.frame(t(data))
data2 <- cbind(genotype, data2)
#cluster ctr+ko but without the 2 mt genes (not good)
work.data <- data2 %>% select(Tmem119, Gabbr1, Gabbr2, Cd164, Clec4a3, Ecscr, Fcrls, Gpr34, Laptm4a, Laptm5, P2ry12, P2ry13, Ppib, Selenok, Selenop, Tmem14c, Trem2, Tspan3, Tspan4, Tspan7, C1qc, Tmsb4x, Rps29, Ftl1)
work.data <- scale(work.data)
set.seed(123)
g <- buildSNNGraph(work.data, k=20, transposed=TRUE, type="jaccard")
#this above gives exactly same result as g <- buildSNNGraph(sce, k=20, type="jaccard") but is faster (and we will have to import it in Seurat anyway)
clust <- igraph::cluster_louvain(g)$membership
table(clust)
ctrko_cl <- cbind(data2, clust)
#PLOT CLUSTERS IN UMAP SPACE
set.seed(123)
ctrko.umap <- umap(work.data)
#plot UMAP with clusters (Fig 6C)
ctrko_cl2 <- ctrko_cl%>%dplyr::mutate_at("clust", as.character)
df <- data.frame(x = ctrko.umap$layout[,1],
                 y = ctrko.umap$layout[,2],
                 clust = ctrko_cl2[, 28])
ggplot(df, aes(x, y, colour = clust)) +
  geom_point(shape=20, size = 0.5)+ ylim(-5,5)+ #we set the limit and do not visualize cluster 2 (only 18 cells), remove the limit to see it
  scale_color_manual(values=c("#F67280", "#3D87F2", "#8F384D", "#99B898", "#547050", "#6C5B7B",  "#C06C84"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
##plot UMAP by genotype (Fig 6B)
ctrko_cl2 <- ctrko_cl%>%dplyr::mutate_at("genotype", as.character)
df4 <- data.frame(x = ctrko.umap$layout[,1],
                  y = ctrko.umap$layout[,2],
                  genotype = ctrko_cl2[, 1])
genotype = ctrko_cl2[, 1]
ggplot(df, aes(x, y, colour = genotype)) +
  geom_point(shape=20, size = 0.5)+ ylim(-5,5)+
  scale_color_manual(values=c("#a8a8a7","#F9A106"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
#plot Gabbr1 scaled expression (Fig 6D upper)
df2 <- data.frame(x = ctrko.umap$layout[,1],
                  y = ctrko.umap$layout[,2],
                  Gabbr1 = scale(ctrko_cl2[, 3]))
ggplot(df2, aes(x, y, colour = Gabbr1)) +
  geom_point(shape=20, size = 0.5)+ ylim(-5,5)+
  scale_color_viridis(limits = c(0, 1), oob = scales::squish)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
#plot Gabbr2 scaled expression (Fig 6D lower)
df2 <- data.frame(x = ctrko.umap$layout[,1],
                  y = ctrko.umap$layout[,2],
                  Gabbr2 = scale(ctrko_cl2[, 4]))
ggplot(df2, aes(x, y, colour = Gabbr2)) +
  geom_point(shape=20, size = 0.5)+ ylim(-5,5)+
  scale_color_viridis(limits = c(0, 1), oob = scales::squish)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))

#plot histogram of Gabbr1 and Gabbr2 to select positive cells and add a col indicating if positive or negative
hist(ctrko_cl$Gabbr1) #Gabbr1
d <- density(ctrko_cl$Gabbr1)
x <- ctrko_cl$Gabbr1
v <- optimize(approxfun(d$x,d$y),interval=c(0,6))$minimum #v is the local minimum (interval set btw peaks) = 0.59, to be more restrictive we use 1
abline(v=v, col="blue")
hist(ctrko_cl$Gabbr2) #Gabbr2
d <- density(ctrko_cl$Gabbr2)
x <- ctrko_cl$Gabbr2
v <- optimize(approxfun(d$x,d$y),interval=c(0,6))$minimum #v is the local minimum (interval set btw peaks) = 0.65, to be more restrictive we use 1
abline(v=v, col="blue")
#add positive or negative condition
ctrko_cl_cond <- dplyr::mutate(ctrko_cl, posneg = ifelse((Gabbr1>= 1 & Gabbr2>=1), "positive", "negative")) 
write.csv(file = 'ctr&ko_cl_cond.csv',ctrko_cl_cond)

#Figure 6E was generated by going back to matlab and using the csv generated below
formaps <- read.csv('ctrko_NN.csv', sep=",", header = TRUE)
formaps <- dplyr::filter(formaps, zdepth >= 1)
formaps <- formaps[!(formaps$C1qc==0 & formaps$Tmsb4x== 0 & formaps$Rps29== 0 & formaps$Ftl1== 0),]
posneg <- ctrko_cl_cond$posneg
info <- formaps$info %>% str_replace_all(c("control" = "ctr", "KO" = "ko"))
formaps <- cbind(formaps[,-30], info, clust, posneg)
formaps <- formaps %>% select(X, Tmem119, Gabbr1, Gabbr2, cell.size, info, clust, posneg)
write.csv(file = 'ctrko_formaps.csv',formaps)
#Figure 6F was generated using an ImageJ script applied on the maps

#Figure 6G-- the specific clusters were manually adjusted in the figure to facilitate direct comparisons
library(Useful2me)
ctrko_cl_cond  <- ctrko_cl_cond%>%dplyr::mutate_at("clust", as.character)
toplot <- subset(ctrko_cl_cond, genotype == "ctr"& (posneg == "positive" & clust %in% c('5', '6')) | (posneg == "negative" & clust %in% c('4', '7')))
ggplot(toplot, aes(clust, C1qc, fill=clust, color=clust)) +
  geom_split_violin(position=position_dodge(0.7))+
  scale_fill_manual(values=c("#99B898", "#547050", "#6C5B7B", "#C06C84"))+
  scale_color_manual(values = c("#99B898", "#547050", "#6C5B7B", "#C06C84"))+
  stat_summary(fun=mean, geom="point", shape=16, size=3, color="#a7a6a5")+ 
  theme_classic()
#Figure 6H (left) 
ggplot(toplot, aes(clust, Tspan3, fill=clust, color=clust)) +
  geom_split_violin(position=position_dodge(0.7))+
  scale_fill_manual(values=c("#99B898", "#547050", "#6C5B7B", "#C06C84"))+
  scale_color_manual(values = c("#99B898", "#547050", "#6C5B7B", "#C06C84"))+
  stat_summary(fun=mean, geom="point", shape=16, size=3, color="#a7a6a5")+ 
  theme_classic()
#Figure 6H (right)
ggplot(toplot, aes(clust, Tspan4, fill=clust, color=clust)) +
  geom_split_violin(position=position_dodge(0.7))+
  scale_fill_manual(values=c("#99B898", "#547050", "#6C5B7B", "#C06C84"))+
  scale_color_manual(values = c("#99B898", "#547050", "#6C5B7B", "#C06C84"))+
  stat_summary(fun=mean, geom="point", shape=16, size=3, color="#a7a6a5")+ 
  theme_classic()

#Find markers for each cluster (Table S2, last sheet)
data <- read.csv('ctr&ko_cl_cond.csv', row.names=1, header=T)
meta <- data[,c('genotype','posneg','clust')]
data <- data[,setdiff(colnames(data),c('genotype','posneg','clust', 'mt_Nd2', 'mt_Co3'))]
data <- t(data)
merfish <- CreateSeuratObject(data,meta=meta)
dir.create("cl_markers")
setwd("/Volumes/GoogleDrive/My Drive/MERFISH_analysis/cl_markers")
Idents(merfish) = 'clust'
cluster1.markers <- FindMarkers(merfish, ident.1 = '1')
cluster2.markers <- FindMarkers(merfish, ident.1 = '2')
cluster3.markers <- FindMarkers(merfish, ident.1 = '3')
cluster4.markers <- FindMarkers(merfish, ident.1 = '4')
cluster5.markers <- FindMarkers(merfish, ident.1 = '5')
cluster6.markers <- FindMarkers(merfish, ident.1 = '6')
cluster7.markers <- FindMarkers(merfish, ident.1 = '7')
write.csv(cluster1.markers, file='ctrkocluster1.markers.csv')
write.csv(cluster2.markers, file='ctrkocluster2.markers.csv')
write.csv(cluster3.markers, file='ctrkocluster3.markers.csv')
write.csv(cluster4.markers, file='ctrkocluster4.markers.csv')
write.csv(cluster5.markers, file='ctrkocluster5.markers.csv')
write.csv(cluster6.markers, file='ctrkocluster6.markers.csv')
write.csv(cluster7.markers, file='ctrkocluster7.markers.csv')
#in Figure 6, clusters were renamed such that an easier direct comparison with scRNA-seq clusters could be made
#cl 7 is 3, cl 6 is 3GG, cl 5 is 4GG, cl 4 is 4, cl 1 is 1, cl 2 is not visualized (small-weird, 18 cells), cl 3 is 2GG
#as described in methods, in the code below we only used positive cells from the GG cluster and negative cells from non-GG clusters, the result is similar without this selection (with, as expected, slightly lower fold changes)
#find DEGs ctr vs ko by cluster (Table S2 and input for Figure 6I)
for( cl in c( '3', '5', '6'))
{
  print(cl)
  cmp = subset(merfish, posneg == 'positive' & clust == cl )
  Idents(cmp) = 'genotype'
  mark = FindMarkers(cmp, ident.1='ko',ident.2='ctr')
  write.csv(mark, file=paste0(cl,'_cmpP.csv'))
}
for( cl in c('1', '4', '7'))
{
  print(cl)
  cmp = subset(merfish, posneg == 'negative' & clust == cl )
  Idents(cmp) = 'genotype'
  mark = FindMarkers(cmp, ident.1='ko',ident.2='ctr')
  write.csv(mark, file=paste0(cl,'_cmpN.csv'))
}

#MERFISH WT control (1 to 1 comparison in Table S2 and used for significance in Figure 6G and H)
#compare pruning genes in ctr cl 4 vs 7 (i.e. 4 vs 3)
DEGs4vs7_ctr <- FindMarkers(subset(merfish, genotype == 'ctr' & posneg == 'negative'), ident.1 = '4', ident.2='7')
write.csv(DEGs4vs7_ctr, file='DEGs4vs3_ctr.csv')
#compare pruning genes in ctr cl 5 vs 6 (i.e. 4GG vs 3GG)
DEGs4GGvs3GG_ctr <- FindMarkers(subset(merfish, genotype == 'ctr' & posneg == 'positive'), ident.1 = '5', ident.2='6')
write.csv(DEGs4GGvs3GG_ctr, file='DEGs4GGvs3GG_ctr.csv')
#compare 4GG vs 4 (Table S2 and significance in Figure 6G and 6H)
cmp = subset(merfish, genotype == 'ctr' & ((posneg == 'positive' & clust=='5')|(posneg == 'negative' & clust=='4')))
DEGs4GGvs4_ctr <- FindMarkers(cmp, ident.1 = '5', ident.2='4')
write.csv(DEGs4GGvs4_ctr, file='DEGs4GGPvs4N_ctr.csv')
#compare 3GG vs 3 (Table S2 and significance in Figure 6H)
cmp = subset(merfish, genotype == 'ctr' & ((posneg == 'positive' & clust=='6')|(posneg == 'negative' & clust=='7')))
DEGs3GGvs3_ctr <- FindMarkers(cmp, ident.1 = '6', ident.2='7')
write.csv(DEGs3GGvs3_ctr, file='DEGs3GGPvs3N_ctr.csv')
#MERFISH KO 1 to 1 comparison in Table S2
cmp <- subset(merfish, genotype == 'ko' & ((posneg == 'positive' & clust=='5')|(posneg == 'negative' & clust=='4')))
DEGs4GGvs4_ko <- FindMarkers(cmp, ident.1 = '5', ident.2='4')
write.csv(DEGs4GGvs4_ko, file='DEGs4GGPvs4N_ko.csv')
