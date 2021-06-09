library(Seurat)
library(monocle3)
library(ggplot2)
library(plyr)
library(dplyr)
library(reshape2)
library(stringr)
library(cowplot)
library(ggforce)
library(dendextend)

sessionInfo()

# setwd('./HHGGNDRXX/')
flist = list.files('.')
sl = lapply(flist[grepl('h5',flist)],function(f){
CreateSeuratObject(Read10X_h5(f),project=stringr::str_replace(f,'.h5',''))
})

 names(sl) = sapply(sl,function(x)levels(x@meta.data$orig.ident))


for(d in grep('WT',names(sl))) {
	sl[[d]]@meta.data$cond = 'WT'
}

for(d in grep('KO',names(sl))) {
	sl[[d]]@meta.data$cond = 'KO'
}

for(d in grep('KO_2',names(sl))) {
	sl[[d]]@meta.data$sample = 'KO_2'
}

for(d in grep('KO_3',names(sl))) {
	sl[[d]]@meta.data$sample = 'KO_3'
}


for(d in grep('WT_1',names(sl))) {
	sl[[d]]@meta.data$sample = 'WT_1'
}

for(d in grep('WT_2',names(sl))) {
	sl[[d]]@meta.data$sample = 'WT_2'
}
for(d in grep('CTR',names(sl))) {
	sl[[d]]@meta.data$sample = 'CTR3'
}


big.data = merge(sl[[1]], y = sl[2:4], add.cell.ids = names(sl), project = "gabab")
big.data[['percent.mt']]=PercentageFeatureSet(big.data, pattern = "^mt-")

keepers = rownames(subset(big.data@meta.data, percent.mt <=10 & nCount_RNA<=20000 &  nFeature_RNA>=800))
big.data = subset(big.data,cells = keepers)

big.data@meta.data$out = 1
keepl=lapply(c('Cx3cr1','Fcrls'),function(g){which(big.data@assays$RNA@counts[g,]>0)})
keepers = purrr::reduce(keepl,intersect)
big.data@meta.data[keepers,'out'] = 0
big.data = subset(big.data,out==0)

keepers = which(big.data@assays$RNA@counts['Gabbr1',]>0)
big.data@meta.data[,'gab'] = '-'
big.data@meta.data[keepers,'gab'] = '+'

big.data@meta.data$gcl = paste(big.data@meta.data$cond,big.data@meta.data$gab,sep='')

big.data = RunPCA(ScaleData(FindVariableFeatures(NormalizeData(big.data))))


cds = new_cell_data_set(big.data@assays$RNA@counts,cell_metadata = big.data@meta.data,gene_metadata = data.frame(row.names=rownames(big.data@assays$RNA@counts),gene_short_name = rownames(big.data@assays$RNA@counts)))
cds@colData$sample = as.factor(cds@colData$sample)
cds = preprocess_cds(cds)
cds = align_cds(cds, alignment_group = 'sample',residual_model_formula_str='~cond',verbose=T )
cds = reduce_dimension(cds)
cds = cluster_cells(cds)
cds = learn_graph(cds)

orig_cds = readRDS('./monocle_emc.RDS')

# original analysis was done with monocle3 cds version 0.2.0
# current cds version 0.2.3.0 has some small differences 
# older preprocessed object is provided for reproducibility 

orig_cds@metadata$cds_version
cds@metadata$cds_version
cowplot::plot_grid(
plot_cells(orig_cds,show_trajectory_graph=F,label_leaves=F,label_branch_points=F,label_cell_groups = F),
plot_cells(cds,show_trajectory_graph=F,label_leaves=F,label_branch_points=F,label_cell_groups = F))

table(clusters(orig_cds),clusters(cds))

# original analysis was done with monocle3 cds version 0.2.0
# current cds version 0.2.3.0 has some small differences 
# older preprocessed object is provided for reproducibility 

wt_sub = subset(big.data, cond=='WT') 

cdsn = new_cell_data_set(wt_sub@assays$RNA@counts,cell_metadata = wt_sub@meta.data,
gene_metadata = data.frame(row.names=rownames(wt_sub@assays$RNA@counts),
gene_short_name = rownames(wt_sub@assays$RNA@counts)))
cdsn =   preprocess_cds(cdsn, num_dim = 50)
cdsn =        align_cds(cdsn, alignment_group = 'sample')
cdsn = reduce_dimension(cdsn)
cdsn =    cluster_cells(cdsn)
cdsn =      learn_graph(cdsn)
orig_wt_sub = readRDS('HHGGNDRXX/monocle_WT.RDS')

table(clusters(orig_wt_sub),clusters(cdsn))
cowplot::plot_grid(plot_cells(cdsn,show_trajectory_graph=F,label_leaves=F,label_branch_points=F),
plot_cells(orig_wt_sub,show_trajectory_graph=F,label_leaves=F,label_branch_points=F))

#Fig5A
plot_cells(orig_wt_sub,show_trajectory_graph=F,label_leaves=F,label_branch_points=F)
plot_cells(orig_wt_sub,show_trajectory_graph=F,label_leaves=F,label_branch_points=F,genes = c('Tmsb4x','P2ry12','C1qc'),norm_method = 'log',scale_to_range = F)

#Fig5B&C
plot_cells(orig_cds,show_trajectory_graph=F,label_leaves=F,label_branch_points=F,color_cells_by = 'cond') +scale_color_manual(values=c('#FDE1B4','#ECECEC'))
plot_cells(orig_cds,show_trajectory_graph=F,label_leaves=F,label_branch_points=F) +scale_color_manual(values=c('#C06C84', '#F67280', '#6C5B7B', '#99B898', '#9F7E69', '#3D87F2', '#B557F2', '#F4D345'))

#Fig5D
integrated_comp = reshape2::melt(100*prop.table(table(clusters(orig_cds,reduction_method = 'UMAP'),orig_cds$cond),margin=1))
colnames(integrated_comp) = c('cluster','condition','cell')
ggplot(integrated_comp, aes(x=cluster, y=cell,fill=condition)) + geom_col() + scale_x_continuous(breaks = 1:8) + 
geom_text(data=subset(integrated_comp,condition=='KO'),aes(y=abs(1.05-cell),label=paste(round(cell,1),'%'),x=cluster)) + 
ylab('Cell proportion') +scale_fill_manual(values=c('#FDE1B4','#ECECEC')) + theme_minimal()


#Fig5G

make_v = function(cl4, gene='Gpr34'){
cko4 =data.frame(t(normalized_counts(cl4[gene,cl4$cond=='KO'],'log',1)))
cko4$cond='cKO'

wt4 =data.frame(t(normalized_counts(cl4[gene,cl4$cond=='WT'],'log',1)))
wt4$cond='WT'

vpl = rbind(wt4,cko4)
vpl$cond = factor(vpl$cond,levels = c("WT", "cKO"))

ggplot(vpl,aes_string(x='cond',y=gene,fill='cond'))+geom_violin() + theme_minimal() + 
stat_summary(fun=median, geom="point")  +scale_fill_manual(values=c('#ECECEC','#FDE1B4'))
# ggsave(file=paste0(gene,"_lognormed_5G",".eps"))
}
cl4 = orig_cds[,clusters(orig_cds) == 4]
make_v(cl4,'C1qb')
make_v(cl4,'Trem2')
make_v(cl4,'P2ry12')
make_v(cl4,'Tmem176b')
make_v(cl4,'Cd53')
make_v(cl4,'Ecscr')
make_v(cl4,'Selenop')
make_v(cl4,'Tspan4')
make_v(cl4,'Tspan3')



#Fig 5H

c4d = c('Atp1b3', 'Tmem176b', 'Trem2', 'Serpinf1', 'Gpr34', 'Cd53', 'P2ry12', 'C1qb', 'Tmem176a', 'Fcgr3', 'Tspan7', 'Selenop', 'Tmem37', 'Serinc3', 'Rnase4', 'C1qc', 'Cd164', 'Tspan4', 'Cd9', 'Ppib', 'Creg1', 'Cd302', 'Tmem14c', 'Selenok', 'mt-Nd2', 'Tmed10', 'Ecscr', 'Serinc1', 'Laptm4a', 'Vsir', 'Clic1', 'Itm2c', 'mt-Nd1', 'Tram1', 'Spcs2', 'Tspan3')
c4u =c('Rps8', 'Gdi2', 'Ssh2', 'Cfl1', 'H3f3b', 'Clta', 'Malat1', 'Mef2c', 'Srgap2', 'Rplp0', 'Cyth4', 'Chd9', 'Bin1', 'Grap', 'Calm2', 'Mef2a', 'Ccnd1', 'Frmd4a', 'Rhoh', 'Nmt1', 'Hnrnpa2b1', 'Irf8', 'Senp2', 'Bhlhe41', 'Fchsd2', 'Fli1', 'Brd9', 'Arhgap5', 'Rps26', 'Fus', 'Tsc22d4', 'Clk1', 'Maf', 'Zcchc6', 'Cacna1d', 'Hk2', 'Lacc1', 'Maml3', 'AC149090.1', 'Akap13', 'Prpf4b', 'Anxa3')
tsm_ids = c4d
wt_agg_cds = cds[tsm_ids,clusters(orig_cds)%in% c(4)]

exprs_mat = as.matrix(scale(t(normalized_counts(wt_agg_cds)[tsm_ids, ])))
my_levels = c('WT','KO')



cl_ord = unlist(lapply(my_levels,function(x){
tmp = order.hclust(exprs_mat[wt_agg_cds$cond==x,] %>% dist %>% hclust)
rownames(exprs_mat[wt_agg_cds$cond==x,])[tmp]
}))


g_ord = do.call(rbind,lapply(my_levels,function(x){
colMeans(exprs_mat[wt_agg_cds$cond==x,])
}))
rownames(g_ord) = my_levels

max_g = sapply(colnames(g_ord),function(x){ my_levels[which.max(g_ord[,x])]} )

g_ord = unlist(lapply(my_levels,function(x){
	if(length(max_g[which(max_g==x)])==1){
		return(names(max_g[which(max_g==x)]))
	}
	if(length(max_g[which(max_g==x)])==0){ return('')}
	names(sort(colMeans(exprs_mat[wt_agg_cds$cond==x,names(max_g[which(max_g==x)])]),decreasing=T))
}))
g_ord = g_ord[which(g_ord!='')]
exprs_mat <- reshape2::melt(exprs_mat)
colnames(exprs_mat) <- c('Cell', 'Gene', 'Zscore')
exprs_mat$Gene <- as.character(exprs_mat$Gene)
scale_min =-3
scale_max =3
exprs_mat[exprs_mat['Zscore']< scale_min,'Zscore'] = scale_min
exprs_mat[exprs_mat['Zscore']> scale_max,'Zscore'] = scale_max
exprs_mat$Cell = factor(exprs_mat$Cell,levels=cl_ord)
exprs_mat$Gene = factor(exprs_mat$Gene,levels=rev(unique(g_ord)))
exprs_mat$Gene = factor(exprs_mat$Gene,levels=rev(tsm_ids))
exprs_mat$Group ='KO'

exprs_mat[rownames(subset(exprs_mat, Cell %in% colnames(wt_agg_cds[,wt_agg_cds$cond=='WT']))),'Group'] = 'WT'

exprs_mat$Group = factor(exprs_mat$Group,levels=my_levels)
ggplot(exprs_mat,aes_string(x='Cell',y='Gene',fill='Zscore')) + geom_tile() + scale_fill_viridis_c() + 
facet_grid(~Group,drop=T,scales='free_x',space='free_x') + 
theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())


cl_4 = orig_cds[,clusters(orig_cds)==4]
cl_4 = reduce_dimension(cl_4)
cl_4 = cluster_cells(cl_4)
cl_4 = learn_graph(cl_4)

#Fig 5i
plot_cells(cl_4, ,show_trajectory_graph=F,label_leaves=F,label_branch_points=F,cell_size = 2) + scale_color_manual(values=c('#6E936A','#FF8C5A'))
tpl =  reshape2::melt(prop.table(table(clusters(cl_4), cl_4$cond),margin = 1))
colnames(tpl) = c('sub_cluster','condition','cells')
tpl$sub_cluster = factor(tpl$sub_cluster) 
ggplot(tpl,aes(x=sub_cluster, y=cells, fill=condition)) + geom_col() + scale_fill_manual(values=c('#FDE1B4','#ECECEC'))  + theme_minimal() 




#S6A
# original object 
S6A = readRDS('e14to.P30.RDS')
plot_cells(S6A, color_cells_by = 'tp', show_trajectory_graph=F,label_leaves=F,label_branch_points=F, label_cell_groups = F)

# load('emc_apr_18.RDA')

#S6A rebuilt, results not exactly similar 

# E14_F_B10 = read.csv('ref_data/GSM3442006_E14_F_B10.dge.txt.gz',sep='\t',row.names=1,stringsAsFactors=F)
# E14_M_B11 = read.csv('ref_data/GSM3442007_E14_M_B11.dge.txt.gz',sep='\t',row.names=1,stringsAsFactors=F)
# E14_F_B12 = read.csv('ref_data/GSM3442008_E14_F_B12.dge.txt.gz',sep='\t',row.names=1,stringsAsFactors=F)
# E14_F_C1=  read.csv('ref_data/GSM3442009_E14_F_C1.dge.txt.gz',sep='\t',row.names=1,stringsAsFactors=F)
# E14_M_B9=  read.csv('ref_data/GSM3442010_E14_M_B9.dge.txt.gz',sep='\t',row.names=1,stringsAsFactors=F)
# E14_F_B6=  read.csv('ref_data/GSM3442011_E14_F_B6.dge.txt.gz',sep='\t',row.names=1,stringsAsFactors=F)
# E14_M_B7=  read.csv('ref_data/GSM3442012_E14_M_B7.dge.txt.gz',sep='\t',row.names=1,stringsAsFactors=F)
# E14_M_B8=  read.csv('ref_data/GSM3442013_E14_M_B8.dge.txt.gz',sep='\t',row.names=1,stringsAsFactors=F)
# P4_M_A4    = read.csv('ref_data/GSM3442016_P4_M_A4.dge.txt.gz',sep='\t',row.names=1,stringsAsFactors=F)
# P4_M_A5    = read.csv('ref_data/GSM3442017_P4_M_A5.dge.txt.gz',sep='\t',row.names=1,stringsAsFactors=F)
# P4_F_A6    = read.csv('ref_data/GSM3442018_P4_F_A6.dge.txt.gz',sep='\t',row.names=1,stringsAsFactors=F)
# P4_F_B3    = read.csv('ref_data/GSM3442019_P4_F_B3.dge.txt.gz',sep='\t',row.names=1,stringsAsFactors=F)
# P4_F_B4    = read.csv('ref_data/GSM3442020_P4_F_B4.dge.txt.gz',sep='\t',row.names=1,stringsAsFactors=F)
# P4_M_B5    = read.csv('ref_data/GSM3442021_P4_M_B5.dge.txt.gz',sep='\t',row.names=1,stringsAsFactors=F)
# P5_M_A1    = read.csv('ref_data/GSM3442014_P5_M_A1.dge.txt.gz',sep='\t',row.names=1,stringsAsFactors=F)
# P5_F_A2    = read.csv('ref_data/GSM3442015_P5_F_A2.dge.txt.gz',sep='\t',row.names=1,stringsAsFactors=F)
# P5_F_NOP_1 = read.csv('ref_data/GSM3442047_P5_female_nopercoll_1.dge.txt.gz',sep='\t',row.names=1,stringsAsFactors=F)
# P5_F_NOP_3 = read.csv('ref_data/GSM3442049_P5_female_nopercoll_3.dge.txt.gz',sep='\t',row.names=1,stringsAsFactors=F)
# P5_F_NOP_2 = read.csv('ref_data/GSM3442048_P5_female_nopercoll_2.dge.txt.gz',sep='\t',row.names=1,stringsAsFactors=F)
# P5_F_P_1   = read.csv('ref_data/GSM3442050_P5_female_percoll_1.dge.txt.gz',sep='\t',row.names=1,stringsAsFactors=F)
# P5_F_P_2   = read.csv('ref_data/GSM3442051_P5_female_percoll_2.dge.txt.gz',sep='\t',row.names=1,stringsAsFactors=F)
# P5_F_P_3   = read.csv('ref_data/GSM3442052_P5_female_percoll_3.dge.txt.gz',sep='\t',row.names=1,stringsAsFactors=F)
# P30_M_2    = read.csv('ref_data/GSM3442023_P30_Male_2.dge.txt.gz',sep='\t',row.names=1,stringsAsFactors=F)
# P30_M_1    = read.csv('ref_data/GSM3442022_P30_Male_1.dge.txt.gz',sep='\t',row.names=1,stringsAsFactors=F)
# P30_M_3    = read.csv('ref_data/GSM3442024_P30_male_3.dge.txt.gz',sep='\t',row.names=1,stringsAsFactors=F)
# P30_M_4    = read.csv('ref_data/GSM3442025_P30_male_4.dge.txt.gz',sep='\t',row.names=1,stringsAsFactors=F)

# ref_data = list('E14_F_B10' = E14_F_B10, 'E14_M_B11' = E14_M_B11, 'E14_F_B12' = E14_F_B12, 'E14_F_C1' = E14_F_C1,
#                 'E14_M_B9' = E14_M_B9, 'E14_F_B6' = E14_F_B6, 'E14_M_B7' = E14_M_B7, 'E14_M_B8' = E14_M_B8,
#                 'P4_M_A4'= P4_M_A4, 'P4_M_A5'= P4_M_A5, 'P4_F_A6'= P4_F_A6, 'P4_F_B3' = P4_F_B3,
#                 'P4_F_B4'= P4_F_B4, 'P4_M_B5'= P4_M_B5 , 'P5_M_A1'= P5_M_A1, 'P5_F_A2'= P5_F_A2,
#                 'P5_F_NOP_1' = P5_F_NOP_1, 'P5_F_NOP_3' = P5_F_NOP_3, 'P5_F_NOP_2' = P5_F_NOP_2,
#                 'P5_F_P_1'= P5_F_P_1,'P5_F_P_2' = P5_F_P_2  , 'P5_F_P_3'= P5_F_P_3,
#                 'P30_M_2'= P30_M_2, 'P30_M_1'= P30_M_1, 'P30_M_3'= P30_M_3 , 'P30_M_4'= P30_M_4)
# ref_data  = lapply(names(ref_data),function(x)CreateSeuratObject(ref_data[[x]]))
# names(ref_data) = c('E14_F_B10', 'E14_M_B11', 'E14_F_B12', 'E14_F_C1', 'E14_M_B9', 'E14_F_B6', 'E14_M_B7', 'E14_M_B8',
#                     'P4_M_A4', 'P4_M_A5', 'P4_F_A6', 'P4_F_B3', 'P4_F_B4', 'P4_M_B5',
#                     'P5_M_A1', 'P5_F_A2', 'P5_F_NOP_1', 'P5_F_NOP_3', 'P5_F_NOP_2', 'P5_F_P_1', 'P5_F_P_2', 'P5_F_P_3',
#                     'P30_M_2', 'P30_M_1', 'P30_M_3', 'P30_M_4')
                   
# ref_data[['P15_1']] = subset(big.data,sample == 'WT_1')
# ref_data[['P15_2']] = subset(big.data,sample == 'WT_2')
# for(x in names(ref_data)){
#     ref_data[[x]]@meta.data$orig.ident = x 
#     ref_data[[x]]@meta.data$batch = x 
#     ref_data[[x]]@meta.data$tp =  stringr::str_split(x,'_',simplify=T)[,1]
# }
# ref_data_cds = merge(ref_data[[1]],ref_data[2:28])
# rm(ref_data);gc();                   
# ref_data_cds = new_cell_data_set(ref_data@assays$RNA@counts,cell_metadata = ref_data@meta.data,
#                                  gene_metadata = data.frame(row.names=rownames(ref_data@assays$RNA@counts),
#                                                             gene_short_name = rownames(ref_data@assays$RNA@counts)))

# ref_data_cds$batch = dplyr::recode(ref_data_cds$orig.ident,'E14_F_B10' = 'E14_F', 'E14_F_B12' = 'E14_F', 'E14_F_B6' = 'E14_F',
#                                    'E14_F_C1' = 'E14_F', 'E14_M_B11' = 'E14_M', 'E14_M_B7' = 'E14_M',
#                                    'E14_M_B8' = 'E14_M', 'E14_M_B9' = 'E14_M','P30_M_1' = 'P30_M','P30_M_2' = 'P30_M',
#                                    'P30_M_3' = 'P30_M', 'P30_M_4' = 'P30_M', 'P4_F_A6' = 'P4_F', 'P4_F_B3' = 'P4_F',
#                                    'P4_F_B4' = 'P4_F', 'P4_M_A4' = 'P4_M', 'P4_M_A5' = 'P4_M', 'P4_M_B5' = 'P4_M',
#                                    'P5_F_A2' = 'P5_F', 'P5_F_NOP_1' = 'P5_F', 'P5_F_NOP_2' = 'P5_F',
#                                    'P5_F_NOP_3' = 'P5_F', 'P5_F_P_1' = 'P5_F_P', 'P5_F_P_2' = 'P5_F_P',
#                                    'P5_F_P_3' = 'P5_F_P', 'P5_M_A1' = 'P5_M', 'P15_1' = 'P15_WT', 'P15_2' = 'P15_WT')
# ref_data_cds$batch = factor(ref_data_cds$batch,levels=c('E14_F','E14_M','P4_F','P4_M',
#                                                         'P5_F', 'P5_F_P', 'P5_M',
#                                                         'P15_WT', 'P30_M'))
# ref_data_cds$batch = droplevels(ref_data_cds$batch)

# ref_data_cds =   preprocess_cds(ref_data_cds, num_dim = 100)
# ref_data_cds =        align_cds(ref_data_cds, alignment_group = 'orig.ident',
#                                 residual_model_formula_str = '~ batch')
# ref_data_cds = reduce_dimension(ref_data_cds,umap.fast_sgd=TRUE, cores=36)
# ref_data_cds =    cluster_cells(ref_data_cds)
# ref_data_cds =      learn_graph(ref_data_cds)
                   

#FIG S6D original 
cdsko = readRDS('FIG_S6D.RDS')
plot_cells(cdsko,label_leaves=F,show_trajectory_graph=F,label_roots=F);

#FIG S6D rebuilt 

# ko_sub = subset(big.data, cond=='KO') 
# cdsko = new_cell_data_set(ko_sub@assays$RNA@counts,cell_metadata = ko_sub@meta.data,gene_metadata = data.frame(row.names=rownames(ko_sub@assays$RNA@counts), gene_short_name = rownames(ko_sub@assays$RNA@counts)))
# cdsko =   preprocess_cds(cdsko, num_dim = 50)
# cdsko =        align_cds(cdsko, alignment_group = 'sample')
# cdsko = reduce_dimension(cdsko)
# cdsko =    cluster_cells(cdsko)
# cdsko =      learn_graph(cdsko)

# plot_cells(cdsko,label_leaves=F,show_trajectory_graph=F,label_roots=F);

# FIG S6 E

cdsn = orig_wt_sub
cds = orig_cds

cds$cluster_a=''
cds@colData[colnames(cdsn),'cluster_a'] = paste0('WT_',clusters(cdsn,reduction_method = 'UMAP')[intersect(colnames(cds),colnames(cdsn))])
cds@colData[colnames(cdsko),'cluster_a'] = paste0('KO_',clusters(cdsko,reduction_method = 'UMAP')[intersect(colnames(cds),colnames(cdsko))])

destination = reshape2::melt(table(cds$cluster_a,clusters(cds,reduction_method = 'UMAP')))
origin= reshape2::melt(table(cds$cluster_a,clusters(cds,reduction_method = 'UMAP')))
prop_destination = reshape2::melt(round(100*prop.table(table(cds$cluster_a,clusters(cds,reduction_method = 'UMAP')),2),2))
prop_origin = reshape2::melt(round(100*prop.table(table(cds$cluster_a,clusters(cds,reduction_method = 'UMAP')),1),2))
prop_origin = subset(prop_origin,value>0)


colnames(origin) = c('condition_cl','integrated_cl','cells')
colnames(destination) = c('condition_cl','integrated_cl','cells')
myrs = aggregate(cells~integrated_cl,data=destination,sum)
myrs$r = 0
myrs$r = c(0,0,.1,.2,.2,.3,.4,.5)
destination$r = sapply(destination$integrated_cl,function(x){myrs[x,'r']})

ggplot(destination) + geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = r , r = 1, amount = cells, fill = as.factor(condition_cl)), stat = 'pie') +
xlab('')+ylab('') + coord_fixed() + scale_fill_manual(values=c('#F9EEDE', '#F7E5BF', '#F4D892', '#F4CF9D', '#F9C470', '#F9A106', '#DD8505','#EEEEEE', '#CCCCCA', '#9B9B9B', '#757574', '#444444')) +theme_void()+theme(plot.title = element_text(hjust = 0.5,size = 7,vjust = 1.1)) + NoLegend()  + facet_wrap(~integrated_cl,nrow=2)
