
#################################### 1. merged all samples 
setwd('~/Analysis/')
library(Seurat)
library(dplyr)
library(cowplot)
OM2 <- Read10X('OM-2/outs/filtered_feature_bc_matrix/')
ONB.199N <- Read10X('ONB-199N/outs/filtered_feature_bc_matrix/')
# get common genes
common_genes <- intersect(rownames(OM2),rownames(ONB.199N_mergedv))

OM2 <- OM2[common_genes,]
OM2 <- CreateSeuratObject(OM2,project = 'OM-2')
OM2[['percent.mt']] <- PercentageFeatureSet(OM2,pattern = '^MT-')
OM2 <- subset(OM2,subset= nFeature_RNA>200 & nFeature_RNA < 8000 & percent.mt < 10 & nCount_RNA > 400)
OM2$group <- 'Normal'
saveRDS(OM2,'~/Analysis/OM-2/OM-2.filtered.rds')

sample.list <- list.files('~/project/',pattern = 'ONB-') ## ONB sample cellranger output 
obj.list <- list()
for(s in sample.list){
  tmp <- Read10X(paste0('~/project/',s,'/outs/filtered_feature_bc_matrix/'))
  tmp <- merge_v(tmp)
  tmp <- tmp[common_genes,]
  tmp <- CreateSeuratObject(tmp,project = s)
  tmp$group <- 'ONB'
  tmp[['percent.mt']] <- PercentageFeatureSet(tmp,pattern = '^MT-')
  obj.list[[s]] <- tmp
}

obj_merged <- merge(x=obj.list[[1]],y=obj.list[2:length(obj.list)])
obj_merged <- merge(x=obj_merged,y=OM2)
obj_merged <- subset(obj_merged,subset= nFeature_RNA>200 & nFeature_RNA < 8000 & percent.mt < 10 & nCount_RNA > 400)

# remove genes which expressed less than 10 cells
obj_merged <- obj_merged[rownames(obj_merged@assays$RNA@data)[Matrix::rowSums(obj_merged@assays$RNA@data > 0) > 10],]
saveRDS(obj_merged,'~/Analysis/1-merged/merged.filtered_low.quality.rds')

RidgePlot(obj_merged,features = c('nCount_RNA','nFeature_RNA','percent.mt'),ncol = 2)

VlnPlot(obj_merged,features = c('nCount_RNA','nFeature_RNA','percent.mt'),ncol = 2,pt.size = 0)

## sample info stat
source('~/bin/my_scripts/R_func/median_statistics_2.R') # function for statistics key features
obj_merged <- readRDS('~/Analysis/1-merged/merged.filtered_low.quality.rds')
sample.info <- CalMedian(obj_merged,features = c('nCount_RNA','nFeature_RNA','percent.mt'))
write.csv(sample.info,'~/Analysis/1-merged/filtered.sample.info.csv')
