################################################## 2. clustering && annotation 
setwd('~/Analysis/')
library(Seurat)
library(dplyr)
library(cowplot)

obj_merged <- readRDS('1-merged/merged.filtered_low.quality.rds')
obj_merged <- NormalizeData(obj_merged)
obj_merged <- FindVariableFeatures(obj_merged,nfeatures = 2000,selection.method = 'vst')
obj_merged <- ScaleData(obj_merged)
obj_merged <- RunPCA(obj_merged)
ElbowPlot(obj_merged,ndims = 50)
obj_merged <- FindNeighbors(obj_merged,dims = 1:15)
obj_merged <- FindClusters(obj_merged,resolution = 0.3)
obj_merged <- RunUMAP(obj_merged,dims = 1:15)

p1 <- DimPlot(obj_merged,label = TRUE)
p2 <- DimPlot(obj_merged,group.by = 'orig.ident')
p3 <- DimPlot(obj_merged,group.by = 'group')
plot_grid(p1,p2,p3,ncol = 2)

FeaturePlot(obj_merged, features = c('PTPRC','EPCAM'),order = TRUE)
immune.markers <- c('PTPRC',
                    'CD3D','CD4','CD8A', #  T
                    'NCAM1','NCR1', # NK
                    'CD14','C1QA', # Macro
                    'CD1C', # DC
                    'S100A12', 'CLEC10A', # Mono
                    'CPA3','TPSAB1', # Mast
                    'CD79A','MZB1' # B/Plasma
)
FeaturePlot(obj_merged, features = immune.markers, ncol = 5,order = TRUE,min.cutoff = 'q3')

major.markers <- c('PTPRC','CD3D','CD14','EPCAM','PECAM1','COL1A2')
FeaturePlot(obj_merged, features = major.markers, ncol = 3,order = TRUE,min.cutoff = 'q3')

##### filter epi & immune double postive cells 
tt2 <- subset(obj_merged, subset = (`EPCAM` >0 | `KRT18`>0) & (`PTPRC` > 0 | `CD3D` > 0) ,slot = 'counts')
DimPlot(obj_merged,cells.highlight = colnames(tt2),sizes.highlight = 0.001)
obj_merged.filtered <- subset(obj_merged, cells= colnames(tt2),invert=TRUE)
obj_merged.filtered <- NormalizeData(obj_merged.filtered)
obj_merged.filtered <- FindVariableFeatures(obj_merged.filtered,nfeatures = 1500,selection.method = 'vst')
obj_merged.filtered <- ScaleData(obj_merged.filtered, vars.to.regress = 'percent.mt')
obj_merged.filtered <- RunPCA(obj_merged.filtered)
ElbowPlot(obj_merged.filtered,ndims = 50)

obj_merged.filtered <- FindNeighbors(obj_merged.filtered,dims = 1:nPCs)
obj_merged.filtered <- FindClusters(obj_merged.filtered,resolution = 0.25)
obj_merged.filtered <- RunUMAP(obj_merged.filtered,dims = 1:nPCs)

p1 <- DimPlot(obj_merged.filtered,label = TRUE)
p2 <- DimPlot(obj_merged.filtered,group.by = 'orig.ident')
p3 <- DimPlot(obj_merged.filtered,group.by = 'group')
plot_grid(p1,p2,p3,ncol = 2)
FeaturePlot(obj_merged.filtered, features = major.markers, ncol = 3,order = TRUE,min.cutoff = 'q3')
FeaturePlot(obj_merged.filtered, features = c('nCount_RNA','nFeature_RNA','percent.mt'), order = TRUE,min.cutoff = 'q3')

major.markers <- c('PTPRC','CD3D','CD3E','CD4','CD8A',
                   'NCAM1','NCR1','CD14','C1QA','C1QC','CD1C','CSF3R','FCGR3B','S100A12','CLEC10A',
                   'CPA3','TPSAB1',
                   'CD79A','MZB1',
                   'EPCAM','KRT18','PECAM1','COL1A2','MYH11')
FeaturePlot(obj_merged.filtered, features = major.markers, ncol = 6,order = TRUE,min.cutoff = 'q3')
# all.markers.ONB.filtered <- FindAllMarkers()

##### fast DE
Rcpp::sourceCpp("fast_de.cpp") # a function for DE which is write by C++
all.markers.ONB.filtered <- t_test_1vA(obj_merged.filtered@assays$RNA@data, as.character(obj_merged.filtered$seurat_clusters))
write.csv(all.markers.ONB.filtered,'2-reduction_clustering_anno/all.markers.csv')
all.markers.ONB.filtered_2 <- all.markers.ONB.filtered %>% dplyr::filter(pct.1>0.3,avg_logFC>0) %>% group_by(cluster) %>% arrange(desc(avg_logFC),.by_group=TRUE)
top20.markers <- all.markers.ONB.filtered_2 %>% group_by(cluster) %>% top_n(20,avg_logFC)
top20.df <- c()
for(i in 0:(length(unique(top20.markers$cluster))-1)){
  tmp <- as.character(top20.markers[top20.markers$cluster==i,]$gene)
  top20.df <- cbind(top20.df,tmp)
}
colnames(top20.df) <- paste0("C_",0:(length(unique(top20.markers$cluster))-1))
write.csv(top20.df,'2-reduction_clustering_anno/Top20.markers.csv')
##### cell annotation（res=0.25）
Idents(obj_merged.filtered) <- obj_merged.filtered$RNA_snn_res.0.25
obj_merged.filtered <- RenameIdents(obj_merged.filtered,`0`='NK/T',`1`='Endo',`2`='Mono',`3`='Macro',
                                    `4`='Epi',`5`='Epi',`6`='SMCs',`7`='Neu',`8`='Epi',`9`='Epi',
                                    `10`='B',`11`='Epi',`12`='Fibro',`13`='Epi',`14`='Epi',`15`='Fibro',
                                    `16`='Epi',`17`='Mast',`18`='Plasma',`19`='Epi',`20`='Epi',`21`='Epi',`22`='Endo')
p4 <- DimPlot(obj_merged.filtered,label = TRUE)
p5 <- DimPlot(obj_merged.filtered,group.by = 'RNA_snn_res.0.25',label = TRUE)
# p3 <- DimPlot(obj_merged.filtered,group.by = 'group')
plot_grid(p4,p5,ncol = 2)

saveRDS(obj_merged.filtered,'2-reduction_clustering_anno/merged.anno_v1.rds')
