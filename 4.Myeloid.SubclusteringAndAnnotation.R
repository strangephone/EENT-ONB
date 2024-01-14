###### Myeloid cell sub-clustering
dir.create('~/Analysis/sub_myeloid')
obj.ONB <- readRDS('2-reduction_clustering_anno/merged.annotation.rds')

################## harmony
library(harmony)
myeloid_subset.harmony <- myeloid_subset
myeloid_subset.harmony <- FindVariableFeatures(myeloid_subset.harmony, nfeatures = '2500',selection.method = 'vst')
myeloid_subset.harmony <- ScaleData(myeloid_subset.harmony,vars.to.regress = c('nFeature_RNA','nCount_RNA','percent.mt'))
myeloid_subset.harmony <- RunPCA(myeloid_subset.harmony)
myeloid_subset.harmony<- RunHarmony(myeloid_subset.harmony, group.by.vars = 'orig.ident',plot_convergence = TRUE,
                                    max.iter.harmony = 20,epsilon.cluster=-Inf,epsilon.harmony = -Inf)

myeloid_subset.harmony <- RunUMAP(myeloid_subset.harmony,dims = 1:12,reduction = 'harmony')
myeloid_subset.harmony <- FindNeighbors(myeloid_subset.harmony,dims = 1:12,reduction = 'harmony')
myeloid_subset.harmony <- FindClusters(myeloid_subset.harmony) # default 0.8
# myeloid_subset.harmony <- FindClusters(myeloid_subset.harmony, resolution = 0.6)
# myeloid_subset.harmony <- FindClusters(myeloid_subset.harmony,resolution = 1.2)

p1 <- DimPlot(myeloid_subset.harmony,label = TRUE,pt.size = 0.6)
p2 <- DimPlot(myeloid_subset.harmony,group.by = 'orig.ident',pt.size = 0.3)
p3 <- DimPlot(myeloid_subset.harmony,group.by = 'group',pt.size = 0.6)
plot_grid(p1,p2,p3,ncol = 2)

FeaturePlot(myeloid_subset.harmony,features = c('CD14','CD68','Mono_Macro6','CD1C','CPA3','CSF3R',
                                                'S100A8','S100A9','LAMP3'),order = TRUE)
markers_1 <- FindAllMarkers(myeloid_subset.harmony,only.pos = TRUE)
top20_sub_myeloid.harmony <- markers_1 %>% group_by(cluster) %>% top_n(20,wt=avg_logFC)
angrycell::mergeplot1(myeloid_subset.harmony,marker_to_plot = c('CD3D','CD14'))

########################## remove mix cells
myeloid_subset.harmony$anno_mix <- ifelse(myeloid_subset.harmony$RNA_snn_res.0.8 %in% c(6,12,14,17),
                                          'Mix','singlets')
DimPlot(myeloid_subset.harmony,group.by = 'anno_mix')
myeloid_subset.harmony <- subset(myeloid_subset.harmony, subset = anno_mix =='singlets')
DimPlot(myeloid_subset.harmony)

######################## dimentionality reduction
myeloid_subset.harmony <- FindVariableFeatures(myeloid_subset.harmony, nfeatures = '2000',selection.method = 'vst')
myeloid_subset.harmony <- ScaleData(myeloid_subset.harmony)
myeloid_subset.harmony <- RunPCA(myeloid_subset.harmony)
myeloid_subset.harmony<- RunHarmony(myeloid_subset.harmony, group.by.vars = 'orig.ident',plot_convergence = TRUE,
                                    max.iter.harmony = 20,epsilon.cluster=-Inf,epsilon.harmony = -Inf)


myeloid_subset.harmony <- RunUMAP(myeloid_subset.harmony,dims = 1:10,reduction = 'harmony')
myeloid_subset.harmony <- FindNeighbors(myeloid_subset.harmony,dims = 1:10,reduction = 'harmony')
myeloid_subset.harmony <- FindClusters(myeloid_subset.harmony)

p1 <- DimPlot(myeloid_subset.harmony,label = TRUE,pt.size = 0.6)
p2 <- DimPlot(myeloid_subset.harmony,group.by = 'orig.ident',pt.size = 0.3)
p3 <- DimPlot(myeloid_subset.harmony,group.by = 'group',pt.size = 0.6)
plot_grid(p1,p2,p3,ncol = 2)
saveRDS(myeloid_subset.harmony, '~/Analysis/sub_myeloid/20211229-myeloid.subtype.v1.rds')

FeaturePlot(myeloid_subset.harmony,features = c('CD14','CD68','Mono_Macro6','CD1C','CPA3','CSF3R',
                                                'S100A8','S100A9','LAMP3'),order = TRUE)
all.markers_sub_myeloid.harmony <- FindAllMarkers(myeloid_subset.harmony,only.pos = TRUE)
top20_sub_myeloid.harmony <- all.markers_sub_myeloid.harmony %>% group_by(cluster) %>% top_n(20,wt=avg_logFC)

################### annotation major type
Idents(myeloid_subset.harmony) <- myeloid_subset.harmony$RNA_snn_res.0.8
myeloid_subset.harmony <- RenameIdents(myeloid_subset.harmony, `0`='Neu',`1`='Macro',`2`='Macro',
                                       `3`='Mono',`4`='Mono',`5`='Macro',`6`='Mast',`7`='Macro',
                                       `8`='DC',`9`='Mono',`10`='Macro',`11`='Neu',`12`='IFN-response',`13`='DC')
myeloid_subset.harmony$Major_2 <- Idents(myeloid_subset.harmony)
DimPlot(myeloid_subset.harmony,label = T,group.by = 'Major_2')
saveRDS(myeloid_subset.harmony,'~/Analysis/sub_myeloid/myeloid.rds')
################### annotation subtype

FeaturePlot(myeloid_subset.harmony,features = c('CD14','FCGR3A'))
# mono
DotPlot(myeloid_subset.harmony,features = c('FCN1','VCAN','CD14','FCGR3A'),idents = c(3,4,9))
# neu
DotPlot(myeloid_subset.harmony, features = c('CSF3R','S100A8','S100A9','FUT4','IFITM2','CMTM2'),idents = c(0,11,12))

Idents(myeloid_subset.harmony) <- myeloid_subset.harmony$seurat_clusters
myeloid_subset.harmony <- RenameIdents(myeloid_subset.harmony,`0`='Neu-CSF3R',`1`='Macro-CD163',`2`='Macro-C1QA',
                                       `3`='Mono-FN1',`4`='Mono-CD14',`5`='Macro-GPNMB',`6`='Mast-CPA3',`7`='Macro-RGS1',
                                       `8`='DC-CD1C',`9`='Mono-CD16',`10`='Macro-SPP1',`11`='Neu-CSF3R',`12`='IFN-response',
                                       `13`='DC-LAMP3')
myeloid_subset.harmony$Subtype <- Idents(myeloid_subset.harmony)
cols <- c(pal_npg("nrc",alpha = 0.5)(10),brewer.pal(12,'Set3')[c(2,8,9,10)])
cols <- cols[c(9,13,2,3,6,8,11,4,5,10,1,7,12)]
myeloid_subset.harmony$Subtype <- factor(myeloid_subset.harmony$Subtype, 
                                         levels= c("DC-CD1C","DC-LAMP3","Macro-CD163","Macro-C1QA","Macro-GPNMB",
                                                   "Macro-RGS1","Macro-SPP1","Mono-FN1","Mono-CD14","Mono-CD16",
                                                   "Neu-CSF3R","Mast-CPA3","IFN-response"),ordered = TRUE)
p_major <- DimPlot(myeloid_subset.harmony,label = TRUE,pt.size = 0.5,group.by = 'Major_2')
p_sub <- DimPlot(myeloid_subset.harmony,label = TRUE,cols = cols,pt.size = 0.5,group.by = 'Subtype',repel = T)

dot.genes <- c('CD1C','FCER1A','CLEC10A','LAMP3','CCR7','IDO1', # DCs
               'CD163','SELENOP','STAB1','C1QA','C1QC','APOE',
               'GPNMB','APOC1','CCL18','RGS1','GPR183','CCL3','SPP1','CXCL5','SLAMF9', # Macro
               'FN1','VCAN','EREG','FCN1','S100A8','S100A9','CSF3R', # mono
               'CPA3','TPSAB1','TPSB2' # Mast
)
p_dot <- angrycell::DotPlot2(object = subset(myeloid_subset.harmony,idents = 'IFN-response',invert=TRUE),
                             features = rev(dot.genes),group.by = 'Subtype')
pdf('~/Analysis/sub_myeloid/plot/myeloid_subtype.major.pdf',width = 7,height = 6)
print(p_major)
dev.off()

saveRDS(myeloid_subset.harmony,'~/Analysis/sub_myeloid/myeloid.sub.anno.rds')


######################### distinguish M1/M2
M1 <- c('FCGR1A','IDO1','SOCS1','CXCL10','CD40','Mono_Macro0','Mono_Macro6','ITGAX','TNF')
M2 <- c('CD163','MRC1','ARG1','SELENOP','DAB1','STAB1')

DotPlot(myeloid_subset.harmony,features = c(M1,M2),idents = grep('Macro-',levels(myeloid_subset.harmony),value = TRUE),
        cols = c('lightgrey','red'),col.min = -1,col.max = 1) +theme(axis.text.x = element_text(angle = 90))
myeloid_subset.harmony <- AddModuleScore(myeloid_subset.harmony,features = list(M1),name = 'Macro_M1')
myeloid_subset.harmony <- AddModuleScore(myeloid_subset.harmony,features = list(M2),name = 'Macro_M2')

FeaturePlot(myeloid_subset.harmony,features = 'Macro_M11',min.cutoff = 0)
FeaturePlot(myeloid_subset.harmony,features = 'Macro_M21',min.cutoff = 0)

macro_myeloid <- subset(myeloid_subset.harmony,idents=grep('Macro',levels(myeloid_subset.harmony),value = TRUE))
macro_myeloid$Subtype <- as.character(macro_myeloid$Subtype) 
ggplot(macro_myeloid@meta.data,aes(x=Subtype,y=Macro_M11,fill=Subtype)) + geom_boxplot() + labs(title='M1 Score')
ggplot(macro_myeloid@meta.data,aes(x=Subtype,y=Macro_M21,fill=Subtype)) + geom_boxplot() + labs(title='M2 Score')