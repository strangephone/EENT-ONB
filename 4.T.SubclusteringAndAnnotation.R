## T cell sub-clustering
dir.create('~/Analysis/sub_NKT')

obj.ONB <- readRDS('2-reduction_clustering_anno/merged.annotation.rds')
NKT_subset <- subset(obj.ONB,subset = major_type=='NK/T')
NKT_subset@meta.data <- NKT_subset@meta.data[,c(1:5,10:14)]
######## harmony
library(harmony)
NKT_subset.harmony <- NKT_subset
NKT_subset.harmony <- FindVariableFeatures(NKT_subset.harmony, nfeatures = '2000',selection.method = 'vst')
NKT_subset.harmony <- ScaleData(NKT_subset.harmony,vars.to.regress = c('nFeature_RNA','nCount_RNA','percent.mt'))
NKT_subset.harmony <- RunPCA(NKT_subset.harmony)
NKT_subset.harmony<- RunHarmony(NKT_subset.harmony, group.by.vars = 'orig.ident',plot_convergence = TRUE,
                                max.iter.harmony = 20,epsilon.cluster=-Inf,epsilon.harmony = -Inf)

NKT_subset.harmony <- RunUMAP(NKT_subset.harmony,dims = 1:10,reduction = 'harmony')
NKT_subset.harmony <- FindNeighbors(NKT_subset.harmony,dims = 1:10,reduction = 'harmony')
NKT_subset.harmony <- FindClusters(NKT_subset.harmony) # default 0.8
NKT_subset.harmony <- FindClusters(NKT_subset.harmony, resolution = 1)
NKT_subset.harmony <- FindClusters(NKT_subset.harmony,resolution = 1.2)

p1 <- DimPlot(NKT_subset.harmony,label = TRUE,pt.size = 0.6)
p2 <- DimPlot(NKT_subset.harmony,group.by = 'orig.ident',pt.size = 0.3)
p3 <- DimPlot(NKT_subset.harmony,group.by = 'group',pt.size = 0.6)
plot_grid(p1,p2,p3,ncol = 2)
Idents(NKT_subset.harmony) <- NKT_subset.harmony$RNA_snn_res.1.2
all.markers_NKT_subset.harmony_res1.2 <- FindAllMarkers(NKT_subset.harmony,only.pos = TRUE)
top20_NKT_subset.harmony <- all.markers_NKT_subset.harmony_res1.2 %>% filter(pct.1>0.3) %>% group_by(cluster) %>% top_n(20,wt=avg_logFC)

DotPlot(NKT_subset.harmony, features = c('CD4','CD8A','CD8B'),dot.scale = 8)
DotPlot(NKT_subset.harmony, features = c('CD3D','CD3E','CD3G','NCAM1','NCR1','FCGR3A','KLRD1'),dot.scale = 6)

angrycell::vioplot2(NKT_subset.harmony,paths = c('CD4','CD8A','CD8B'))
angrycell::vioplot2(NKT_subset.harmony,paths = c('NCAM1','NCR1'),jitter = T)

#### confirm CD4T,CD8T,NK
Idents(NKT_subset.harmony) <- NKT_subset.harmony$RNA_snn_res.1.2
#(1) feature plot
FeaturePlot(NKT_subset.harmony,features = c('CD3D','CD3E','CD4','CD8A','CD8B','NCAM1','NCR1','FCGR3A'),
            order = TRUE,ncol = 3,min.cutoff = 'q3')
#(2) vlolin  plot
angrycell::vioplot2(NKT_subset.harmony,paths = c('CD3D','CD3E','CD4','CD8A','CD8B','NCAM1','NCR1','FCGR3A'),
                    text.x.size=12,strip.text.size = 8, cell.order = levels(NKT_subset.harmony),
                    add.line=F,add.ave.point=F,jitter = TRUE)
#(3) dot plot 
DotPlot(NKT_subset.harmony, features = c('CD3D','CD3E','CD4','CD8A','CD8B','NCAM1','NCR1','FCGR3A'),scale = 6)

### cell score 
naive <- c('CCR7','SELL','LEF1')
NKT_subset.harmony <- AddModuleScore(NKT_subset.harmony,features = list(naive),name = 'Naive')
FeaturePlot(NKT_subset.harmony,features = 'Naive1',min.cutoff = 0.5,pt.size = 0.5,order = T)

Treg_suppressive <- c('CTLA4','TIGIT','ICOS')
NKT_subset.harmony <- AddModuleScore(NKT_subset.harmony,features = list(Treg_suppressive),name = 'Treg_suppressive')
FeaturePlot(NKT_subset.harmony,features = 'Treg_suppressive1',min.cutoff = 0.5,pt.size = 0.5,order = T)

nk <- c('FCGR3A','NCAM1','NCR1','KLRD1')
NKT_subset.harmony <- AddModuleScore(NKT_subset.harmony,features = list(nk),name = 'NK')
FeaturePlot(NKT_subset.harmony,features = 'NK1',min.cutoff = 0.5,pt.size = 0.5,order = T)

gdT_1 <- arrange(C15_16,avg_logFC) %>% head(10)
gdT_1 <- as.character(rownames(gdT_1))

gdT_2 <- arrange(C15_16,avg_logFC) %>% tail(10)
gdT_2 <- as.character(rownames(gdT_2))

NKT_subset.harmony <- AddModuleScore(NKT_subset.harmony,features = list(gdT_1),name = 'gdT_1')
FeaturePlot(NKT_subset.harmony,features = 'gdT_11',min.cutoff = 0.5,pt.size = 0.5,order = T)

NKT_subset.harmony <- AddModuleScore(NKT_subset.harmony,features = list(gdT_2),name = 'gdT_2')
FeaturePlot(NKT_subset.harmony,features = 'gdT_21',min.cutoff = 0.5,pt.size = 0.5,order = T)

gdT <- c('CD3D', 'CD3E', 'TRDC', 'TRGC1', 'TRGC2')
# gdT <- as.character(rownames(gdT))
NKT_subset.harmony <- AddModuleScore(NKT_subset.harmony,features = list(gdT),name = 'gdT')
FeaturePlot(NKT_subset.harmony,features = 'gdT1',min.cutoff = 'q80',pt.size = 0.5,order = T)
##### annotation

Idents(NKT_subset.harmony) <- NKT_subset.harmony$RNA_snn_res.1.2
NKT_subset.harmony <- RenameIdents(NKT_subset.harmony,`0`='CD4-Tcm',`1`='CD8-Tem',`2`='NK-1',
                                   `3`='CD8-Tem',`4`='CD8-Teff',`5`='CD8-Trm',
                                   `6`='CD4-Tn',`7`='CD8-Teff',`8`='CD4-Treg_resting',`9`='NK-2',
                                   `10`='CD4-Treg_suppressive',`11`='LowQuality',`12`='NK-1',`13`='CD4-Tcm',
                                   `14`='CD8-Tex',`15`='CD8-GZMK',`16`='CD8-Teff')
NKT_subset.harmony$celltype <- Idents(NKT_subset.harmony)
p1 <- DimPlot(NKT_subset.harmony,label = TRUE,group.by = 'RNA_snn_res.1.2',pt.size = 0.6)
pm <- DimPlot(NKT_subset.harmony,label = TRUE,group.by = 'celltype',pt.size = 0.6,repel = T)
plot_grid(p1,pm)

saveRDS(NKT_subset.harmony,'~/Analysis/sub_NKT/NKT.subtype.harmony.rds')

NKT_subset.harmony$subtype <- as.character(Idents(NKT_subset.harmony))
NKT_subset.harmony <- FindClusters(NKT_subset.harmony,resolution = 2)
DimPlot(NKT_subset.harmony,group.by = 'RNA_snn_res.2',label = T)
NKT_subset.harmony@meta.data$subtype2 <- ifelse(NKT_subset.harmony$RNA_snn_res.2==23,'CD4-Tfh',NKT_subset.harmony@meta.data$subtype)
p1 <- DimPlot(NKT_subset.harmony,label = TRUE,group.by = 'RNA_snn_res.1.2',pt.size = 0.6)
pm <- DimPlot(NKT_subset.harmony,label = TRUE,group.by = 'subtype2',pt.size = 0.6,repel = T)
plot_grid(p1,pm,ncol = 2)
########## featureplot(gene set)
# naive
naive <- c('CCR7','LEF1','SELL')
NKT_subset.harmony <- AddModuleScore(NKT_subset.harmony,features = list(naive),name = 'naive.scores')
p1 <- FeaturePlot(NKT_subset.harmony,features = 'naive.scores1',order = TRUE,cols = c('grey','#E41A1C'),
                  min.cutoff = 'q50',max.cutoff = 2,pt.size = 0.1) + labs(title = expression(atop('Naive T cells','(CCR7,LEF1,SELL)')))
# NK
NK.genes <- c('NCAM1','NCR1','TYROBP','KLRD1','KLRF1')
NKT_subset.harmony <- AddModuleScore(NKT_subset.harmony,features = list(NK.genes),name = 'NK.scores')
p2 <- FeaturePlot(NKT_subset.harmony,features = 'NK.scores1',order = TRUE,cols = c('grey','#E41A1C'),
                  min.cutoff = 'q50',max.cutoff = 2,pt.size = 0.1) + labs(title = expression(atop('NK','(NCAM1,NCR1,TYROBP,KLRD1,KLRF1)')))
# CD4 Treg
Treg.genes <- c('CD4','FOXP3','IL2RA','TNFRSF4')
NKT_subset.harmony <- AddModuleScore(NKT_subset.harmony,features = list(Treg.genes),name = 'Treg.scores')
p3 <- FeaturePlot(NKT_subset.harmony,features = 'Treg.scores1',order = TRUE,cols = c('grey','#E41A1C'),
                  min.cutoff = 'q50',max.cutoff = 2,pt.size = 0.1) + labs(title = expression(atop('Tregs','(CD4,FOXP3,IL2RA,TNFRSF4)')))
# CD4 Tfh
cd4.Tfh.genes <- c('CD4','IL21','ICOS','CXCL13')
NKT_subset.harmony <- AddModuleScore(NKT_subset.harmony,features = list(cd4.Tfh.genes),name = 'CD4Tfh.scores')
p4 <- FeaturePlot(NKT_subset.harmony,features = 'CD4Tfh.scores1',order = TRUE,cols = c('grey','#E41A1C'),
                  min.cutoff = 'q50',max.cutoff = 2,pt.size = 0.1) + labs(title = expression(atop('CD4+ Tfh','(CD4,IL21,ICOS,CXCL13)')))
# CD4 Tcm
cd4.Tcm.genes <- c('CD4','IL7R','ANXA1')
NKT_subset.harmony <- AddModuleScore(NKT_subset.harmony,features = list(cd4.Tcm.genes),name = 'CD4Tcm.scores')
p5 <- FeaturePlot(NKT_subset.harmony,features = 'CD4Tcm.scores1',order = TRUE,cols = c('grey','#E41A1C'),
                  min.cutoff = 'q50',max.cutoff = 2,pt.size = 0.1) + labs(title = expression(atop('CD4+ Tcm','(CD4,IL7R,ANXA1)')))
# CD8 Tem
cd8.Tem.genes <- c('CD8A','CD8B','GZMK','CCL4')
NKT_subset.harmony <- AddModuleScore(NKT_subset.harmony,features = list(cd8.Tem.genes),name = 'CD8Tem.scores')
p6 <- FeaturePlot(NKT_subset.harmony,features = 'CD8Tem.scores1',order = TRUE,cols = c('grey','#E41A1C'),
                  min.cutoff = 'q50',max.cutoff = 2,pt.size = 0.1) + labs(title = expression(atop('CD8+ Tem','(CD8A,CD8B,GZMK,CCL4)')))
# CD8 Trm
cd8.Trm.genes <- c('CD8A','CD8B','ZNF683','ITGA1')
NKT_subset.harmony <- AddModuleScore(NKT_subset.harmony,features = list(cd8.Trm.genes),name = 'CD8Trm.scores')
p7 <- FeaturePlot(NKT_subset.harmony,features = 'CD8Trm.scores1',order = TRUE,cols = c('grey','#E41A1C'),
                  min.cutoff = 'q50',max.cutoff = 2,pt.size = 0.1) + labs(title = expression(atop('CD8+ Trm','(CD8A,CD8B,ZNF683,ITGA1)')))
# CD8 Teff
cd8.Teff.genes <- c('CD8A','CD8B','GZMH','CX3CR1','FGFBP2')
NKT_subset.harmony <- AddModuleScore(NKT_subset.harmony,features = list(cd8.Teff.genes),name = 'CD8Teff.scores')
p8 <- FeaturePlot(NKT_subset.harmony,features = 'CD8Teff.scores1',order = TRUE,cols = c('grey','#E41A1C'),
                  min.cutoff = 'q50',max.cutoff = 2,pt.size = 0.1) + labs(title = expression(atop('CD8+ Teff','(CD8A,CD8B,GZMH,CX3CR1,FGFBP2)')))
# CD8 Tex
cd8.Tex.genes <- c('CD8A','CD8B','PDCD1','CTLA4','LAG3')
NKT_subset.harmony <- AddModuleScore(NKT_subset.harmony,features = list(cd8.Tex.genes),name = 'CD8Tex.scores')
p9 <- FeaturePlot(NKT_subset.harmony,features = 'CD8Tex.scores1',order = TRUE,cols = c('grey','#E41A1C'),
                  min.cutoff = 0,max.cutoff = 2,pt.size = 0.1) + labs(title = expression(atop('CD8+ Tex','(CD8A,CD8B,PDCD1,CTLA4,LAG3)')))
Cairo::CairoPDF(file = '~/Analysis/sub_NKT/avg.Featureplot.pdf',width = 20,height = 15)
plot_grid(p1,p2,p3,p4,p5,p6,p7,p8,p9,ncol = 3)
dev.off()
saveRDS(NKT_subset.harmony,'~/Analysis/sub_NKT/NKT.subtype.harmony_20211216.rds')

cols <- pal_d3("category20")(20)
color_palette <- read.table('~/bin/my_scripts/R_func/color_palette.txt')
color_palette <- as.character(color_palette$V1)[1:14]
color_palette <- color_palette[c(1:4,6,5,7:14)]
color_palette <- color_palette[c(1:8,13,10:12,9)]

names(color_palette) <- levels(NKT_subset.harmony$SubType)
Cairo::CairoPDF('~/Analysis/sub_NKT/1_cell_anno/NKT.anno.v1.pdf',width = 11,height = 10)
DimPlot(NKT_subset.harmony,group.by = 'SubType',cols = color_palette,pt.size = 0.5) + 
  theme(legend.position = "bottom")
dev.off()
