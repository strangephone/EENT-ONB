################################# SCENIC data preparation
NKT_subset.harmony <- readRDS('~/Analysis/sub_NKT/NKT.subtype.rds')
counts <- as.matrix(t(NKT_subset.harmony[['RNA']]@counts))
write.table(counts,'~/Analysis/sub_NKT/2_pyscenic/NKT.counts.tsv',sep = '\t',quote = F)
meta <- NKT_subset.harmony@meta.data
write.csv(meta,'~/Analysis/sub_NKT/2_pyscenic/NKT.meta.csv')
################################ selected genes heatmap
genes.heatmap <- c('TCF7','SELL','LEF1','CCR7', # naive
                   'IL2RA','FOXP3','IKZF2',
                   'CD28','TNFRSF14','ICOS','TNFRSF9',
                   'ZNF683','ID2','HOPX','EOMES','TBX21','ZEB2','HIF1A','TOX',
                   'IL2','GZMA','GZMK','IFNG','GNLY','PRF1','GZMB','NKG7',
                   'LAG3','PDCD1','HAVCR2','TIGIT','CTLA4',
                   'NCAM1','NCR1','XCL1','FGFBP2','PTGDS'
)
Idents(NKT_subset.harmony) <- NKT_subset.harmony$SubType
tt <- AverageExpression(NKT_subset.harmony,features = genes.heatmap,assays = 'RNA')
data_heatmap <- tt$RNA[,names(tt[[1]])[1:12]] # remove lowquality
library(pheatmap)
library(ComplexHeatmap)
anno_col <- data.frame(colnames(data_heatmap),CellType=colnames(data_heatmap),row.names = 1)
celltype_colors <- color_palette[1:12]
names(celltype_colors) <- colnames(data_heatmap)
pdf('~/Analysis/sub_NKT/3_heatmap/heatmap.v1.pdf',width = 11,height = 9)
pheatmap(data_heatmap,cluster_rows = F,cluster_cols = F,scale = 'row',breaks = seq(-1,3,length.out = 100),
         border_color = NA,gaps_col = c(5,10),gaps_row = c(4,7,11,19,27,32),show_colnames = F,
         annotation_col = anno_col,annotation_colors = list(CellType=celltype_colors)
) 
dev.off()
pdf('~/Analysis/sub_NKT/3_heatmap/heatmap.v1.legend.pdf',width = 8,height = 9)
ggplot() + annotate("point", x=1:12,y=1,shape=19, color=color_palette[1:12],size=5) +
  #  annotate("point", x=1:12,y=2,shape=16, color=color_palette[1:12],size=5) +
  #  annotate("point", x=1:12,y=3,shape=20, color=color_palette[1:12],size=5) +
  theme_classic()
dev.off()

### plot top10 regulon specificity scores
library(ggrepel)
rss.tab <- read.csv('~/Analysis/sub_NKT/2_pyscenic/jupyter/NKT.rss.csv',row.names = 1,check.names = F)
rss.tab <- t(rss.tab)

rss.plot.list <- list()
for(c in colnames(rss.tab)){
  temp <- rss.tab[,c]
  temp <- sort(temp,decreasing = TRUE) 
  rss.df <- data.frame(score=temp,rank=1:length(temp))
  rss.df$name <- rownames(rss.df)
  rss.df$group <- c(rep('g1',10),rep('g2',nrow(rss.df)-10))
  
  # ylims <- c(floor(min(rss.df$score)* 100.0) / 100.0, ceiling(max(rss.df$score)* 100.0) / 100.0)
  rss.plot.list[[c]] <-   ggplot(rss.df,aes(x=rank,y=score,color=group)) +geom_point(shape=20,size=1) + theme_classic() + 
    theme(axis.line = element_line(size=1),axis.text = element_text(size = 10,face = 'bold'),legend.position = 'none') + 
    geom_label_repel(data = rss.df[1:10,],aes(x=rank,y=score,label=name),nudge_x = 30) + 
    xlab('Rank') + ylab('Regulon Specificity Score(RSS)') + labs(title = c)
  CairoPNG(paste0('~/Analysis/sub_NKT/2_pyscenic/RSS_regulon_sepcificity_score/',c,'.RSS.png'),width = 1000,height=958)
  print(rss.plot.list[[c]])
  dev.off()
}

saveRDS(NKT_subset.harmony,'~/Analysis/sub_NKT/NKT.subtype.rds')