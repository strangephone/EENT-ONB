########################################## Trajectory 
library(monocle)
library(Seurat)
library(dplyr)
myeloid_subset.harmony <- readRDS('~/Analysis/sub_NKT/myeloid.sub.anno.rds')
## Mono/Macrophage
Mono_Macro.monocle <- subset(myeloid_subset.harmony,idents = grep('^Mono|^Macro',levels(myeloid_subset.harmony),value = TRUE))

expr_matrix <- Mono_Macro.monocle@assays$RNA@counts
sample.info <- Mono_Macro.monocle@meta.data
gene_anno <- data.frame(gene_short_name=rownames(Mono_Macro.monocle))
rownames(gene_anno) <- gene_anno$gene_short_name

pd <- new("AnnotatedDataFrame", data = sample.info)
fd <- new("AnnotatedDataFrame", data = gene_anno)
cds_Mono_Macro <- newCellDataSet(expr_matrix,phenoData = pd,featureData = fd,expressionFamily = negbinomial.size())

cds_Mono_Macro <- estimateSizeFactors(cds_Mono_Macro)
cds_Mono_Macro <- estimateDispersions(cds_Mono_Macro)

cds_Mono_Macro <- detectGenes(cds_Mono_Macro, min_expr = 0.1) #这一操作会在fData(cds)中添加一列num_cells_expressed
print(head(fData(cds_Mono_Macro)))
expressed_genes <- row.names(subset(fData(cds_Mono_Macro),num_cells_expressed >= 10))
cds_Mono_Macro <- cds_Mono_Macro[expressed_genes,]

### monoce hvg 
disp_table <- dispersionTable(cds_Mono_Macro)
disp.genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)

ordered_genes <- disp.genes$gene_id
cds_Mono_Macro <- setOrderingFilter(cds_Mono_Macro, ordered_genes)
# plot_ordering_genes(cds_Mono_Macro)

cds_Mono_Macro <- reduceDimension(cds_Mono_Macro, max_components = 2,method = 'DDRTree')
cds_Mono_Macro <- orderCells(cds_Mono_Macro)
saveRDS(cds_Mono_Macro,'~/Analysis/sub_NKT/Trajectory/Mono_Macro/Mono_Macro.monocle.hvg.rds')

cds_Mono_Macro <- readRDS('~/Analysis/sub_NKT/Trajectory/Mono_Macro/Mono_Macro.monocle.hvg.rds')
plot_cell_trajectory(cds_Mono_Macro,color_by = 'Subtype',cell_size = 0.5) + scale_color_d3() 
plot_cell_trajectory(cds_Mono_Macro,color_by = 'Subtype',cell_size = 0.5) + scale_color_d3() + facet_wrap(~Subtype,ncol = 4)

### seurat hvg 
Mono_Macro.monocle <- FindVariableFeatures(Mono_Macro.monocle,nfeatures = 1000)
ordered_genes <- VariableFeatures(Mono_Macro.monocle)
cds_Mono_Macro <- setOrderingFilter(cds_Mono_Macro, ordered_genes)
# plot_ordering_genes(cds_Mono_Macro)

cds_Mono_Macro <- reduceDimension(cds_Mono_Macro, max_components = 2,method = 'DDRTree')
cds_Mono_Macro <- orderCells(cds_Mono_Macro)
saveRDS(cds_Mono_Macro,'~/Analysis/sub_NKT/Trajectory/Mono_Macro/Mono_Macro.seurat.hvg.rds')
cds_Mono_Macro <- readRDS('~/Analysis/sub_NKT/Trajectory/Mono_Macro/Mono_Macro.seurat.hvg.rds')
plot_cell_trajectory(cds_Mono_Macro,color_by = 'Subtype',cell_size = 0.5) + scale_color_d3() 
plot_cell_trajectory(cds_Mono_Macro,color_by = 'Subtype',cell_size = 0.5) + scale_color_d3() + facet_wrap(~Subtype,ncol = 4)

### seurat FindMarkers 
deg_Mono.Macro <- FindAllMarkers(Mono_Macro.monocle,only.pos = TRUE)
write.csv(deg_Mono.Macro,'~/Analysis/sub_NKT/Trajectory/Mono_Macro/Mono_Macro.Seurat.deg.csv')
ordered_genes <- unique(deg_Mono.Macro[deg_Mono.Macro$avg_logFC>=0.5,]$gene)
cds_Mono_Macro <- setOrderingFilter(cds_Mono_Macro, ordered_genes)

cds_Mono_Macro <- reduceDimension(cds_Mono_Macro, max_components = 2,method = 'DDRTree')
cds_Mono_Macro <- orderCells(cds_Mono_Macro)
saveRDS(cds_Mono_Macro,'~/Analysis/sub_NKT/Trajectory/Mono_Macro/Mono_Macro.seurat.FindMarkers.rds')
cds_Mono_Macro <- readRDS('~/Analysis/sub_NKT/Trajectory/Mono_Macro/Mono_Macro.seurat.FindMarkers.rds')
plot_cell_trajectory(cds_Mono_Macro,color_by = 'Subtype',cell_size = 0.5) + scale_color_d3() 
plot_cell_trajectory(cds_Mono_Macro,color_by = 'Subtype',cell_size = 0.5) + scale_color_d3() + facet_wrap(~Subtype,ncol = 4)
