## data preparing for cellphonedb
library(Seurat)
# ONB_tumor cells 
ONB_tumor.harmony <- readRDS('~/Analysis/malignant/malignant.rds')
ONB_tumor.harmony$celltype_cpdb <- ONB_tumor.harmony$anno_2

# T/NK cells
NKT_subset.harmony <- readRDS('~/Analysis/sub_NKT/NKT.subtype.harmony.rds')
NKT_subset.harmony$celltype_cpdb <- NKT_subset.harmony$Major_2
NKT_subset.cpdb <- subset(NKT_subset.harmony, subset = celltype_cpdb != 'LowQuality')

# Myeloid cells
myeloid_subset.harmony <- readRDS('~/Analysis/sub_myeloid/myeloid.rds')
myeloid_subset.harmony$celltype_cpdb <- myeloid_subset.harmony$Major_2
myeloid_subset.cpbd <- subset(myeloid_subset.harmony, subset = celltype_cpdb != 'IFN-response')
# DimPlot(myeloid_subset.cpbd)

# other cells
ONB_merged <- readRDS('2-reduction_clustering_anno/merged.rds')
## B cell/CAF/Endo/SMCs
Idents(ONB_merged) <- ONB_merged$major_type
others.cpdb <- subset(ONB_merged, idents = c('B','Plasma','Fibro','SMCs','Endo'))
others.cpdb$celltype_cpdb <- as.character(others.cpdb$major_type)

################ merge cell types
ONB_cpdb.merged <- merge(x=ONB_tumor.harmony,y=list(NKT_subset.cpdb,myeloid_subset.cpbd,others.cpdb))
# remove OM-2
ONB_cpdb.merged <- subset(ONB_cpdb.merged, subset = group=='ONB')
ONB_cpdb.merged@meta.data <- ONB_cpdb.merged@meta.data[,c(1:4,36)]
write.table(ONB_cpdb.merged[['RNA']]@data,'~/Analysis/interaction/TME.counts.txt',quote = F,sep = '\t')
write.table(ONB_cpdb.merged@meta.data[,'celltype_cpdb', drop=F],'~/Analysis/interaction/TME.meta.txt',quote = F,sep = '\t')
saveRDS(ONB_cpdb.merged,'~/Analysis/interaction/CellPhoneDB_data.rds')

ONB_cpdb.merged@meta.data$celltype_cpdb <- gsub('CD8|CD4','T',ONB_cpdb.merged@meta.data$celltype_cpdb)
ONB_cpdb.merged@meta.data$celltype_cpdb <- gsub('Basal|Neural-1|Neural-2|Mesenchymal-1|Mesenchymal-2|Mesenchymal-3','Malignant',ONB_cpdb.merged@meta.data$celltype_cpdb)
write.table(ONB_cpdb.merged@meta.data[,'celltype_cpdb', drop=F],'~/Analysis/interaction/TME.meta.merged.txt',quote = F,sep = '\t')
saveRDS(ONB_cpdb.merged,'~/Analysis/interaction/20220104-CellPhoneDB_data.mergedT_malignant.rds')
