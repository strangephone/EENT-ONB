############################################## 3. infercnv 
### infercnv data prepare
library(Seurat)
geneOrder <- read.table('geneOrdering.txt',header = F)
########################################## fibroblast as ref
OM2 <- readRDS('OM2.rds')
Idents(OM2) <- OM2$celltype_2
OM2_fibro <- names(OM2$orig.ident[OM2$celltype_2=='Fibro'])
ref_fibro <- sample(OM2_fibro,800)
OM2_fibro <- OM2_fibro[!(OM2_fibro%in% ref_fibro)]
spike_in_fibro <- sample(OM2_fibro,500)
length(intersect(ref_fibro,spike_in_fibro))

obj_tmp <- readRDS('2-reduction_clustering_anno/merged.anno_v1.rds')
Idents(obj_tmp) <- obj_tmp$major_type

inter_gene <- intersect(rownames(obj_tmp),as.character(geneOrder$V1))
all(inter_gene %in% rownames(OM2)) # --> TRUE
obj_tmp <- subset(obj_tmp,features = inter_gene)
OM2 <- subset(OM2,features=inter_gene)
all(rownames(OM2) %in% rownames(obj_tmp)) ## --> TRUE

Epi <- subset(obj_tmp, idents = c('Epi'))
Epi$Groups <- 'Epithelial'

all(rownames(OM2)==rownames(Epi))

ref_fibro.obj <- subset(OM2,cells=ref_fibro) # extract ref cells (800)
ref_fibro.obj$Groups <- 'OM2_fibroblast'
spike_in_fibro.obj <- subset(OM2,cells=spike_in_fibro) # extract spike in cells (500)
spike_in_fibro.obj$Groups <- 'Spike-in fibroblast'

################### ONB samples infcnv data preparation
samples <- list.files('~/project/',pattern = 'ONB')

for(s in samples){
  tmp <- subset(Epi, cells = names(Epi$orig.ident[Epi$orig.ident==s]))
  input <- merge(tmp,list(spike_in_fibro.obj,ref_fibro.obj))
  cellcounts <- as.matrix(input[['RNA']]@counts)
  write.table(cellcounts,paste0('~/Analysis/3-infercnv/sample_split/',s,'.counts.mtx.txt'),sep = '\t',quote = F)
  
  meta <- data.frame(colnames(input),input$Groups)
  write.table(meta, paste0('~/Analysis/3-infercnv/sample_split/',s,'.meta.txt'),sep = '\t',quote = F,row.names = F)
}
### run infercnv (eg. ONB-199N) 
dyn.load('~/Software/JAGS-4.3.0/lib/libjags.so.4')
library(infercnv)
library(Seurat)
infercnv_obj = CreateInfercnvObject(raw_counts_matrix="ONB-199N.counts.mtx.txt", # 可以直接提供矩阵对象
                                    annotations_file="ONB-199N.meta.txt",
                                    delim="\t",
                                    gene_order_file="geneOrdering.txt",
                                    ref_group_names= 'OM2_fibroblast')

# perform infercnv operations to reveal cnv signal
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                             out_dir="ONB-199N.denoise",  # 输出文件夹
                             cluster_by_groups=F,   # 是否根据细胞注释文件的分组对肿瘤细胞进行分群，目的是通过cnv分群区分良恶性细胞，以及不同的变异类型，所以选择False
                             k_obs_groups = 2,
                             num_threads = 40,
                             denoise=T, #去噪
                             HMM=F, #选择T会有HMM CNV得分判定结果，选择F会加快运算速度
                             output_format = "pdf")
#### replot CNV
cnv.path <- '~/Analysis/3-infercnv/sample_split'
## ONB-311, k=3
ONB.311 <- readRDS(paste0(cnv.path,'/ONB-311.denoise/run.final.infercnv_obj'))
plot_cnv(ONB.311,output_filename = 'ONB-311',cluster_by_groups = F,
         out_dir = paste0(cnv.path,'/0-replot/ONB-311'),k_obs_groups = 3,output_format = 'pdf')

## ONB-860 ,k =3
ONB.860 <- readRDS(paste0(cnv.path,'/ONB-860.denoise/run.final.infercnv_obj'))
plot_cnv(ONB.860,output_filename = 'ONB-860',cluster_by_groups = F,
         out_dir = paste0(cnv.path,'/0-replot/ONB-860'),k_obs_groups = 3,output_format = 'pdf')

## ONB-946, k = 5
ONB.946 <- readRDS(paste0(cnv.path,'/ONB-946.denoise/run.final.infercnv_obj'))
plot_cnv(ONB.946,output_filename = 'ONB-946',cluster_by_groups = F,
         out_dir = paste0(cnv.path,'/0-replot/ONB-946'),k_obs_groups = 5,output_format = 'pdf')
## ONB-983N, k = 3
ONB.983N <- readRDS(paste0(cnv.path,'/ONB-983N.denoise/run.final.infercnv_obj'))
plot_cnv(ONB.983N,output_filename = 'ONB-983N',cluster_by_groups = F,
         out_dir = paste0(cnv.path,'/0-replot/ONB-983N'),k_obs_groups = 3,output_format = 'pdf')

## ONB-599, k = 5
ONB.599 <- readRDS(paste0(cnv.path,'/ONB-599.denoise/run.final.infercnv_obj'))
plot_cnv(ONB.599,output_filename = 'ONB-599',cluster_by_groups = F,
         out_dir = paste0(cnv.path,'/0-replot/ONB-599'),k_obs_groups = 5,output_format = 'pdf')


############################ distinguish basal and neural
library(harmony)
sub_malignant <- readRDS('~/Analysis/malignant/malignant.v1.rds')
ONB_tumor.harmony <- sub_malignant
ONB_tumor.harmony@meta.data <- ONB_tumor.harmony@meta.data[,1:5]
ONB_tumor.harmony <- NormalizeData(ONB_tumor.harmony)
ONB_tumor.harmony <- FindVariableFeatures(ONB_tumor.harmony,nfeatures = 1500,selection.method = 'vst')
ONB_tumor.harmony <- ScaleData(ONB_tumor.harmony,vars.to.regress = c('nFeature_RNA','nCount_RNA','percent.mt'))
ONB_tumor.harmony <- RunPCA(ONB_tumor.harmony)
ONB_tumor.harmony<- RunHarmony(ONB_tumor.harmony, group.by.vars = 'orig.ident',plot_convergence = TRUE,
                               kmeans_init_nstart=20, kmeans_init_iter_max=100)

ONB_tumor.harmony <- RunUMAP(ONB_tumor.harmony,dims = 1:nPCs,reduction = 'harmony')
ONB_tumor.harmony <- FindNeighbors(ONB_tumor.harmony,dims = 1:nPCs,reduction = 'harmony')
ONB_tumor.harmony <- FindClusters(ONB_tumor.harmony,resolution = 0.2)

p1 <- DimPlot(ONB_tumor.harmony,label = TRUE)
p2 <- DimPlot(ONB_tumor.harmony,group.by = 'orig.ident')
plot_grid(p1,p2,ncol = 1)
FeaturePlot(ONB_tumor.harmony,features = 'UBE2C')
# all.markers_ONB.harmony <- FindAllMarkers(ONB_tumor.harmony,only.pos = TRUE)

all.markers_ONB_tumor.harmony <- FindAllMarkers(ONB_tumor.harmony,only.pos = TRUE)
write.csv(all.markers_ONB_tumor.harmony,'~/Analysis/malignant/all.markers_ONB.tumor.csv')
top20_ONB_tumor.harmony <- all.markers_ONB_tumor.harmony %>% group_by(cluster) %>% top_n(20,wt=avg_logFC)

cell.distribution <- as.data.frame.array(table(ONB_tumor.harmony$orig.ident,ONB_tumor.harmony$seurat_clusters))
write.table(cell.distribution,'~/Analysis/malignant/cell.distrubution.res0.3.txt',quote = F,sep = '\t')

basal <- c('IFT57','MTHFD1L','ALX1','FZD6','CELSR1','GRHL2','KIF20B','LMO4','TEAD2','BMP4','PRICKLE1',
           'ARHGEF2','BIRC5','PBK','HMGA2','CDK2','AURKB','CDK1','CDC6','CCNA2','HELLS','RBBP8','LZTS2',
           'KNTC1','LRRCC1','ANLN','CEP55')
neural <- c('UNCX','NEUROD1','NKX2-2','SOX9','FABP7','S100B','SOX10','KCNQ2','PCDHAC2','PCDHAC1','FZD9','MYT1',
            'GNG8','HES6','OLFM1','CACNA1B','STX1A','SYN1','SNAP25','VAMP2','NRXN2','RAB3A','PTPRN2','STXBP1',
            'LIN7B')
ONB_tumor.harmony <- AddModuleScore(ONB_tumor.harmony,features = list(basal),name = 'Basal')
ONB_tumor.harmony <- AddModuleScore(ONB_tumor.harmony,features = list(neural),name = 'Neural')

p_basal <- FeaturePlot(ONB_tumor.harmony,features = 'Basal1',min.cutoff = 0,order = T)
p_neural <- FeaturePlot(ONB_tumor.harmony,features = 'Neural1',min.cutoff = 0,order = T)
plot_grid(p_basal,p_neural,ncol = 2)