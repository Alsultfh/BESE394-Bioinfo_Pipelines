
options(width = 210)
library(Seurat)
library(ggplot2)
library(Signac)
library(JASPAR2020)
library(tidygraph)
library(igraph)
library(EnsDb.Mmusculus.v79)
library(BSgenome.Mmusculus.UCSC.mm10)
library(TFBSTools)
library(motifmatchr)
library(Matrix)
library(irlba)
library(chromVAR)


## iF IN DOCKER WE NEED:
# BiocManager::install("biovizBase")
# install.packages("Matrix", type = "source")
# install.packages("irlba", type = "source")
#install.packages('ggseqlogo')

# Multiome https://www.10xgenomics.com/datasets/mouse-brain-nuclei-isolated-with-chromium-nuclei-isolation-kit-saltyez-protocol-and-10x-complex-tissue-dp-ct-sorted-and-ct-unsorted-1-standard
# wget https://cf.10xgenomics.com/samples/cell-arc/2.0.2/M_Brain_Chromium_Nuc_Isolation_vs_SaltyEZ_vs_ComplexTissueDP/M_Brain_Chromium_Nuc_Isolation_vs_SaltyEZ_vs_ComplexTissueDP_filtered_feature_bc_matrix.tar.gz
# wget https://cf.10xgenomics.com/samples/cell-arc/2.0.2/M_Brain_Chromium_Nuc_Isolation_vs_SaltyEZ_vs_ComplexTissueDP/M_Brain_Chromium_Nuc_Isolation_vs_SaltyEZ_vs_ComplexTissueDP_atac_fragments.tsv.gz
# wget https://cf.10xgenomics.com/samples/cell-arc/2.0.2/M_Brain_Chromium_Nuc_Isolation_vs_SaltyEZ_vs_ComplexTissueDP/M_Brain_Chromium_Nuc_Isolation_vs_SaltyEZ_vs_ComplexTissueDP_atac_fragments.tsv.gz.tbi


plot_QC_features <- function(seurat_object, pdf_path){
  pdf(pdf_path, width =10)
      # shows the distribution of the transcripts per cells
      print(VlnPlot(
          object = seurat_object,
          features = c("nCount_RNA", "nCount_ATAC", 'nFeature_RNA',  'mitoPct', 
                      'FRiP', 'log10GenesPerUMI', "TSS.enrichment", 
                      "nucleosome_signal"),
          pt.size = 0, ncol =4))

      DefaultAssay(seurat_object) <- "ATAC"
      print(FragmentHistogram(object = seurat_object, region = 'chr1-1-10000000', group.by = 'nucleosome_group'))
      #print(TSSPlot(seurat_object, group.by = 'high.tss') + NoLegend())
      print(DensityScatter(seurat_object, x = 'nCount_ATAC', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE))
      print(ggplot(seurat_object@meta.data, aes(color=orig.ident, x=nFeature_RNA, fill= orig.ident)) + 
          geom_density(alpha = 0.2) + 
          theme_classic() +
          scale_x_log10() + 
          geom_vline(xintercept = 300) + ggtitle('GENESpercell'))
      # Visualize the distribution of genes detected per cell via boxplot
      print(ggplot(seurat_object@meta.data, aes(x=orig.ident, y=log10(nFeature_RNA), fill=orig.ident)) + 
          geom_boxplot() + 
          theme_classic() +
          theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
          theme(plot.title = element_text(hjust=0.5, face="bold")) +
          ggtitle("NCells vs NGenes"))
      # correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
      print(ggplot(seurat_object@meta.data, aes(x=nCount_RNA, y=nFeature_RNA, colour=mitoPct, group=orig.ident)) + 
          geom_point() + 
          scale_colour_gradient(low = "gray90", high = "black", limits=c(0,100)) +
          stat_smooth(method=lm) +
          scale_x_log10() + 
          scale_y_log10() + 
          theme_classic() +
          geom_vline(xintercept = 500) +
          geom_hline(yintercept = 6000) +
          geom_hline(yintercept = 250) + ggtitle(paste0('UMIvsGENESpercell  Ncell:: ', ncol(seurat_object))))
  dev.off()

}


### PARAMETERS
basepath <- ("/ibex/scratch/gonzalac/BESE394a/W5_multiome/")
#basepath <- ("/home/gonzalac/data_W5/")
### 1. Load RNA seq and ATAC data  ---------------------------------------------------------------
### all data are in /home/gonzalac/data_W5 directory (singularity can not access scratch)
counts <- Read10X(paste0(basepath,"Data/filtered_feature_bc_matrix"))
fragpath <- paste0(basepath,'Data/M_Brain_Chromium_Nuc_Isolation_vs_SaltyEZ_vs_ComplexTissueDP_atac_fragments.tsv.gz')
str(counts)

# get gene annotations for hg38
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevels(annotation) <- paste0('chr', seqlevels(annotation)) ##format seqnames containing chr correctly


### 2. Creacte Seurat object  ---------------------------------------------------------------
# first create it with the RNA data counts[['Gene Expression']]
# format correctly the rownames of the ATACseq data 
# then load the ATACseq data as an extra assay

mouse_brain <- CreateSeuratObject(
  counts = counts[['Gene Expression']],
  assay = "RNA"
)

atac_counts <- counts$Peaks ### Rownames look like chrX:162621438-162622286
grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-")) ## Create a GRanges object with the chrX and ranges 
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts) ## filter for only standard chrom (80880/80935)
atac_counts <- atac_counts[as.vector(grange.use), ]

mouse_brain[["ATAC"]] <- CreateChromatinAssay(
  counts = atac_counts,
  sep = c(":", "-"),
  fragments = fragpath,
  annotation = annotation
)

str(mouse_brain)
head(mouse_brain@meta.data)
head(mouse_brain@assays$ATAC@ranges)

### 3. Add necesary cols to Metadata to do QC  ---------------------------------------------------------------
DefaultAssay(mouse_brain) <- "ATAC"
mouse_brain <- NucleosomeSignal(mouse_brain)
mouse_brain <- TSSEnrichment(mouse_brain) ##I removed the , fast=FALSE argument
total_fragments <- CountFragments(fragments = fragpath)
rownames(total_fragments) <- total_fragments$CB
mouse_brain$fragments <- total_fragments[colnames(mouse_brain), "frequency_count"]
mouse_brain <- FRiP(object = mouse_brain, assay = 'ATAC', total.fragments = 'fragments')
mouse_brain$nucleosome_group <- ifelse(mouse_brain$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
mouse_brain$high.tss <- ifelse(mouse_brain$TSS.enrichment > 3, 'High', 'Low')

DefaultAssay(mouse_brain) <- "RNA"
mouse_brain$mitoPct <- PercentageFeatureSet(mouse_brain, pattern = "^mt-")
mouse_brain$RPSPct  <- PercentageFeatureSet(object = mouse_brain, pattern = "^Rp[sl]")
mouse_brain$log10GenesPerUMI <- log10(mouse_brain$nFeature_RNA) / log10(mouse_brain$nCount_RNA)

pdf_path <- paste0(basepath, ('Plots/5k_mouse_brain_GEX_QC_Pre.pdf'))
plot_QC_features(mouse_brain, pdf_path)
saveRDS(mouse_brain, paste0(basepath, 'Data/mouse_brain_multiome_Pre_QC.rds'))
mouse_brain <- readRDS(paste0(basepath, 'Data/mouse_brain_multiome_Pre_QC.rds'))

mouse_brain
summary(mouse_brain@meta.data)
sd(mouse_brain@meta.data)
5*(sd(mouse_brain@meta.data$nCount_RNA)) ##[1] 35075.9
5*(sd(mouse_brain@meta.data$nFeature_RNA)) ##[1] 7761.97
5*(sd(mouse_brain@meta.data$nCount_ATAC)) ##[1] 22594.5
5*(sd(mouse_brain@meta.data$nFeature_ATAC)) ##[1] 10253.58
5*(sd(mouse_brain@meta.data$TSS.enrichment)) ##[1] 7.319532
5*(sd(mouse_brain@meta.data$FRiP)) ##[1] 0.7048801

sd(mouse_brain@meta.data$nCount_RNA) ##[1] 7015.181
sd(mouse_brain@meta.data$nFeature_RNA) ##[1] 1552.394
sd(mouse_brain@meta.data$nCount_ATAC) ##[1] 4518.899
sd(mouse_brain@meta.data$nFeature_ATAC) ##[1] 2050.715
sd(mouse_brain@meta.data$TSS.enrichment) ##[1] 2.439844
sd(mouse_brain@meta.data$FRiP) ##[1] 0.23496

#          orig.ident      nCount_RNA      nFeature_RNA    nCount_ATAC    
#  SeuratProject:23990   Min.   :     4   Min.   :    4   Min.   :    17  
#                        1st Qu.:  2291   1st Qu.: 1363   1st Qu.:  1975  
#                        Median :  6018   Median : 2716   Median :  4140  
#                        Mean   :  7543   Mean   : 2747   Mean   :  5121  
#                        3rd Qu.: 10330   3rd Qu.: 3766   3rd Qu.:  7131  
#                        Max.   :130578   Max.   :12583   Max.   :110332  
#  nFeature_ATAC   nucleosome_signal nucleosome_percentile TSS.enrichment    
#  Min.   :    6   Min.   :0.0000    Min.   :0.0000        Min.   : 0.00304  
#  1st Qu.: 1019   1st Qu.:0.5858    1st Qu.:0.2500        1st Qu.: 3.26340  
#  Median : 2054   Median :0.6712    Median :0.5000        Median : 4.39007  
#  Mean   : 2495   Mean   :   Inf    Mean   :0.5001        Mean   : 4.85064  
#  3rd Qu.: 3503   3rd Qu.:0.7692    3rd Qu.:0.7500        3rd Qu.: 5.82311  
#  Max.   :28969   Max.   :   Inf    Max.   :1.0000        Max.   :63.73626  
#  TSS.percentile   fragments           FRiP        nucleosome_group  
#  Min.   :0.00   Min.   :    26   Min.   :0.1397   Length:23990      
#  1st Qu.:0.25   1st Qu.:  3796   1st Qu.:0.2933   Class :character  
#  Median :0.50   Median :  8764   Median :0.4435   Mode  :character  
#  Mean   :0.50   Mean   : 13497   Mean   :0.4900                     
#  3rd Qu.:0.75   3rd Qu.: 18066   3rd Qu.:0.6441             
#  Max.   :1.00   Max.   :529795   Max.   :1.3989                     
#    high.tss            mitoPct            RPSPct        log10GenesPerUMI
#  Length:23990       Min.   : 0.0000   Min.   : 0.0000   Min.   :0.7888  
#  Class :character   1st Qu.: 0.8524   1st Qu.: 0.8315   1st Qu.:0.8899  
#  Mode  :character   Median : 1.7136   Median : 1.2663   Median :0.9086  
#                     Mean   : 3.2920   Mean   : 1.8132   Mean   :0.9102  
#                     3rd Qu.: 3.9331   3rd Qu.: 2.2123   3rd Qu.:0.9334  
#                     Max.   :42.2763   Max.   :15.3965   Max.   :1.0000  


### 4. Filter Cells  ---------------------------------------------------------------
mouse_brain <- subset(
  x = mouse_brain,
  subset = nCount_ATAC < 22594.5 &
    nCount_RNA < 35075.9 &
    nCount_ATAC > 1000 &
    nCount_RNA > 1000 &
    nFeature_ATAC < 10253.58 &
    nFeature_RNA < 7761.97 &
    nucleosome_signal < 2 &
    TSS.enrichment > 1 & 
    FRiP > 0.2 &
    mitoPct < 15
)

plot_QC_features(mouse_brain, paste0(basepath,'Plots/5k_mouse_brain_GEX_QC_Post2.pdf'))
saveRDS(mouse_brain, paste0(basepath, 'Data/mouse_brain_multiome_Post_QC.rds'))


mouse_brain <- readRDS(paste0(basepath, 'Data/mouse_brain_multiome_Post_QC.rds'))


### 5. Process the data (normalization, scaling, PCA, clustering, UMAP, ....)   ---------------------------------------------------------------
### Normalization
### Find Var Features
### Scaling
### PCA
### FindNeigh FindClust
### UMAP
### Azimuth
DefaultAssay(mouse_brain) <- "RNA"
mouse_brain <- NormalizeData(mouse_brain)
mouse_brain <- FindVariableFeatures(mouse_brain, 
            selection.method = "vst", 
            nfeatures = 2000)
mouse_brain <- ScaleData(mouse_brain)
mouse_brain <- RunPCA(mouse_brain)
mouse_brain <- FindNeighbors(mouse_brain, dims = 1:30)
mouse_brain <- FindClusters(mouse_brain, 
            resolution = 0.4, 
            algorithm = 3, 
            cluster.name="RNA_clusters_03")
mouse_brain <- RunUMAP(mouse_brain, dims=1:30, 
            reduction = "pca", 
            reduction.name = "rna_umap")
mouse_brain <- Azimuth::RunAzimuth(mouse_brain, 
            reference = "mousecortexref")

head(mouse_brain@meta.data)

### 6. Visualize UMAP  ---------------------------------------------------------------
pdf(paste0(basepath, 'Plots/UMAP_RNA.pdf'))
DimPlot(mouse_brain, reduction = "rna_umap", group.by='RNA_clusters_03', 
        label = TRUE, repel = TRUE) + NoLegend()
DimPlot(mouse_brain, reduction = "rna_umap", group.by = "predicted.subclass", 
        label = TRUE, label.size = 3) + NoLegend()
dev.off()

DefaultAssay(mouse_brain) <- "ATAC"
mouse_brain <- FindTopFeatures(mouse_brain, min.cutoff = 5)
mouse_brain <- RunTFIDF(mouse_brain)

# install.packages('Matrix')
# install.packages("irlba")
mouse_brain <- RunSVD(mouse_brain)

pdf(paste0(basepath, 'Plots/Correlation.pdf'))
DepthCor(mouse_brain)
dev.off()

mouse_brain <- RunUMAP(mouse_brain, dims=1:30, 
                reduction = "lsi", 
                reduction.name = "atac_umap")

pdf(paste0(basepath, 'Plots/UMAP_ATAC.pdf'))
DimPlot(mouse_brain, reduction = "atac_umap", label = TRUE, repel = TRUE) + NoLegend()
DimPlot(mouse_brain, reduction = "atac_umap", group.by = "predicted.subclass",
label = TRUE, label.size = 3) + NoLegend()
dev.off()

saveRDS(mouse_brain, paste0(basepath, 'Data/mouse_brain_multiome.rds'))

### 7. Integrate  ---------------------------------------------------------------
mouse_brain <- FindMultiModalNeighbors(
  object = mouse_brain,
  reduction.list = list("pca", "lsi"), 
  dims.list = list(1:50, 1:40),
  verbose = TRUE
)

# build a joint UMAP visualization
mouse_brain <- RunUMAP(
  object = mouse_brain, 
  reduction.name = "wnn.umap", dims = 1:30,
  assay = "RNA",
  verbose = TRUE
)

saveRDS(mouse_brain, paste0(basepath, 'Data/mouse_brain_multiome_Integrated.rds'))

pdf(paste0(basepath, 'Plots/Umap_integrated.pdf'))
p1 <- DimPlot(mouse_brain, reduction = "rna_umap", label = TRUE, group.by = "predicted.subclass",
        repel = TRUE) + NoLegend()
p2 <- DimPlot(mouse_brain, reduction = "atac_umap", label = TRUE, group.by = "predicted.subclass",
        repel = TRUE) + NoLegend()
p3 <- DimPlot(mouse_brain, reduction = "wnn.umap", label = TRUE, group.by = "predicted.subclass",
        repel = TRUE) + NoLegend()
print(p1)
print(p2)
print(p3)
print(p1 | p2 | p3)
dev.off()

pdf(paste0(basepath, 'Plots/Umap_integrated_weights.pdf'))
Idents(mouse_brain) <- "predicted.subclass"
VlnPlot(mouse_brain, features = c("RNA.weight", "ATAC.weight"),  pt.size = 0, ncol = 1)
dev.off()

### 8. DA (ATAC peaks)  ---------------------------------------------------------------
#mouse_brain <- readRDS(paste0(basepath, 'Data/mouse_brain_multiome_Integrated.rds'))

DefaultAssay(mouse_brain) <- 'ATAC'
Idents(mouse_brain) <- "predicted.subclass"
unique(Idents(mouse_brain))

### DA peaks of inhibitory neurons Lamp5 against Sst,Sst Chodl,Pvalb,Vip
da_peaks_Lamp5_Sst <- FindMarkers(
  object = mouse_brain,
  ident.1 = c("Lamp5"), 
  ident.2 = c("Sst"),
  test.use = 'LR',
  latent.vars = 'nCount_ATAC'
)
saveRDS(da_peaks_Lamp5_Sst, paste0(basepath, 'Data/DA_peaks_Lamp5_Sst.rds'))


da_peaks_Lamp5_Sstchodl <- FindMarkers(
  object = mouse_brain,
  ident.1 = c("Lamp5"), 
  ident.2 = c("Sst Chodl"),
  test.use = 'LR',
  latent.vars = 'nCount_ATAC'
)
saveRDS(da_peaks_Lamp5_Sstchodl, paste0(basepath, 'Data/DA_peaks_Lamp5_Sstchodl.rds'))


da_peaks_Lamp5_Pvalb <- FindMarkers(
  object = mouse_brain,
  ident.1 = c("Lamp5"), 
  ident.2 = c("Pvalb"),
  test.use = 'LR',
  latent.vars = 'nCount_ATAC'
)
saveRDS(da_peaks_Lamp5_Pvalb, paste0(basepath, 'Data/DA_peaks_Lamp5_Pvalb.rds'))


da_peaks_Lamp5_Vip <- FindMarkers(
  object = mouse_brain,
  ident.1 = c("Lamp5"), 
  ident.2 = c("Vip"),
  test.use = 'LR',
  latent.vars = 'nCount_ATAC'
)
saveRDS(da_peaks_Lamp5_Vip, paste0(basepath, 'Data/DA_peaks_Lamp5_Vip.rds'))

#da_peaks <- readRDS(paste0(basepath, 'Data/DA_peaks.rds'))


#### Plots -------------------------------------------------------------------------
mouse_brain <- readRDS(paste0(basepath, 'Data/mouse_brain_multiome_Integrated.rds'))

### read da peaks
da_peaks_Lamp5_Sst <- readRDS(paste0(basepath, 'Data/DA_peaks_Lamp5_Sst.rds'))
da_peaks_Lamp5_Sstchodl <- readRDS(paste0(basepath, 'Data/DA_peaks_Lamp5_Sstchodl.rds'))
da_peaks_Lamp5_Pvalb <- readRDS(paste0(basepath, 'Data/DA_peaks_Lamp5_Pvalb.rds'))
da_peaks_Lamp5_Vip <- readRDS(paste0(basepath, 'Data/DA_peaks_Lamp5_Vip.rds'))

da_peaks <- da_peaks_Lamp5_Sst
da_peaks <- da_peaks[order(da_peaks$avg_log2FC, decreasing = TRUE), ]
head(da_peaks)

pdf(paste0(basepath, 'Plots/da_peaks_Lamp5_Sst.pdf'), width=12, height=20)
cowplot::plot_grid(
  VlnPlot(
    object = mouse_brain,
    assay = 'ATAC',
    features = rownames(da_peaks)[1:5],
    pt.size = 0.1,
    group.by='predicted.subclass', 
    ncol=3
  ),
  VlnPlot(
    object = mouse_brain,
    assay = 'RNA',
    features = ClosestFeature(mouse_brain, rownames(da_peaks)[1:5])$gene_name,
    pt.size = 0.1,
    group.by='predicted.subclass', 
    ncol=3
  ),
  FeaturePlot(
    object = mouse_brain,
    reduction = 'wnn.umap', 
    order=TRUE,
    features = rownames(da_peaks)[1:5],
    pt.size = 0.1,
    max.cutoff = 'q95',
    ncol=3
  ) & NoLegend(),
  FeaturePlot(
    object = mouse_brain,
    reduction = 'wnn.umap', 
    order=TRUE,
    features = ClosestFeature(mouse_brain, rownames(da_peaks)[1:5])$gene_name,
    pt.size = 0.1,
    max.cutoff = 'q95',
    ncol=3
  ) & NoLegend(),
nrow=4)
dev.off()


# Annotate peaks with his closest feature
open_Lamp5 <- rownames(da_peaks[da_peaks$avg_log2FC > 3, ])
open_Sst <- rownames(da_peaks[da_peaks$avg_log2FC < 3, ])
closest_Lamp5 <- ClosestFeature(mouse_brain, open_Lamp5)
closest_Sst <- ClosestFeature(mouse_brain, open_Sst)

# https://www.cellsignal.com/pathways/neuronal-and-glial-cell-markers
# Visualize the coverage of the peaks
pdf(paste0(basepath, 'Plots/Coverage_selected.pdf'), height=12)
CoveragePlot(
  object = mouse_brain,
  region = c('Olig1', 
             'Gfap'),
  extend.upstream = 1000,
  extend.downstream = 1000,
  ncol = 1
)

dev.off()
head(rownames(mouse_brain[["ATAC"]]))


################################################################################
# Motif analysis
################################################################################

# Motif analysis with the DA peaks
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

# add motif information
mouse_brain <- AddMotifs(
  object = mouse_brain,
  genome = BSgenome.Mmusculus.UCSC.mm10,
  pfm = pfm
)


DefaultAssay(mouse_brain) <- "ATAC"
pwm_set <- getMatrixSet(x = JASPAR2020, opts = list(species = 'Mus musculus', all_versions = FALSE))
motif.matrix <- CreateMotifMatrix(features = granges(mouse_brain), pwm = pwm_set, genome = 'mm10', use.counts = FALSE)
motif.object <- CreateMotifObject(data = motif.matrix, pwm = pwm_set)
mouse_brain <- SetAssayData(mouse_brain, assay = 'ATAC', layer = 'motifs', new.data = motif.object)



# Note that this step can take 30-60 minutes 
mouse_brain <- RunChromVAR(
  object = mouse_brain,
  genome = BSgenome.Mmusculus.UCSC.mm10
)



### read da peaks
da_peaks_Lamp5_Sst <- readRDS(paste0(basepath, 'Data/DA_peaks_Lamp5_Sst.rds'))
da_peaks_Lamp5_Sstchodl <- readRDS(paste0(basepath, 'Data/DA_peaks_Lamp5_Sstchodl.rds'))
da_peaks_Lamp5_Pvalb <- readRDS(paste0(basepath, 'Data/DA_peaks_Lamp5_Pvalb.rds'))
da_peaks_Lamp5_Vip <- readRDS(paste0(basepath, 'Data/DA_peaks_Lamp5_Vip.rds'))

list <- c(da_peaks_Lamp5_Sst,da_peaks_Lamp5_Sstchodl,da_peaks_Lamp5_Pvalb,da_peaks_Lamp5_Vip)

wkg_peaks <- da_peaks_Lamp5_Vip 
wkg_peaks_name <- deparse(substitute(da_peaks_Lamp5_Vip))
extracted_name <- sub("da_peaks_", "", wkg_peaks_name)
wkg_name <- paste0(basepath, 'Plots/motifs_',extracted_name,'.pdf')

top_peaks <- rownames(wkg_peaks[wkg_peaks$p_val < 0.005, ])

enriched.motifs <- FindMotifs(
  object = mouse_brain,
  features = top_peaks
)

pdf(wkg_name)
MotifPlot(
  object = mouse_brain,
  motifs = head(rownames(enriched.motifs))
)
dev.off()




################################################################################
# Generate a RNA activity matrix based on the ATAC-seq data
################################################################################

# We can create a proxy of the gene expression from the ATAC-seq data using  GeneActivity function. 
# We can also create a proxy of the gene expression from the ATAC-seq data using the chromVAR package. 
# This package uses the motif accessibility to infer the gene expression. 
# We can use the motif models from JASPAR2020 to perform this analysis.

gene.activities <- GeneActivity(mouse_brain)
mouse_brain[['RNA_ACTIVITY']] <- CreateAssayObject(counts = gene.activities)
mouse_brain <- NormalizeData(
  object = mouse_brain,
  assay = 'RNA_ACTIVITY',
  normalization.method = 'LogNormalize',
  scale.factor = median(mouse_brain$nCount_RNA)
)

DefaultAssay(mouse_brain) <- 'RNA'
RNA_plot <- FeaturePlot(
    object = mouse_brain,
    order=TRUE,
    features =c("Sema5a","Dennd4a","Nkain1"),
    pt.size = 0.1,
    max.cutoff = 'q95',
    ncol = 3
  )& NoLegend()

DefaultAssay(mouse_brain) <- 'RNA_ACTIVITY'
RNA_Activity_plot <- FeaturePlot(
    object = mouse_brain,
    order=TRUE,
    features =c("Sema5a","Dennd4a","Nkain1"),
    pt.size = 0.1,
    max.cutoff = 'q95',
    ncol = 3
  )& NoLegend()

pdf(paste0(basepath, 'Plots/RNA_comparison.pdf'), height=12)
cowplot::plot_grid(
  RNA_plot,
  RNA_Activity_plot,
nrow=2)
dev.off()
