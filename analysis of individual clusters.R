### The cells in my single cell analysis do not cluster according to their treatment groups. In order to differences between treatment groups in individual clusters, I will separate raw count data from individual clusters (as identified by prior clusterning analysis) and then perform normalisation, dimensionaloty reduction and differential expression analysis between cells in different treatment groups (WT, PDK)

library(Seurat)
library(ggplot2)
library(scTransform)
library(SingleCellExperiment)

# From prior analysis are extracted cells which correspond to certain cluster (in this case cluster 6)
c_0 <- data.frame(cluster_number = seu@meta.data$seurat_clusters, cells=colnames(seu@assays$SCT), class=seu@meta.data$class)
c_0 <- dplyr::filter(c_0, cluster_number == 6)
c_0 <- counts[, (colnames(counts) %in% c_0$cells)]
dim(c_0)

# Create an object which contains information about the treament group of each cell
new_class <- data.frame(class=new$class, cells=rownames(new))
new_class <- dplyr::filter(new_class, cells %in% colnames(c_0))
dim(new_class)

# Create Seurat object and normalise data using scTransform
seu.c_0 <- CreateSeuratObject(counts = c_0, min.cells = 3, min.features = 200)
seu.c_0 <- SCTransform(object = seu.c_0, verbose = FALSE)


# Dimensionality reduction and clustering 
set.seed(123)
seu.c_0 <- RunPCA(seu.c_0, verbose = FALSE)
seu.c_0 <- FindNeighbors(seu.c_0, dims = 1:30)
seu.c_0 <- FindClusters(seu.c_0, resolution = 1, verbose = FALSE)
seu.c_0 <- RunUMAP(seu.c_0, dims = 1:30)
DimPlot(seu.c_0, label = TRUE)

# Add metadata to the seurat object. This metadata contains information about the treament group. Double-check if the cells in the metadata object are the same as the cells in the seurat object
new_class <- as.matrix(new_class)
all.equal(new_class[,2], colnames(seu.c_0@assays$SCT))
new_class <- data.frame(new_class= new_class[,1], row.names=new_class[,2])
seu.c_0 <- AddMetaData(object=seu.c_0, metadata=new_class)

# Find differentially expressed genes between treatment groups (PDK, WT)
markers <- FindMarkers(seu.c_0, group.by='new_class', ident.1 = 'PDK', ident.2='WT', logfc.threshold=0.001)


## visualisation of results with volcano plot
# add a column of NAs
markers$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
markers$diffexpressed[markers$avg_logFC > 0.2 & markers$p_val < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
markers$diffexpressed[markers$avg_logFC < -0.2 & markers$p_val < 0.05] <- "DOWN"

# Re-plot but this time color the points with "diffexpressed"
p <- ggplot(data=markers, aes(x=avg_logFC, y=-log10(p_val), col=diffexpressed)) + geom_point() + theme_minimal()

# Add lines as before...
p2 <- p + geom_vline(xintercept=c(-0.2, 0.2), col="red") +
        geom_hline(yintercept=-log10(0.05), col="red")

# Create a new column "delabel" to de, that will contain the name of genes differentially expressed (NA in case they are not)
markers$delabel <- NA
markers $delabel[markers$diffexpressed != "NO"] <- rownames(markers)[markers$diffexpressed != "NO"]

library(ggrepel)
# plot adding up all layers we have seen so far
ggplot(data=markers, aes(x=avg_logFC, y=-log10(p_val), col=diffexpressed, label=delabel)) +
        geom_point() + 
        theme_minimal() +
        ggtitle("Cluster 6") +
        geom_text_repel() +
        scale_color_manual(values=c("blue", "black", "red")) +
        geom_vline(xintercept=c(-0.2, 0.2), col="red") +
        geom_hline(yintercept=-log10(0.05), col="red")