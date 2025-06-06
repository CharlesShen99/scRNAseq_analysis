# scRNAseq_analysis
seurat_obj <- CreateSeuratOBject(counts = raw_counts, meta.data = metadata)
# This retrieve raw counts matrix
counts <- GetAssayData(seurat_obj, assay = "RNA", layer = "counts")
# Normalize data
seurat_obj <- NormalizeData(seurat_obj)
data <- GetAssayData(seurat_obj, assay = "RNA", slot = "data")
# Scale data so each gene mean = 0 and variance is 1
seurat_obj <- ScaleData(seurat_obj)
scaled_data <- GetAssayData(seurat_obj, slot = "scale.data")
# Run PCA
seurat_obj <- RunPCA(seurat_obj)
