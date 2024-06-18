#------ By time point Script - By Miquel Gratacós - Supervised by Laia Bassaganyas

## ----pkgs, message=FALSE, warning=FALSE, include=FALSE------------------------------------------------------------------
library(Seurat)
library(Matrix)
library(singleCellTK)
library(dplyr)
library(patchwork)
library(ggplot2)
library(purrr)
library(xlsx)
library(rio)
library(SCINA)
library(openxlsx)


## ----message=FALSE, warning=FALSE, eval = FALSE-------------------------------------------------------------------------
matrix_list <- readRDS("./01_scripts/saved_objects/matrix_list.rds")
days <- names(matrix_list)


# Create Seurat Objects

for (i in seq_along(matrix_list)) {
	raw_object <- CreateSeuratObject(counts = matrix_list[[i]], project = days[i], min.cells = 5, min.features = 200)
	assign(paste0(days[i], "_raw"), raw_object)  # els Seurat objects estan guardats com a J0_raw, etc. Important options like min.cells 5 and min.features as 200 must be considered before starting the analysis                               
}


# QC ---------------------------------------------------------------------------
J0_raw$percent.mt <- PercentageFeatureSet(J0_raw, pattern = "^MT-") # Here, mitochondrial Genes are being identified
J0_raw$percent.rb <- PercentageFeatureSet(J0_raw, pattern = "^RP[SL]") # Ribosomes Genes are being identified
J0_counts <- J0_raw@assays$RNA$counts
tsb<- which(grepl("-TsB", rownames(J0_counts))) # -TsB marker genes must be deleted
cmo3 <- which(grepl("CMO3", rownames(J0_counts))) # CM03 marker genes must be deleted
J0_counts <- J0_counts[-c(cmo3, tsb),]
J0_raw <- subset(J0_raw, features = rownames(J0_counts))


J2_raw$percent.mt <- PercentageFeatureSet(J2_raw, pattern = "^MT-")
J2_raw$percent.rb <- PercentageFeatureSet(J2_raw, pattern = "^RP[SL]")
J2_counts <- J2_raw@assays$RNA$counts
tsb<- which(grepl("-TsB", rownames(J2_counts)))
cmo3 <- which(grepl("CMO3", rownames(J2_counts)))
J2_counts <- J2_counts[-c(cmo3, tsb),]
J2_raw <- subset(J2_raw, features = rownames(J2_counts))


J7_raw$percent.mt <- PercentageFeatureSet(J7_raw, pattern = "^MT-")
J7_raw$percent.rb <- PercentageFeatureSet(J7_raw, pattern = "^RP[SL]")
J7_counts <- J7_raw@assays$RNA$counts
tsb<- which(grepl("-TsB", rownames(J7_counts)))
cmo3 <- which(grepl("CMO3", rownames(J7_counts)))
J7_counts <- J7_counts[-c(cmo3, tsb),]
J7_raw <- subset(J7_raw, features = rownames(J7_counts))


J14_raw$percent.mt <- PercentageFeatureSet(J14_raw, pattern = "^MT-")
J14_raw$percent.rb <- PercentageFeatureSet(J14_raw, pattern = "^RP[SL]")
J14_counts <- J14_raw@assays$RNA$counts
tsb<- which(grepl("-TsB", rownames(J14_counts)))
cmo3 <- which(grepl("CMO3", rownames(J14_counts)))
J14_counts <- J14_counts[-c(cmo3, tsb),]
J14_raw <- subset(J14_raw, features = rownames(J14_counts))


J21_raw$percent.mt <- PercentageFeatureSet(J21_raw, pattern = "^MT-")
J21_raw$percent.rb <- PercentageFeatureSet(J21_raw, pattern = "^RP[SL]")
J21_counts <- J21_raw@assays$RNA$counts
tsb<- which(grepl("-TsB", rownames(J21_counts)))
cmo3 <- which(grepl("CMO3", rownames(J21_counts)))
J21_counts <- J21_counts[-c(cmo3, tsb),]
J21_raw <- subset(J21_raw, features = rownames(J21_counts))

# QC
head(J0_raw@meta.data, 5)
# violin plots control
VlnPlot(J0_raw, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"),
				ncol = 2)

## scatterplots ------------------------------------------------------------
dir.create("./01_scripts/output/plots/QC")
objects <- lst(J0_raw, J2_raw, J7_raw, J14_raw, J21_raw) # this lst is from purrr package instead of basic list(), to keep the names object names in the list
output_path <- "./01_scripts/output/plots/QC/"

for (i in seq_along(objects)) {
	# Create plots
	plot1 <- FeatureScatter(objects[[i]], feature1 = "nCount_RNA", feature2 = "percent.mt")
	plot2 <- FeatureScatter(objects[[i]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
	plot3 <- FeatureScatter(objects[[i]], feature1 = "nCount_RNA", feature2 = "percent.rb")
	plot4 <- FeatureScatter(objects[[i]], feature1 = "percent.rb", feature2 = "percent.mt")
	
	# Create a PDF file to save the plots
	pdf_file <- paste0(output_path, "QCplots_", names(objects)[i], ".pdf")
	pdf(file = pdf_file)
	
	# Save plots to the PDF file
	plot(plot1 + plot2 + plot3 + plot4)
	
	# Close the PDF file
	dev.off()
} 

## Own plots ##-------------------------------------------------------------
as_tibble(
	J0_raw[[]],
	rownames="Cell.Barcode") -> J0.qc.metrics

head(J0.qc.metrics)

J0.qc.metrics %>%
	arrange(percent.mt) %>%
	ggplot(aes(nCount_RNA,nFeature_RNA,colour=percent.mt)) + 
	geom_point(size=0.7) + 
	scale_color_gradientn(colors=c("black","blue","green2","red","yellow")) +
	ggtitle("All QC metrics")
## -------------------------------------------------------------------------

# ------ sub-setting after looking the plots ---------------------------------##

J0_sub <- subset(J0_raw, subset = nFeature_RNA > 4000 & 
								 	nFeature_RNA < 10000 & 
								 	percent.mt < 7 &
								 	nCount_RNA < 150000 &
								 	percent.rb > 10 & 
								 	percent.rb < 30)

J2_sub <- subset(J2_raw, subset = nFeature_RNA > 3000 &
								 	nFeature_RNA < 10500 & 
								 	percent.mt < 12 &
								 	percent.rb > 10 & 
								 	percent.rb < 30)

J7_sub <- subset(J7_raw, subset = nFeature_RNA > 3000 &
								 	nFeature_RNA < 10500 & 
								 	percent.mt < 12 &
								 	percent.rb > 8 & 
								 	percent.rb < 30)

J14_sub <- subset(J14_raw, subset = nFeature_RNA > 2000 &
										nFeature_RNA < 9500 & 
										percent.mt < 10 & 
										percent.rb > 8 & 
										percent.rb < 35)

J21_sub <- subset(J21_raw, subset = nFeature_RNA > 1500 &
										nFeature_RNA < 10000 & 
										percent.mt < 12 & 
										percent.rb > 8 & 
										percent.rb < 35)

#FeatureScatter(J21_raw, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

saveRDS(J0_sub, "./01_scripts/saved_objects/J0_sub.rds")
saveRDS(J2_sub, "./01_scripts/saved_objects/J2_sub.rds")
saveRDS(J7_sub, "./01_scripts/saved_objects/J7_sub.rds")
saveRDS(J14_sub, "./01_scripts/saved_objects/J14_sub.rds")
saveRDS(J21_sub, "./01_scripts/saved_objects/J21_sub.rds")

### Normalization, Variable regions, PCA ---------------------------------------

J0_sub <- readRDS("./01_scripts/saved_objects/J0_sub.rds")
J2_sub <- readRDS("./01_scripts/saved_objects/J2_sub.rds")
J7_sub <- readRDS("./01_scripts/saved_objects/J7_sub.rds")
J14_sub <- readRDS("./01_scripts/saved_objects/J14_sub.rds")
J21_sub <- readRDS("./01_scripts/saved_objects/J21_sub.rds")

sub <- lst(J0_sub, J2_sub, J7_sub, J14_sub, J21_sub)

for (i in seq_along(sub)) {
	# Normalize data
	sub[[i]] <- NormalizeData(sub[[i]])
	# FindVariableFeatures
	sub[[i]] <- FindVariableFeatures(sub[[i]], nfeatures = 5000)  # nfeatures will be 5000 for all
	assign(paste0("J", c(0, 2, 7, 14, 21)[i]), sub[[i]]) # Rename the object( (like J0, J2 and so on)
}


# plot top15 variablefeatures per hour ---------------------------------------
hours <- lst(J0, J2, J7, J14, J21)

pdf("./01_scripts/output/plots/VariableFeaturePlots_ok.pdf")
for (i in seq_along(hours)) {
	# Get the top 15 variable features
	top15 <- head(VariableFeatures(hours[[i]]), 15)
	# Create the variable feature plot
	plot <- VariableFeaturePlot(hours[[i]]) + ggtitle(names(hours[i]))
	# Label the top 15 features
	plot <- LabelPoints(plot = plot, points = top15, repel = FALSE)
	print(plot)
}
dev.off()

# First scale, without regressing anything --------------------------- 
data_list <- list(J0, J2, J7, J14, J21)
names_list <- c("J0", "J2", "J7", "J14", "J21")

# Iterate over the list
for (i in seq_along(data_list)) {
	# Scale the data
	data_list[[i]] <- ScaleData(data_list[[i]])
	
	# Run PCA
	pca <- RunPCA(data_list[[i]])
	
	# Visualize VizDimLoadings
	VizDimLoadings(pca, dims = 1:2, reduction = "pca")
}
# --------------------------------------------------------------------

# Incorporating cell cycle genes to regress, but regressing out just the difference between
# G2M and S phases (profileration). Cycling (G2M,S) vs non-cycling (G1) will be maintained
s.genes <- cc.genes$s.genes # cc.genes are included in Seurat packages
g2m.genes <- cc.genes$g2m.genes

# First, assign cell-cycle scores to all J's
hours <- lst(J0, J2, J7, J14, J21)
add_scores <- lapply(names(hours), function(name) {
	updated <- CellCycleScoring(hours[[name]], s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
})
add_scores <- setNames(add_scores, paste0(names(hours), "_scored"))
# new objects as separate variables in the global environment (J0_scored, etc)
list2env(add_scores, envir = .GlobalEnv)


# re-scaling regressing by cell cycle difference (S/G2M vs G1), and MT gene expression -----
rescale <- function(data) {
	# Calculate CC.Dif
	data$CC.Dif <- data$S.Score - data$G2M.Score
	rescaled <- ScaleData(data, features = rownames(data), vars.to.regress = c("CC.Dif", "percent.mt"))
	return(rescaled)
}

rescaled_hours <- lapply(list(J0_scored, J2_scored, J7_scored, J14_scored, J21_scored), rescale)

J0_rescaled <- rescaled_hours[[1]]
J2_rescaled <- rescaled_hours[[2]]
J7_rescaled <- rescaled_hours[[3]]
J14_rescaled <- rescaled_hours[[4]]
J21_rescaled <- rescaled_hours[[5]]

# Save fianlly scaled objects
saveRDS(J0_rescaled, "./01_scripts/saved_objects/J0_rescaled.rds")
saveRDS(J2_rescaled, "./01_scripts/saved_objects/J2_rescaled.rds")
saveRDS(J7_rescaled, "./01_scripts/saved_objects/J7_rescaled.rds")
saveRDS(J14_rescaled, "./01_scripts/saved_objects/J14_rescaled.rds")
saveRDS(J21_rescaled, "./01_scripts/saved_objects/J21_rescaled.rds")


# RunPCA in the rescaled objects -----------------------------------------------
J0_rescaled <- readRDS("./01_scripts/saved_objects/J0_rescaled.rds")
J2_rescaled <- readRDS("./01_scripts/saved_objects/J2_rescaled.rds")
J7_rescaled <- readRDS("./01_scripts/saved_objects/J7_rescaled.rds")
J14_rescaled <- readRDS("./01_scripts/saved_objects/J14_rescaled.rds")
J21_rescaled <- readRDS("./01_scripts/saved_objects/J21_rescaled.rds")


J0.pca.res <- RunPCA(J0_rescaled, features=VariableFeatures(J0_rescaled))
print(J0.pca.res[["pca"]], dims = 1:5, nfeatures = 5) 
VizDimLoadings(J0.pca.res, dims = 1:2, reduction = "pca")

J2.pca.res <- RunPCA(J2_rescaled, features=VariableFeatures(J2_rescaled))
print(J2.pca.res[["pca"]], dims = 1:5, nfeatures = 5) 
VizDimLoadings(J2.pca.res, dims = 1:2, reduction = "pca")    

J7.pca.res <- RunPCA(J7_rescaled, features=VariableFeatures(J7_rescaled))
print(J7.pca.res[["pca"]], dims = 1:5, nfeatures = 5,) 
VizDimLoadings(J7.pca.res, dims = 1:2, reduction = "pca")    

J14.pca.res <- RunPCA(J14_rescaled,features=VariableFeatures(J14_rescaled))
print(J14.pca.res[["pca"]], dims = 1:5, nfeatures = 5) 
VizDimLoadings(J14.pca.res, dims = 1:2, reduction = "pca")    

J21.pca.res <- RunPCA(J21_rescaled,features=VariableFeatures(J21_rescaled))
print(J21.pca.res[["pca"]], dims = 1:5, nfeatures = 5) 
VizDimLoadings(J21.pca.res, dims = 1:2, reduction = "pca")    

# Plots PCA cell cycle -----------------------------------------------------
data_list <- list(J0.pca.res, J2.pca.res, J7.pca.res, J14.pca.res, J21.pca.res)
names_list <- c("J0","J2","J7","J14","J21")


# Iterate over the list
for (i in seq_along(data_list)) {
	# PCA plot is created with the classification of each cell cycle
	pca_plot <- DimPlot(data_list[[i]], reduction = "pca", cols = c("S" = "indianred3", "G1" = "steelblue4", "G2M" = "orange")) +
		labs(title = paste(names_list[i], "- PCA with cell cycle classification"))
	
	# We save the plot as image
	png(paste0("./01_scripts/output/plots/", names_list[i], "_PCA_cell_cycle.png"),
			width = 1600, height = 900, res = 300)
	print(pca_plot)
	dev.off()
}


## fi normalization, scale and PCA -------------------------------------------##  

#### Dimensionality Reduction ------------------------------------------------##

# Determine the ‘dimensionality’ for each hour -------------------------------                    

# Iterate over the list
for (i in seq_along(data_list)) {
	# Generate Elbow plot 
	ElbowPlot(data_list[[i]])
	
	# Generate DimHeatmap
	DimHeatmap(data_list[[i]], dims = 1:21, cells = 500)
}

## Defining cell clusters ------------------------------------------------------

find_clusters <- function(pca.res) {
	# Find neighbors
	cl <- FindNeighbors(pca.res, dims = 1:20)
	# Find clusters
	cl <- FindClusters(cl, resolution = 0.5)
	return(cl)
}

cl_objects <- lapply(list(J0.pca.res, J2.pca.res, J7.pca.res, J14.pca.res, J21.pca.res), find_clusters)

J0.cl <- cl_objects[[1]]
J2.cl <- cl_objects[[2]]
J7.cl <- cl_objects[[3]]
J14.cl <- cl_objects[[4]]
J21.cl <- cl_objects[[5]]


## Run UMAP --------------------------------------------------------------------
saved.seed <- 8293
set.seed(saved.seed)
#-------------------------------------------------------------------------------

cluster_objects <- list(J0.cl, J2.cl, J7.cl, J14.cl, J21.cl)
for (i in seq_along(cluster_objects)) {
	umap_result <- RunUMAP(cluster_objects[[i]], dims = 1:15,
												 seed.use = saved.seed)
	cluster_objects[[i]] <- umap_result
	assign(paste0("J", c(0, 2, 7, 14, 21)[i], ".cl"), umap_result, envir = .GlobalEnv)
}

saveRDS(J0.cl, "./01_scripts/saved_objects/J0_umap.rds")
saveRDS(J2.cl, "./01_scripts/saved_objects/J2_umap.rds")
saveRDS(J7.cl, "./01_scripts/saved_objects/J7_umap.rds")
saveRDS(J14.cl, "./01_scripts/saved_objects/J14_umap.rds")
saveRDS(J21.cl, "./01_scripts/saved_objects/J21_umap.rds")
# ---------------------------------------------------#
J0.cl <- readRDS("./01_scripts/saved_objects/J0_umap.rds")
J2.cl <- readRDS("./01_scripts/saved_objects/J2_umap.rds")
J7.cl <- readRDS("./01_scripts/saved_objects/J7_umap.rds")
J14.cl <- readRDS("./01_scripts/saved_objects/J14_umap.rds")
J21.cl <- readRDS("./01_scripts/saved_objects/J21_umap.rds")


# UMAP plots with seurat clusters --------------------------------------------
data_list <- list(J0.cl,J2.cl,J7.cl,J14.cl,J21.cl)
names_list <- c("J0","J2","J7","J14","J21")

# Iterate over the list
for (i in seq_along(data_list)) {
	# UMAP plot is created
	umap_plot <- DimPlot(data_list[[i]], reduction = "umap") + 
		labs(title = paste(names_list[i], "- cell clusters"))
	
	# Save plot as image
	png(paste0("./01_scripts/output/plots/", names_list[i], "_plots/", names_list[i], "_UMAP_seurat-clusters.png"), 
			width = 1600, height = 1300, res = 300)
	print(umap_plot)
	dev.off()
}


## Examining the properties (QC metrics) of the clusters) ----------------------

# basic QC -------------------------------------------------------------------
Idents(J21.cl) <- "seurat_clusters"  
plot1 <- VlnPlot(J21.cl,features="nCount_RNA") # number of reads
plot2 <- VlnPlot(J21.cl,features="nFeature_RNA") # number of genes
plot3 <- FeaturePlot(J21.cl, features = "percent.mt") # percent.mt vs nCount_RNA
plot4 <- FeaturePlot(J21.cl, features = "percent.rb") # percent.rb vs nCount_RNA
plot1 + plot2 + plot3 + plot4
# ----------------------------------------------------------------------------

# Cell cycle issues ----------------------------------------------------------
cols <- c("S" = "indianred3",
				 "G1" = "steelblue4",
				 "G2M" = "orange")

# Iterate over the list
for (i in seq_along(data_list)) {
	# Assign the labels of each Cell Cycle
	Idents(data_list[[i]]) <- "Phase"
	
	# Create the Cell Cycle UMAP with cluster by each cell cycle
	um_plot <- DimPlot(data_list[[i]], reduction = "umap", cols = cols) +
		labs(title = paste(names_list[i], "- UMAP with cell cycle classification"))
	
	# Save graphic as image for each timepoint
	png(paste0("./01_scripts/output/plots/", names_list[i], "_plots/", names_list[i], "_UMAP_cellcycle.png"), 
			width = 1600, height = 1300, res = 300)
	print(um_plot)
	dev.off()
}

#### Differentially expressed features (cluster biomarkers) ------------------##

Idents(J0.cl) <- "seurat_clusters"
Idents(J2.cl) <- "seurat_clusters"
Idents(J7.cl) <- "seurat_clusters"
Idents(J14.cl) <- "seurat_clusters"

write_to_excel <- function(data, filename, sheetName) {
	if (!file.exists(filename)) {
		write.xlsx(data, filename, sheetName = sheetName, rowNames = FALSE)
	} else {
		wb <- loadWorkbook(filename)
		addWorksheet(wb, sheetName)
		writeData(wb, sheet = sheetName, x = data)
		saveWorkbook(wb, filename, overwrite = TRUE)
	}
}

J0.markers <- FindAllMarkers(J0.cl,
														 logfc.threshold = 0.15,
														 min.pct = 0.20,
														 only.pos = TRUE)
sig.J0.markers <- J0.markers %>%
	group_by(cluster) %>%
	filter(avg_log2FC > 0.99 & p_val_adj < 0.05) %>%
	arrange(-avg_log2FC, .by_group = TRUE)
write_to_excel(as.data.frame(sig.J0.markers), "./01_scripts/output/tables/sig.markers.by.timepoint.xlsx", sheetName = "J0")

J2.markers <- FindAllMarkers(J2.cl,
														 logfc.threshold = 0.15,
														 min.pct = 0.20,
														 only.pos = TRUE)
sig.J2.markers <- J2.markers %>%
	group_by(cluster) %>%
	filter(avg_log2FC > 0.99 & p_val_adj < 0.05) %>%
	arrange(-avg_log2FC, .by_group = TRUE)
write_to_excel(as.data.frame(sig.J2.markers), "./01_scripts/output/tables/sig.markers.by.timepoint.xlsx", sheetName = "J2")

J7.markers <- FindAllMarkers(J7.cl,
														 logfc.threshold = 0.15,
														 min.pct = 0.20,
														 only.pos = TRUE)
sig.J7.markers <- J7.markers %>%
	group_by(cluster) %>%
	filter(avg_log2FC > 0.99 & p_val_adj < 0.05) %>%
	arrange(-avg_log2FC, .by_group = TRUE)
write_to_excel(as.data.frame(sig.J7.markers), "./01_scripts/output/tables/sig.markers.by.timepoint.xlsx", sheetName = "J7")

J14.markers <- FindAllMarkers(J14.cl,
														 logfc.threshold = 0.15,
														 min.pct = 0.20,
														 only.pos = TRUE)
sig.J14.markers <- J14.markers %>%
	group_by(cluster) %>%
	filter(avg_log2FC > 0.99 & p_val_adj < 0.05) %>%
	arrange(-avg_log2FC, .by_group = TRUE)
write_to_excel(as.data.frame(sig.J14.markers), "./01_scripts/output/tables/sig.markers.by.timepoint.xlsx", sheetName = "J14")

# Visual control: plot markers per cluster -----------------------------
j14.best.gene.cluster <- sig.J14.markers %>%
	group_by(cluster) %>%
	slice(1) %>%
	pull(gene)
VlnPlot(J14.cl,features=j14.best.gene.cluster)
# ----------------------------------------------------------------------


#### Automated Cell Type Annotation with SCINA -------------------------------##

## signatures from CancerSEA (14 functional states) ----------------------------

# import gene signatures -----------------------------------------------------      
file <- "03_other_files/14signatures.xlsx"
signatures <- import_list(file)
gene.signatures <- lapply(signatures, function(df) df$GeneName)

J0.data <- J0.cl@assays$RNA$scale.data 

J0.scina.results <- SCINA(J0.data,
													gene.signatures, 
													max_iter = 100, 
													convergence_n = 10, 
													convergence_rate = 0.999, 
													sensitivity_cutoff = 0.9, 
													rm_overlap=FALSE, 
													allow_unknown=TRUE) 

J0.cl$scina_labels <- J0.scina.results$cell_labels


J2.data <- J2.cl@assays$RNA$scale.data 

J2.scina.results <- SCINA(J2.data,
													gene.signatures, 
													max_iter = 100, 
													convergence_n = 10, 
													convergence_rate = 0.999, 
													sensitivity_cutoff = 0.9, 
													rm_overlap=FALSE, 
													allow_unknown=TRUE) 

J2.cl$scina_labels <- J2.scina.results$cell_labels


J7.data <- J7.cl@assays$RNA$scale.data 

J7.scina.results <- SCINA(J7.data,
													gene.signatures, 
													max_iter = 100, 
													convergence_n = 10, 
													convergence_rate = 0.999, 
													sensitivity_cutoff = 0.9, 
													rm_overlap=FALSE, 
													allow_unknown=TRUE) 

J7.cl$scina_labels <- J7.scina.results$cell_labels


J14.data <- J14.cl@assays$RNA$scale.data 

J14.scina.results <- SCINA(J14.data,
													 gene.signatures, 
													 max_iter = 100, 
													 convergence_n = 10, 
													 convergence_rate = 0.999, 
													 sensitivity_cutoff = 0.9, 
													 rm_overlap=FALSE, 
													 allow_unknown=TRUE) 

J14.cl$scina_labels <- J14.scina.results$cell_labels
# save final objects ---------------------------------------------------------
saveRDS(J0.cl, "./01_scripts/saved_objects/J0_scina_annot.rds")
saveRDS(J2.cl, "./01_scripts/saved_objects/J2_scina_annot.rds")
saveRDS(J7.cl, "./01_scripts/saved_objects/J7_scina_annot.rds")
saveRDS(J14.cl, "./01_scripts/saved_objects/J14_scina_annot.rds")

J0.cl <- readRDS("./01_scripts/saved_objects/J0_scina_annot.rds")
J2.cl <- readRDS("./01_scripts/saved_objects/J2_scina_annot.rds")
J7.cl <- readRDS("./01_scripts/saved_objects/J7_scina_annot.rds")
J14.cl <- readRDS("./01_scripts/saved_objects/J14_scina_annot.rds")

# save the metadata relating seurat clusters and scina_labels per jour -------
metadata.0 <- J0.cl@meta.data[, c("seurat_clusters", "scina_labels")]
label_cluster.0 <- as.data.frame(table(metadata.0$seurat_clusters, metadata.0$scina_labels)) %>%
	arrange(Var1)
colnames(label_cluster.0) <- c("Seurat_cluster", "Predicted_signature", "Frequency")

metadata.2 <- J2.cl@meta.data[, c("seurat_clusters", "scina_labels")]
label_cluster.2 <- as.data.frame(table(metadata.2$seurat_clusters, metadata.2$scina_labels)) %>%
	arrange(Var1)
colnames(label_cluster.2) <- c("Seurat_cluster", "Predicted_signature", "Frequency")

metadata.7 <- J7.cl@meta.data[, c("seurat_clusters", "scina_labels")]
label_cluster.7 <- as.data.frame(table(metadata.7$seurat_clusters, metadata.7$scina_labels)) %>%
	arrange(Var1)
colnames(label_cluster.7) <- c("Seurat_cluster", "Predicted_signature", "Frequency")

metadata.14 <- J14.cl@meta.data[, c("seurat_clusters", "scina_labels")]
label_cluster.14 <- as.data.frame(table(metadata.14$seurat_clusters, metadata.14$scina_labels)) %>%
	arrange(Var1)
colnames(label_cluster.14) <- c("Seurat_cluster", "Predicted_signature", "Frequency")


# add dge markers identified with seurat clusters ----------------------------
# open markers file  
sig.markers <- "./01_scripts/output/tables/sig.markers.by.timepoint.xlsx"
data <- import_list(sig.markers)
list2env(data,envir=.GlobalEnv) # creates individual dataframes for each sheet

genes.clusters.0 <- J0 %>%
	group_by(cluster) %>%
	summarise(sig.seurat.markers = paste(gene, collapse = ","))
genes.clusters.2 <- J2 %>%
	group_by(cluster) %>%
	summarise(sig.seurat.markers = paste(gene, collapse = ","))
genes.clusters.7 <- J7 %>%
	group_by(cluster) %>%
	summarise(sig.seurat.markers = paste(gene, collapse = ","))
genes.clusters.14 <- J14 %>%
	group_by(cluster) %>%
	summarise(sig.seurat.markers = paste(gene, collapse = ","))


label_cluster.0 <- merge(label_cluster.0, genes.clusters.0, by.x="Seurat_cluster", by.y="cluster", all.x = TRUE)
label_cluster.2 <- merge(label_cluster.2, genes.clusters.2, by.x="Seurat_cluster", by.y="cluster", all.x = TRUE)
label_cluster.7 <- merge(label_cluster.7, genes.clusters.7, by.x="Seurat_cluster", by.y="cluster", all.x = TRUE)
label_cluster.14 <- merge(label_cluster.14, genes.clusters.14, by.x="Seurat_cluster", by.y="cluster", all.x = TRUE)


write.xlsx(label_cluster.0, "./01_scripts/output/tables/seurat_cluster_scina_labels.by.timepoint.xlsx", sheetName= "J0", row.names = F)
write.xlsx(label_cluster.2, "./01_scripts/output/tables/seurat_cluster_scina_labels.by.timepoint.xlsx", sheetName= "J2", row.names = F, append = T)
write.xlsx(label_cluster.7, "./01_scripts/output/tables/seurat_cluster_scina_labels.by.timepoint.xlsx", sheetName= "J7", row.names = F, append = T)
write.xlsx(label_cluster.14, "./01_scripts/output/tables/seurat_cluster_scina_labels.by.timepoint.xlsx", sheetName= "J14", row.names = F, append = T)

# Define the colors of each cluster
cols <- c("Stemness" = 'darkorange3', "Quiescence" = "blue4", "Proliferation" = "orange1",
					"Metastasis" = "turquoise", "Invasion" = "plum", "Inflammation" = "khaki2", "Hypoxia" = "green", 
					"EMT" = "darkgreen", "DNArepair" = "mediumvioletred", "DNAdamage" = "palevioletred4", 
					"Differentiation" = "antiquewhite3", "Cell.cycle" = "slateblue1", "Apoptosis" = "gray41", "Angiogenesis" = "red") 


Idents(J0.cl) <- "scina_labels"
J0.annot <- DimPlot(J0.cl, reduction = "umap", cols = cols) + 
	labs(title = "J0 - 14 cancer signatures")
png("./01_scripts/output/plots/J0_UMAP_scina-14sig.png", width=1750,height=1300, res=300)
J0.annot
dev.off()


Idents(J2.cl) <- "scina_labels"
J2.annot <- DimPlot(J2.cl, reduction = "umap", cols = cols) + 
	labs(title = "J2 - 14 cancer signatures")
png("./01_scripts/output/plots/J2_UMAP_scina-14sig.png", width=1750,height=1300, res=300)
J2.annot
dev.off()


Idents(J7.cl) <- "scina_labels"
J7.annot <- DimPlot(J7.cl, reduction = "umap", cols = cols) + 
	labs(title = "J7 - 14 cancer signatures")
png("./01_scripts/output/plots/J7_UMAP_scina-14sig.png", width=1750,height=1300, res=300)
J7.annot
dev.off()


Idents(J14.cl) <- "scina_labels"
J14.annot <- DimPlot(J14.cl, reduction = "umap", cols = cols) + 
	labs(title = "J14 - 14 cancer signatures")
png("./01_scripts/output/plots/J14_UMAP_scina-14sig.png", width=1750,height=1300, res=300)
J14.annot
dev.off()


# Proportion plot -------------------------------------------------

labels <- names(gene.signatures)

label_counts.j0 <- table(as.data.frame(J0.cl@meta.data$scina_labels))
label_proportions.j0 <- prop.table(label_counts.j0) * 100

label_counts.j2 <- table(as.data.frame(J2.cl@meta.data$scina_labels))
label_proportions.j2 <- prop.table(label_counts.j2) * 100

label_counts.j7 <- table(as.data.frame(J7.cl@meta.data$scina_labels))
label_proportions.j7 <- prop.table(label_counts.j7) * 100

label_counts.j14 <- table(as.data.frame(J14.cl@meta.data$scina_labels))
label_proportions.j14 <- prop.table(label_counts.j14) * 100

plot_data.j0 <- data.frame(label = names(label_proportions.j0), proportion = as.numeric(label_proportions.j0))
colnames(plot_data.j0) <- c("label", "J0")
plot_data.j2 <- data.frame(label = names(label_proportions.j2), proportion = as.numeric(label_proportions.j2))
colnames(plot_data.j2) <- c("label", "J2")
plot_data.j7 <- data.frame(label = names(label_proportions.j7), proportion = as.numeric(label_proportions.j7))
colnames(plot_data.j7) <- c("label", "J7")
plot_data.j14 <- data.frame(label = names(label_proportions.j14), proportion = as.numeric(label_proportions.j14))
colnames(plot_data.j14) <- c("label", "J14")

library(str2str)
all.plot_data <- cbind_fill(plot_data.j0, plot_data.j2, plot_data.j7, plot_data.j14)
all.plot_data <- all.plot_data[,c(1,2,4,6,8)]

plot_data_long <- tidyr::pivot_longer(all.plot_data, cols = -label, names_to = "Sample", values_to = "Proportion")
plot_data_long$Sample <- factor(plot_data_long$Sample, levels = unique(plot_data_long$Sample))

library(cowplot)
prop_pl <- ggplot(plot_data_long, aes(x = Sample, y = Proportion, fill = label)) +
	geom_bar(stat="identity", alpha = 0.6) + xlab("") + ylab("Percentage") +
	theme_cowplot() + 
	theme(axis.text.y = element_text(size=11),
				axis.text.x = element_text(size=11),
				axis.title.x = element_text(size=12),
				axis.title.y = element_text(size=12),
				legend.title = element_blank(),
				legend.text = element_text(size=12)) +
	scale_fill_manual(values=cols)
prop_pl

png("./01_scripts/output/plots/Proprortion_Plot_by-time-point_scina_labels.png",width=1750,height=1300, res=300)
prop_pl
dev.off()

##------------------------

label_counts.j0_2 <- table(as.data.frame(J0.cl@meta.data$seurat_clusters))
label_proportions.j0_2 <- prop.table(label_counts.j0_2) * 100

label_counts.j2_2 <- table(as.data.frame(J2.cl@meta.data$seurat_clusters))
label_proportions.j2_2 <- prop.table(label_counts.j2_2) * 100

label_counts.j7_2 <- table(as.data.frame(J7.cl@meta.data$seurat_clusters))
label_proportions.j7_2 <- prop.table(label_counts.j7_2) * 100

label_counts.j14_2 <- table(as.data.frame(J14.cl@meta.data$seurat_clusters))
label_proportions.j14_2 <- prop.table(label_counts.j14_2) * 100

plot_data.j0 <- data.frame(label = names(label_proportions.j0_2), proportion = as.numeric(label_proportions.j0_2))
colnames(plot_data.j0) <- c("label", "J0")
plot_data.j2 <- data.frame(label = names(label_proportions.j2_2), proportion = as.numeric(label_proportions.j2_2))
colnames(plot_data.j2) <- c("label", "J2")
plot_data.j7 <- data.frame(label = names(label_proportions.j7_2), proportion = as.numeric(label_proportions.j7_2))
colnames(plot_data.j7) <- c("label", "J7")
plot_data.j14 <- data.frame(label = names(label_proportions.j14_2), proportion = as.numeric(label_proportions.j14_2))
colnames(plot_data.j14) <- c("label", "J14")

library(str2str)
all.plot_data <- cbind_fill(plot_data.j0, plot_data.j2, plot_data.j7, plot_data.j14)
all.plot_data <- all.plot_data[,c(1,2,4,6,8)]

plot_data_long <- tidyr::pivot_longer(all.plot_data, cols = -label, names_to = "Sample", values_to = "Proportion")
plot_data_long$Sample <- factor(plot_data_long$Sample, levels = unique(plot_data_long$Sample))

plot_data_long$label <- as.factor(plot_data_long$label)

library(cowplot)
prop_pl <- ggplot(plot_data_long, aes(x = Sample, y = Proportion, fill = label)) +
	geom_bar(stat="identity", alpha = 0.6) + xlab("") + ylab("Percentage") +
	theme_cowplot() + 
	theme(axis.text.y = element_text(size=11),
				axis.text.x = element_text(size=11),
				axis.title.x = element_text(size=12),
				axis.title.y = element_text(size=12),
				legend.title = element_blank(),
				legend.text = element_text(size=12))
prop_pl

png("./01_scripts/output/plots/Proprortion_Plot_by-time-point_seurat_clusters.png",width=1750,height=1300, res=300)
prop_pl
dev.off()

## Proportion plot phase ------------------------

cols <- c("S" = "indianred3",
					"G1" = "steelblue4",
					"G2M" = "orange")

label_counts.j0_p <- table(as.data.frame(J0.cl@meta.data$Phase))
label_proportions.j0_p <- prop.table(label_counts.j0_p) * 100

label_counts.j2_p <- table(as.data.frame(J2.cl@meta.data$Phase))
label_proportions.j2_p <- prop.table(label_counts.j2_p) * 100

label_counts.j7_p <- table(as.data.frame(J7.cl@meta.data$Phase))
label_proportions.j7_p <- prop.table(label_counts.j7_p) * 100

label_counts.j14_p <- table(as.data.frame(J14.cl@meta.data$Phase))
label_proportions.j14_p <- prop.table(label_counts.j14_p) * 100

plot_data.j0 <- data.frame(label = names(label_proportions.j0_p), proportion = as.numeric(label_proportions.j0_p))
colnames(plot_data.j0) <- c("label", "J0")
plot_data.j2 <- data.frame(label = names(label_proportions.j2_p), proportion = as.numeric(label_proportions.j2_p))
colnames(plot_data.j2) <- c("label", "J2")
plot_data.j7 <- data.frame(label = names(label_proportions.j7_p), proportion = as.numeric(label_proportions.j7_p))
colnames(plot_data.j7) <- c("label", "J7")
plot_data.j14 <- data.frame(label = names(label_proportions.j14_p), proportion = as.numeric(label_proportions.j14_p))
colnames(plot_data.j14) <- c("label", "J14")

all.plot_data <- cbind_fill(plot_data.j0, plot_data.j2, plot_data.j7, plot_data.j14)
all.plot_data <- all.plot_data[,c(1,2,4,6,8)]

plot_data_long <- tidyr::pivot_longer(all.plot_data, cols = -label, names_to = "Sample", values_to = "Proportion")
plot_data_long$Sample <- factor(plot_data_long$Sample, levels = unique(plot_data_long$Sample))


library(cowplot)
prop_pl <- ggplot(plot_data_long, aes(x = Sample, y = Proportion, fill = label)) +
	geom_bar(stat="identity", alpha = 0.6) + xlab("") + ylab("Percentage") +
	theme_cowplot() + 
	theme(axis.text.y = element_text(size=11),
				axis.text.x = element_text(size=11),
				axis.title.x = element_text(size=12),
				axis.title.y = element_text(size=12),
				legend.title = element_blank(),
				legend.text = element_text(size=12))+
	scale_fill_manual(values=cols)
prop_pl

png("./01_scripts/output/plots/Proprortion_Plot_by-time-point_cellcycle.png",width=1750,height=1300, res=300)
prop_pl
dev.off()



## LIANA ANALYSIS -----------------------------------------

library(tidyverse)
library(magrittr)
library(liana)
library(OmnipathR)

consensus_omni <- select_resource("Consensus")[[1]] %>%
	glimpse

op_resource <- curated_ligand_receptor_interactions(
	cellphonedb = TRUE,
	cellinker = TRUE,
	talklr = TRUE,
	signalink = TRUE) 


liana_object_J0 <- liana_wrap(J0.cl,
													 resource = 'custom',
													 external_resource = op_resource,
													 idents_col = "seurat_clusters",
													 min_cells = 45,
													 supp_columns = c("ligand.expr", "receptor.expr", "ligand.stat", "receptor.stat", "ligand.pval", "receptor.pval", "ligand.FDR", "receptor.FDR"),
													 assay = "RNA",
													 expr_prop = 0.2)

liana_object_J0_scina <- liana_wrap(J0.cl,
															resource = 'custom',
															external_resource = op_resource,
															idents_col = "scina_labels",
															min_cells = 45,
															supp_columns = c("ligand.expr", "receptor.expr", "ligand.stat", "receptor.stat", "ligand.pval", "receptor.pval", "ligand.FDR", "receptor.FDR"),
															assay = "RNA",
															expr_prop = 0.2)

liana_object_J2 <- liana_wrap(J2.cl,
															resource = 'custom',
															external_resource = op_resource,
															idents_col = "seurat_clusters",
															min_cells = 45,
															supp_columns = c("ligand.expr", "receptor.expr", "ligand.stat", "receptor.stat", "ligand.pval", "receptor.pval", "ligand.FDR", "receptor.FDR"),
															assay = "RNA",
															expr_prop = 0.2)

liana_object_J2_scina <- liana_wrap(J2.cl,
																		resource = 'custom',
																		external_resource = op_resource,
																		idents_col = "scina_labels",
																		min_cells = 45,
																		supp_columns = c("ligand.expr", "receptor.expr", "ligand.stat", "receptor.stat", "ligand.pval", "receptor.pval", "ligand.FDR", "receptor.FDR"),
																		assay = "RNA",
																		expr_prop = 0.2)

liana_object_J7 <- liana_wrap(J7.cl,
															resource = 'custom',
															external_resource = op_resource,
															idents_col = "seurat_clusters",
															min_cells = 45,
															supp_columns = c("ligand.expr", "receptor.expr", "ligand.stat", "receptor.stat", "ligand.pval", "receptor.pval", "ligand.FDR", "receptor.FDR"),
															assay = "RNA",
															expr_prop = 0.2)

liana_object_J7_scina <- liana_wrap(J7.cl,
															resource = 'custom',
															external_resource = op_resource,
															idents_col = "scina_labels",
															min_cells = 45,
															supp_columns = c("ligand.expr", "receptor.expr", "ligand.stat", "receptor.stat", "ligand.pval", "receptor.pval", "ligand.FDR", "receptor.FDR"),
															assay = "RNA",
															expr_prop = 0.2)

liana_object_J14 <- liana_wrap(J14.cl,
															resource = 'custom',
															external_resource = op_resource,
															idents_col = "seurat_clusters",
															min_cells = 45,
															supp_columns = c("ligand.expr", "receptor.expr", "ligand.stat", "receptor.stat", "ligand.pval", "receptor.pval", "ligand.FDR", "receptor.FDR"),
															assay = "RNA",
															expr_prop = 0.2)

liana_object_J14_scina <- liana_wrap(J14.cl,
															 resource = 'custom',
															 external_resource = op_resource,
															 idents_col = "scina_labels",
															 min_cells = 45,
															 supp_columns = c("ligand.expr", "receptor.expr", "ligand.stat", "receptor.stat", "ligand.pval", "receptor.pval", "ligand.FDR", "receptor.FDR"),
															 assay = "RNA",
															 expr_prop = 0.2)

saveRDS(liana_object_J0, "./01_scripts/saved_objects/liana_object_J0.RDS")
saveRDS(liana_object_J2, "./01_scripts/saved_objects/liana_object_J2.RDS")
saveRDS(liana_object_J7, "./01_scripts/saved_objects/liana_object_J7.RDS")
saveRDS(liana_object_J14, "./01_scripts/saved_objects/liana_object_J14.RDS")

saveRDS(liana_object_J0_scina, "./01_scripts/saved_objects/liana_object_J0_scina.RDS")
saveRDS(liana_object_J2_scina, "./01_scripts/saved_objects/liana_object_J2_scina.RDS")
saveRDS(liana_object_J7_scina, "./01_scripts/saved_objects/liana_object_J7_scina.RDS")
saveRDS(liana_object_J14_scina, "./01_scripts/saved_objects/liana_object_J14_scina.RDS")

# identify interactions of interest
J0_int <- liana_object_J0 %>%
	liana_aggregate() %>%
	# only keep interactions concordant between methods
	filter(aggregate_rank <= 0.05) %>% 
	# keep top 20 interactions (regardless of cell type)
	distinct_at(c("source","target","ligand.complex","receptor.complex","aggregate_rank")) %>% 
	head(50)

# identify interactions of interest
J0_int_scina <- liana_object_J0_scina %>%
	liana_aggregate() %>%
	# only keep interactions concordant between methods
	filter(aggregate_rank <= 0.05) %>% 
	# keep top 20 interactions (regardless of cell type)
	distinct_at(c("source","target","ligand.complex","receptor.complex","aggregate_rank")) %>% 
	head(50)

# identify interactions of interest
J2_int <- liana_object_J2 %>%
	liana_aggregate() %>%
	# only keep interactions concordant between methods
	filter(aggregate_rank <= 0.05) %>% 
	# keep top 20 interactions (regardless of cell type)
	distinct_at(c("source","target","ligand.complex","receptor.complex","aggregate_rank")) %>% 
	head(50)

# identify interactions of interest
J2_int_scina <- liana_object_J2_scina %>%
	liana_aggregate() %>%
	# only keep interactions concordant between methods
	filter(aggregate_rank <= 0.05) %>% 
	# keep top 20 interactions (regardless of cell type)
	distinct_at(c("source","target","ligand.complex","receptor.complex","aggregate_rank")) %>% 
	head(50)

# identify interactions of interest
J7_int <- liana_object_J7 %>%
	liana_aggregate() %>%
	# only keep interactions concordant between methods
	filter(aggregate_rank <= 0.05) %>% 
	# keep top 20 interactions (regardless of cell type)
	distinct_at(c("source","target","ligand.complex","receptor.complex","aggregate_rank")) %>% 
	head(50)

# identify interactions of interest
J7_int_scina <- liana_object_J7_scina %>%
	liana_aggregate() %>%
	# only keep interactions concordant between methods
	filter(aggregate_rank <= 0.05) %>% 
	# keep top 20 interactions (regardless of cell type)
	distinct_at(c("source","target","ligand.complex","receptor.complex","aggregate_rank")) %>% 
	head(50)

# identify interactions of interest
J14_int <- liana_object_J14 %>%
	liana_aggregate() %>%
	# only keep interactions concordant between methods
	filter(aggregate_rank <= 0.05) %>% 
	# keep top 20 interactions (regardless of cell type)
	distinct_at(c("source","target","ligand.complex","receptor.complex","aggregate_rank")) %>% 
	head(50)

# identify interactions of interest
J14_int_scina <- liana_object_J14_scina %>%
	liana_aggregate() %>%
	# only keep interactions concordant between methods
	filter(aggregate_rank <= 0.05) %>% 
	# keep top 20 interactions (regardless of cell type)
	distinct_at(c("source","target","ligand.complex","receptor.complex","aggregate_rank")) %>% 
	head(50)

dir.create("./01_scripts/output/tables/cell-cell_interaction")

write.xlsx(J0_int, "./01_scripts/output/tables/cell-cell_interaction/J0_liana_interactions.xlsx")
write.xlsx(J2_int, "./01_scripts/output/tables/cell-cell_interaction/J2_liana_interactions.xlsx")
write.xlsx(J7_int, "./01_scripts/output/tables/cell-cell_interaction/J7_liana_interactions.xlsx")
write.xlsx(J14_int, "./01_scripts/output/tables/cell-cell_interaction/J14_liana_interactions.xlsx")

write.xlsx(J0_int_scina, "./01_scripts/output/tables/cell-cell_interaction/J0_liana_interactions_scina.xlsx")
write.xlsx(J2_int_scina, "./01_scripts/output/tables/cell-cell_interaction/J2_liana_interactions_scina.xlsx")
write.xlsx(J7_int_scina, "./01_scripts/output/tables/cell-cell_interaction/J7_liana_interactions_scina.xlsx")
write.xlsx(J14_int_scina, "./01_scripts/output/tables/cell-cell_interaction/J14_liana_interactions_scina.xlsx")

# SEURAT CLUSTERS CELL-CELL INTERACTION ------------------------

liana_object_J0<-readRDS("./01_scripts/saved_objects/liana_object_J0.RDS")
liana_object_J2<-readRDS("./01_scripts/saved_objects/liana_object_J2.RDS")
liana_object_J7<-readRDS("./01_scripts/saved_objects/liana_object_J7.RDS")
liana_object_J14<-readRDS("./01_scripts/saved_objects/liana_object_J14.RDS")

# identify interactions of interest J0 BASED ON TABLE
J0_liana_object_top<- liana_object_J0 %>%
	liana_aggregate() %>%
	# only keep interactions concordant between methods
	filter(aggregate_rank <= 0.05)

library(circlize)

chord_freq(J0_liana_object_top)
J0_heatmap_l_seurat<-heat_freq(J0_liana_object_top)

J0_l_dotplot<-liana_dotplot(J0_liana_object_top,
							ntop = 10,
							source_groups = c("3", "5","2","6"),
							target_groups = c("1","2","5","6"))

dir.create("./01_scripts/output/plots/cell-cell_interaction_plots/by_timepoints")

png("./01_scripts/output/plots/cell-cell_interaction_plots/by_timepoints/J0_seurat_clusters_liana.png",
		width=5000,height=3500, res=300)

J0_l_dotplot

dev.off()

png("./01_scripts/output/plots/cell-cell_interaction_plots/by_timepoints/J0_seurat_clusters_heatmap.png",
		width=5000,height=3500, res=300)

J0_heatmap_l_seurat

dev.off()


# identify interactions of interest J2 BASED ON TABLE
J2_liana_object_top<- liana_object_J2 %>%
	liana_aggregate() %>%
	# only keep interactions concordant between methods
	filter(aggregate_rank <= 0.05)

chord_freq(J2_liana_object_top)
J2_heatmap_l_seurat<-heat_freq(J2_liana_object_top)

J2_l_dotplot<-liana_dotplot(J2_liana_object_top,
							ntop = 10,
							source_groups = c("3", "6","4","1","2"),
							target_groups = c("0","3","6","5","1","4"))

png("./01_scripts/output/plots/cell-cell_interaction_plots/by_timepoints/J2_seurat_clusters_liana.png",
		width=5700,height=3500, res=300)

J2_l_dotplot

dev.off()

png("./01_scripts/output/plots/cell-cell_interaction_plots/by_timepoints/J2_seurat_clusters_heatmap.png",
		width=5000,height=3500, res=300)

J2_heatmap_l_seurat

dev.off()

# identify interactions of interest J7 BASED ON TABLE
J7_liana_object_top<- liana_object_J7 %>%
	liana_aggregate() %>%
	# only keep interactions concordant between methods
	filter(aggregate_rank <= 0.05)

chord_freq(J7_liana_object_top)
J7_heatmap_l_seurat<-heat_freq(J7_liana_object_top)

J7_l_dotplot<-liana_dotplot(J7_liana_object_top,
							ntop = 10,
							source_groups = c("5", "4","1","6","2","0"),
							target_groups = c("6","1","2","0","3"))

png("./01_scripts/output/plots/cell-cell_interaction_plots/by_timepoints/J7_seurat_clusters_liana.png",
		width=5700,height=3500, res=300)

J7_l_dotplot

dev.off()

png("./01_scripts/output/plots/cell-cell_interaction_plots/by_timepoints/J7_seurat_clusters_heatmap.png",
		width=5000,height=3500, res=300)

J7_heatmap_l_seurat

dev.off()

# identify interactions of interest J14 BASED ON TABLE
J14_liana_object_top<- liana_object_J14 %>%
	liana_aggregate() %>%
	# only keep interactions concordant between methods
	filter(aggregate_rank <= 0.05)

chord_freq(J14_liana_object_top)
J14_heatmap_l_seurat<-heat_freq(J14_liana_object_top)

J14_l_dotplot<-liana_dotplot(J14_liana_object_top,
							ntop = 10,
							source_groups = c("5"),
							target_groups = c("5","2","3","1"))

png("./01_scripts/output/plots/cell-cell_interaction_plots/by_timepoints/J14_seurat_clusters_liana.png",
		width=4000,height=3500, res=300)

J14_l_dotplot

dev.off()

png("./01_scripts/output/plots/cell-cell_interaction_plots/by_timepoints/J14_seurat_clusters_heatmap.png",
		width=5000,height=3500, res=300)

J14_heatmap_l_seurat

dev.off()


# SCINA CELL-CELL INTERACTION -------------------------

liana_object_J0_scina<-readRDS("./01_scripts/saved_objects/liana_object_J0_scina.RDS")
liana_object_J2_scina<-readRDS("./01_scripts/saved_objects/liana_object_J2_scina.RDS")
liana_object_J7_scina<-readRDS("./01_scripts/saved_objects/liana_object_J7_scina.RDS")
liana_object_J14_scina<-readRDS("./01_scripts/saved_objects/liana_object_J14_scina.RDS")

# identify interactions of interest J0 BASED ON TABLE
J0_liana_object_top<- liana_object_J0_scina %>%
	liana_aggregate() %>%
	# only keep interactions concordant between methods
	filter(aggregate_rank <= 0.05)

chord_freq(J0_liana_object_top)
J0_heatmap_l_scina<-heat_freq(J0_liana_object_top)

J0_s_dotplot<-liana_dotplot(J0_liana_object_top,
							ntop = 10,
							source_groups = c("Inflammation","Proliferation"),
							target_groups = c("DNAdamage","Cell.cycle","Inflammation","Quiescence","Hypoxia"))

png("./01_scripts/output/plots/cell-cell_interaction_plots/by_timepoints/J0_scina_clusters_liana.png",
		width=7500,height=3500, res=300)

J0_s_dotplot

dev.off()

png("./01_scripts/output/plots/cell-cell_interaction_plots/by_timepoints/J0_scina_clusters_heatmap.png",
		width=5000,height=3500, res=300)

J0_heatmap_l_scina

dev.off()

# identify interactions of interest J2 BASED ON TABLE
J2_liana_object_top<- liana_object_J2_scina %>%
	liana_aggregate() %>%
	# only keep interactions concordant between methods
	filter(aggregate_rank <= 0.05)

chord_freq(J2_liana_object_top)
J2_heatmap_l_scina<-heat_freq(J2_liana_object_top)

J2_s_dotplot<-liana_dotplot(J2_liana_object_top,
							ntop = 10,
							source_groups = c("Inflammation"),
							target_groups = c("Inflammation","Cell.cycle","DNAdamage","Quiescence","Cell.cycle","Hypoxia"))

png("./01_scripts/output/plots/cell-cell_interaction_plots/by_timepoints/J2_scina_clusters_liana.png",
		width=5000,height=3500, res=300)

J2_s_dotplot

dev.off()

png("./01_scripts/output/plots/cell-cell_interaction_plots/by_timepoints/J2_scina_clusters_heatmap.png",
		width=5000,height=3500, res=300)

J2_heatmap_l_scina

dev.off()

# identify interactions of interest J7 BASED ON TABLE
J7_liana_object_top<- liana_object_J7_scina %>%
	liana_aggregate() %>%
	# only keep interactions concordant between methods
	filter(aggregate_rank <= 0.05)

chord_freq(J7_liana_object_top)
J7_heatmap_l_scina<-heat_freq(J7_liana_object_top)

J7_s_dotplot<-liana_dotplot(J7_liana_object_top,
							ntop = 10,
							source_groups = c("Hypoxia","Inflammation","Apoptosis"),
							target_groups = c("Inflammation","Cell.cycle","Stemness","Metastasis"))

png("./01_scripts/output/plots/cell-cell_interaction_plots/by_timepoints/J7_scina_clusters_liana.png",
		width=9000,height=3500, res=300)

J7_s_dotplot

dev.off()

png("./01_scripts/output/plots/cell-cell_interaction_plots/by_timepoints/J7_scina_clusters_heatmap.png",
		width=5000,height=3500, res=300)

J7_heatmap_l_scina

dev.off()

# identify interactions of interest J14 BASED ON TABLE
J14_liana_object_top<- liana_object_J14_scina %>%
	liana_aggregate() %>%
	# only keep interactions concordant between methods
	filter(aggregate_rank <= 0.05)

chord_freq(J14_liana_object_top)
J14_heatmap_l_scina<-heat_freq(J14_liana_object_top)

J14_s_dotplot<-liana_dotplot(J14_liana_object_top,
							ntop = 10,
							source_groups = c("Metastasis","Apoptosis"),
							target_groups = c("Metastasis","Cell.cycle","Stemness"))

png("./01_scripts/output/plots/cell-cell_interaction_plots/by_timepoints/J14_scina_clusters_liana.png",
		width=6500,height=3500, res=300)

J14_s_dotplot

dev.off()

png("./01_scripts/output/plots/cell-cell_interaction_plots/by_timepoints/J14_scina_clusters_heatmap.png",
		width=5000,height=3500, res=300)

J14_heatmap_l_scina

dev.off()