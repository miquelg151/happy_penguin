## INTEGRATION ANALYSIS


#pkgs
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


J0 <- readRDS("./01_scripts/saved_objects/J0_rescaled.rds")
J2 <- readRDS("./01_scripts/saved_objects/J2_rescaled.rds")
J7 <- readRDS("./01_scripts/saved_objects/J7_rescaled.rds")
J14 <- readRDS("./01_scripts/saved_objects/J14_rescaled.rds")

# add cc.diff in each object to be able to regress it after merging ----------
J0$CC.Dif <- J0$S.Score - J0$G2M.Score
J2$CC.Dif <- J2$S.Score - J2$G2M.Score
J7$CC.Dif <- J7$S.Score - J7$G2M.Score
J14$CC.Dif <- J14$S.Score - J14$G2M.Score


obj <- merge (x=J0, y =list(J2, J7, J14))
obj[["RNA"]] <- as(obj[["RNA"]], Class = "Assay5")

dir.create("./01_scripts/output/plots/integration/")

obj_unintegrated <- NormalizeData(obj)
obj_unintegrated <- FindVariableFeatures(obj_unintegrated, nfeatures = 5000)
obj_unintegrated <- ScaleData(obj_unintegrated, features = rownames(obj_unintegrated), vars.to.regress = c("percent.mt", "CC.Dif"))
# save scaled object (very long to be generated)
saveRDS(obj_unintegrated, "./01_scripts/saved_objects/object.scaled.rds")
# --------------------------------------------------------------------------

obj_unintegrated <- readRDS("./01_scripts/saved_objects/object.scaled.rds")
obj.pca <- RunPCA(obj_unintegrated)


obj1 <- FindNeighbors(obj.pca)
obj1 <- FindClusters(obj1, resolution = 0.5, cluster.name = "unintegrated_clusters")
obj1 <- RunUMAP(obj1, dims = 1:20)
obj1$orig.ident <- factor(x  = obj1$orig.ident, levels = c("J0", "J2","J7", "J14"))
p1 <- DimPlot(obj1, split.by = "orig.ident")+
	ggtitle("unintegrated split by timepoint")
p1.1 <- DimPlot(obj1, group.by = "orig.ident")+
	ggtitle("unintegrated group by timepoint")



# Harmony integration

obj2 <- IntegrateLayers(object = obj.pca, method = HarmonyIntegration,
											 orig.reduction = "pca", new.reduction = "integrated.harmony",
											 verbose = F)
obj2 <- FindNeighbors(obj2, reduction = "integrated.harmony", dims = 1:20)
obj2 <- FindClusters(obj2, resolution = 0.5, cluster.name = "harmony_clusters")
obj2 <- RunUMAP(obj2, reduction = "integrated.harmony", dims = 1:20, reduction.name = "umap.harmony")
obj2$orig.ident <- factor(x = obj2$orig.ident, levels = c("J0", "J2", "J7", "J14"))
p2 <- DimPlot(obj2, reduction = "umap.harmony", group.by = "orig.ident") +
	ggtitle("harmony integrated group by timepoint")
p2.1 <- DimPlot(obj2, reduction = "umap.harmony", split.by = "orig.ident") +
	ggtitle("harmony integrated split by timepoint")


saveRDS(obj2, "./01_scripts/saved_objects/object.int.umap.rds")
obj2 <- readRDS("./01_scripts/saved_objects/object.int.umap.rds")


pdf("./01_scripts/output/plots/integration/integration_UMAPs.pdf",
		paper = "A4")

p1.1
p2

dev.off()

png("./01_scripts/output/plots/integration/unintegrated_integration_UMAPs_split_by.png",
		width=2500,height=1500, res=300)

p1

dev.off()

png("./01_scripts/output/plots/integration/integrated_integration_UMAPs_split_by.png",
		width=2500,height=1500, res=300)
p2.1

dev.off()

# rejoin layers before annotation and differential gene expression analysis -
obj2<-readRDS("./01_scripts/saved_objects/object.int.umap.rds")
obj.int.join <- JoinLayers(obj2)

saveRDS(obj.int.join, "./01_scripts/saved_objects/object.int.joinedLayers.rds")

obj.markers <- FindAllMarkers(obj.int.join,
														 logfc.threshold = 0.15,
														 min.pct = 0.20,
														 only.pos = TRUE)
sig.obj.markers <- obj.markers %>%
	group_by(cluster) %>%
	filter(avg_log2FC > 0.99 & p_val_adj < 0.05) %>%
	arrange(-avg_log2FC, .by_group = TRUE)
write.xlsx(as.data.frame(sig.obj.markers), "./01_scripts/output/tables/sig.markers.integration.xlsx")


#### Automated Cell Type Annotation with SCINA -------------------------------##

## signatures from CancerSEA (14 functional states) ----------------------------

# import gene signatures -----------------------------------------------------      
obj.int.join <- readRDS("./01_scripts/saved_objects/object.int.joinedLayers.rds")

file <- "03_other_files/14signatures.xlsx"
signatures <- import_list(file)
gene.signatures <- lapply(signatures, function(df) df$GeneName)

obj.int.data <- obj.int.join@assays$RNA$scale.data 

obj.int.scina <- SCINA(obj.int.data,
											 gene.signatures, 
											 max_iter = 100, 
											 convergence_n = 10, 
											 convergence_rate = 0.999, 
											 sensitivity_cutoff = 0.9, 
											 rm_overlap=FALSE, 
											 allow_unknown=TRUE) 

obj.int.join$scina_labels <- obj.int.scina$cell_labels

metadata <- obj.int.join@meta.data[, c("seurat_clusters", "scina_labels")]
label_cluster <- as.data.frame(table(obj.int.join$seurat_clusters, obj.int.join$scina_labels)) %>%
	arrange(Var1)
colnames(label_cluster) <- c("Seurat_cluster", "Predicted_signature", "Frequency")

# add dge markers identified with seurat clusters ----------------------------
# open markers file
library(readxl)
sig.markers <- read_excel("./01_scripts/output/tables/sig.markers.integration.xlsx")


genes.clusters <-  sig.markers%>%
	group_by(cluster) %>%
	summarise(sig.seurat.markers = paste(gene, collapse = ","))

label_cluster <- merge(label_cluster, genes.clusters, by.x="Seurat_cluster", by.y="cluster", all.x = TRUE)

write.xlsx(label_cluster, "./01_scripts/output/tables/seurat_cluster_scina_labels.integrated.xlsx", row.names = F)

# Define the colors of each cluster
cols <- c("Stemness" = 'darkorange3', "Quiescence" = "blue4", "Proliferation" = "orange1",
					"Metastasis" = "turquoise", "Invasion" = "plum", "Inflammation" = "khaki2", "Hypoxia" = "green", 
					"EMT" = "darkgreen", "DNArepair" = "mediumvioletred", "DNAdamage" = "palevioletred4", 
					"Differentiation" = "antiquewhite3", "Cell.cycle" = "slateblue1", "Apoptosis" = "gray41", "Angiogenesis" = "red") 


Idents(obj.int.join) <- "scina_labels"
obj.annot <- DimPlot(obj.int.join, reduction = "umap.harmony", cols = cols, split.by = "orig.ident") + 
	labs(title = "J0 - 14 cancer signatures")
png("./01_scripts/output/plots/integration/integrative_UMAP_scina-14sig.png", width=2500,height=1500, res=300)
obj.annot
dev.off()

## LIANA ANALYSIS
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

liana_object <- liana_wrap(obj.int.join,
													 resource = 'custom',
													 external_resource = op_resource,
													 idents_col = "seurat_clusters",
													 min_cells = 45,
													 supp_columns = c("ligand.expr", "receptor.expr", "ligand.stat", "receptor.stat", "ligand.pval", "receptor.pval", "ligand.FDR", "receptor.FDR"),
													 assay = "RNA",
													 expr_prop = 0.2)

liana_object_scina <- liana_wrap(obj.int.join,
													 resource = 'custom',
													 external_resource = op_resource,
													 idents_col = "scina_labels",
													 min_cells = 45,
													 supp_columns = c("ligand.expr", "receptor.expr", "ligand.stat", "receptor.stat", "ligand.pval", "receptor.pval", "ligand.FDR", "receptor.FDR"),
													 assay = "RNA",
													 expr_prop = 0.2)

saveRDS(liana_object, "./01_scripts/saved_objects/liana_object_integrated.RDS")
saveRDS(liana_object_scina, "./01_scripts/saved_objects/liana_object_integrated_scina.RDS")


liana_object<-readRDS("./01_scripts/saved_objects/liana_object_integrated.RDS")
liana_object_scina <- readRDS("./01_scripts/saved_objects/liana_object_integrated_scina.RDS")


# identify interactions of interest
liana_object_top <- liana_object %>%
	liana_aggregate() %>%
	# only keep interactions concordant between methods
	filter(aggregate_rank <= 0.05) %>% 
	# keep top 20 interactions (regardless of cell type)
	distinct_at(c("source","target","ligand.complex","receptor.complex")) %>% 
	head(50)

# identify interactions of interest
liana_object_top_scina <- liana_object_scina %>%
	liana_aggregate() %>%
	# only keep interactions concordant between methods
	filter(aggregate_rank <= 0.05) %>% 
	# keep top 20 interactions (regardless of cell type)
	distinct_at(c("source","target","ligand.complex","receptor.complex","aggregate_rank")) %>% 
	head(50)

write.xlsx(liana_object_top, "./01_scripts/output/tables/cell-cell_interaction/integrative_analysis_interactions.xlsx")
write.xlsx(liana_object_top_scina, "./01_scripts/output/tables/cell-cell_interaction/integrative_analysis_interactions_scina.xlsx")


# identify interactions of interest
liana_object<- liana_object %>%
	liana_aggregate() %>%
	# only keep interactions concordant between methods
	filter(aggregate_rank <= 0.05)

##liana_dotplot(liana_object,
##							ntop = 10,
##							source_groups = c("0"))
##liana_dotplot(liana_object,
##							ntop = 10,
##							source_groups = c("1"))
##liana_dotplot(liana_object,
##							ntop = 10,
##							source_groups = c("2"))
##liana_dotplot(liana_object,
##							ntop = 10,
##							source_groups = c("3"))
##liana_dotplot(liana_object,
##							ntop = 10,
##							source_groups = c("4"))
##liana_dotplot(liana_object,
##							ntop = 10,
##							source_groups = c("5"))
##liana_dotplot(liana_object,
##							ntop = 10,
##							source_groups = c("6"))
##liana_dotplot(liana_object,
##							ntop = 10,
##							source_groups = c("7"))
##liana_dotplot(liana_object,
##							ntop = 10,
##							source_groups = c("8"))


heatmap_l<- heat_freq(liana_object)

p_l_dotplot<-liana_dotplot(liana_object,
							ntop = 10,
							source_groups = c("6","8","4","3"),
							target_groups = c("1","8","3","7"))

dir.create("./01_scripts/output/plots/cell-cell_interaction_plots")

png("./01_scripts/output/plots/cell-cell_interaction_plots/seurat_clusters_liana.png",
		width=4500,height=3500, res=300)

p_l_dotplot

dev.off()

png("./01_scripts/output/plots/cell-cell_interaction_plots/seurat_clusters_heatmap.png",
		width=4500,height=3500, res=300)

heatmap_l

dev.off()

# identify interactions of interest
liana_object_scina<- liana_object_scina %>%
	liana_aggregate() %>%
	# only keep interactions concordant between methods
	filter(aggregate_rank <= 0.05)

heatmap_scina<- heat_freq(liana_object_scina)
				 
p_l_dotplot<-liana_dotplot(liana_object_scina,
							ntop = 10,
							source_groups = c("Inflammation"),
							target_groups = c("Differentiation","DNAdamage","Cell.cycle","Hypoxia","Stemness"))

png("./01_scripts/output/plots/cell-cell_interaction_plots/scina_clusters_liana.png",
		width=5000,height=3500, res=300)

p_l_dotplot

dev.off()

png("./01_scripts/output/plots/cell-cell_interaction_plots/scina_clusters_liana_heatmap.png",
		width=5000,height=3500, res=300)

heatmap_scina

dev.off()


#------------------------------- DISCUSSION PLOT

# revCSC --> ANXA1, CLU, LY6A (NO APAREIX), KLF4
# proCSC --> LGR5, OLFM4, EPHB2, LRIG1

library(SCpubr)
my.col =c("lightgrey","#6D53DC","#DC5356")


p<-FeaturePlot(obj2,
							 features = c("ANXA1",
							 						 "CLU",
							 						 "KLF4",
							 						 "YAP1",
							 						 "NFKB1",
							 						 "TGFB1"),
						cols = my.col,
						pt.size = 1,
							 split.by = "orig.ident",
						order = T) &
	theme(legend.position = "right")

dir.create("./01_scripts/output/plots/integration/genes_discussion")
png("./01_scripts/output/plots/integration/genes_discussion/revCSC.png",
		width=5000,height=3500, res=300)

p

dev.off()

p2<-FeaturePlot(obj2,
							 features = c("EPHB2",
							 						 "LRIG1",
							 						 "CTNNB1",
							 						 "MAPK1",
							 						 "MYC"),
						cols= my.col,	
							 split.by = "orig.ident",
						order = T)&
	theme(legend.position = "right")

FeaturePlot(obj2,
						features = c("EPHB2",
												 "LRIG1",
												 "CTNNB1",
												 "MAPK1",
												 "MYC"),
						pt.size = 1,
						split.by = "orig.ident",
						order = T)&
	theme(legend.position = "right")

p2

png("./01_scripts/output/plots/integration/genes_discussion/proCSC.png",
		width=5000,height=3500, res=300)

p2

dev.off()

p3<-FeaturePlot(obj2,
						features = c("ANXA1",
												 "TGFB1"),
						pt.size = 1,
						cols = my.col,
						split.by = "orig.ident",
						order = T)&
	theme(legend.position = "right")
p3
