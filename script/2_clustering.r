suppressMessages(library(Seurat))
suppressMessages(library(ggplot2))
suppressMessages(library(patchwork))
suppressMessages(library(cowplot))
library(tidyverse)
library(clustree)
setwd('~/Desktop/project/Ciona_ST/')
set.seed(123)
result_dir <- 'result/result1/clustering/'

nc_merge <- readRDS('result/result1/preprocessing/nc_merge.rds')
nc_merge <- SCTransform(nc_merge, assay = "Spatial", verbose = FALSE)
nc_merge.list <- SplitObject(nc_merge, split.by = 'donor')
features <- SelectIntegrationFeatures(object.list = nc_merge.list)
ciona_nr.anchors <- FindIntegrationAnchors(object.list = nc_merge.list, anchor.features = features)
ciona_nc.combined <- IntegrateData(anchorset = ciona_nr.anchors)
ciona_nc.combined <- ScaleData(ciona_nc.combined, verbose = FALSE)
ciona_nc.combined <- RunPCA(ciona_nc.combined, npcs = 50, verbose = FALSE)


###Plot PCA distribution across individuals
p1 <- DimPlot(ciona_nc.combined, reduction = 'pca', group.by = 'donor') +
    theme_cowplot() +
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          axis.title = element_text(face = "bold",size = rel(1)),
          plot.title = element_blank(),
          legend.position = c(0.8,0.2),
          legend.margin=margin(t = 0, unit='cm'))

ggsave(filename = paste0(result_dir, 'pca_after_cca_individual.pdf'), p1, width=3.2, height=3)

p2 <- DimPlot(ciona_nc.combined, reduction = 'pca', group.by = 'orig.ident') +
    theme_cowplot() +
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          axis.title = element_text(face = "bold",size = rel(1)),
          plot.title = element_blank(),
          legend.position = c(0.7,0.2),
          legend.margin=margin(t = 0, unit='cm'))

ggsave(filename = paste0(result_dir, 'pca_after_cca_slide.pdf'), p2, width=3.2, height=3)




##
ciona_nc.combined <- RunUMAP(ciona_nc.combined, dims = 1:18)

ciona_nc.combined <- FindNeighbors(ciona_nc.combined, dims = 1:18)

resolutions <- c(0.1,0.2,0.4,0.8)
for (i in resolutions){
     ciona_nc.combined <- FindClusters(ciona_nc.combined,resolution = i)
}



saveRDS(ciona_nc.combined, paste0(result_dir, 'ciona_nc.combined.rds'))

### Vis


p3 <- DimPlot(
    ciona_nc.combined,
    reduction = "umap",
    group.by = "integrated_snn_res.0.2"
) +
  theme_cowplot() +
  theme(
    axis.title = element_text(face = "bold", size = rel(1)),
    axis.text  = element_blank(),
    axis.ticks = element_blank(),
    axis.line  = element_blank(),
    plot.title = element_blank(),
    legend.margin = margin(t = 0, unit = "cm")
  )

ggsave(filename = paste0(result_dir, 'umap_clustering_res2.pdf'), p3, width=3.3, height=3)

anno_ids <- c('Body wall muscle', 'Cerebral ganglion', 'Neural gland duct + \n Dorsal strand', 'Neural gland', 'Ciliated funnel')
names(anno_ids) <- levels(ciona_nc.combined)
ciona_nc.combined <- RenameIdents(ciona_nc.combined, anno_ids)
ciona_nc.combined$cell_type <- Idents(ciona_nc.combined)

p3 <- DimPlot(
    ciona_nc.combined,
    reduction = "umap",
    group.by = 'cell_type'
) +
  theme_cowplot() +
  theme(
    axis.title = element_text(face = "bold", size = rel(1)),
    axis.text  = element_blank(),
    axis.ticks = element_blank(),
    axis.line  = element_blank(),
    plot.title = element_blank(),
    legend.margin = margin(t = 0, unit = "cm"),
    legend.key.height = unit(0.8, "cm"), 
    legend.spacing.y = unit(1, "cm"),

  )
ggsave(filename = paste0(result_dir, 'umap_anno.pdf'), p3, width=4.8, height=3)


nc_merge <- AddMetaData(nc_merge,ciona_nc.combined@meta.data)
Idents(nc_merge) <- nc_merge$cell_type
p1 <- SpatialDimPlot(nc_merge, images = c('slide1'), crop = FALSE, pt.size.factor = 1.3,image.alpha = 0.4) 
ggsave(filename = paste0(result_dir, 'spatial_dim_anno_slide1.pdf'), p1, width=4.8, height=3)
p2 <- SpatialDimPlot(nc_merge, images = c('slide2'), crop = FALSE, pt.size.factor = 1.3,image.alpha = 0.4) 
ggsave(filename = paste0(result_dir, 'spatial_dim_anno_slide2.pdf'), p2, width=4.8, height=3)
p3 <- SpatialDimPlot(nc_merge, images = c('slide3'), crop = FALSE, pt.size.factor = 1.3,image.alpha = 0.4) 
ggsave(filename = paste0(result_dir, 'spatial_dim_anno_slide3.pdf'), p3, width=4.8, height=3)
p4 <- SpatialDimPlot(nc_merge, images = c('slide2'), crop = FALSE, pt.size.factor = 1.3,image.alpha = 0.4) 
ggsave(filename = paste0(result_dir, 'spatial_dim_anno_slide4.pdf'), p4, width=4.8, height=3)


ciona_nc.0.2.markers <- FindAllMarkers(ciona_nc.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
ciona_nc.0.2.markers <- ciona_nc.0.2.markers[ciona_nc.0.2.markers$p_val_adj < 0.05,]


ciona_human_genes <- read.table('result/ciona_human_homo.tsv', sep = '\t', header = TRUE)
ciona_human_genes$gene_model <- paste0('KH2013:',ciona_human_genes$gene_model)
colnames(ciona_human_genes)[2] <- 'human_homolog'
ciona_nc.0.2.markers <- left_join(ciona_nc.0.2.markers, ciona_human_genes[,c('gene_model', 'human_homolog')], by = c('gene' = 'gene_model'))
ciona_nc.0.2.markers.top10 <- ciona_nc.0.2.markers %>% group_by(cluster) %>% slice_max(n=10, order_by = avg_log2FC)

p1 <- DoHeatmap(
  ciona_nc.combined,
  features = ciona_nc.0.2.markers.top10$gene,
  group.by = "cell_type",
  label = FALSE
) +
  theme(
    axis.ticks.x = element_blank()
  )

ggsave(
  filename = paste0(result_dir, "heatmap_top10.pdf"),
  plot = p1,
  width = 5, height =6
)

ciona_nc.0.2.markers$cluster <- gsub("[\r\n]+", " ", ciona_nc.0.2.markers$cluster)

write.table(ciona_nc.0.2.markers, paste0(result_dir, 'res2_cluster.markers.txt'),quote = FALSE,sep = ',', col.names = TRUE, row.names = FALSE)

saveRDS(ciona_nc.combined, paste0(result_dir, 'ciona_nc.combined.rds'))
saveRDS(nc_merge, paste0(result_dir, 'nc_merge_bf_cca.rds'))


p1 <- SpatialDimPlot(nc_merge, group.by = 'integrated_snn_res.0.4', images = 'slide2', crop = FALSE, pt.size.factor = 1.3,image.alpha = 0.8) + 
      guides(fill=guide_legend(title="Clusters",override.aes = list(size = 5)))

ggsave(
  filename = paste0(result_dir, "slide2_spatial_dim_res04.pdf"),
  plot = p1,
  width = 5, height =3
)

Idents(ciona_nc.combined) <- ciona_nc.combined$integrated_snn_res.0.4

bodywallmuscle.markers <- FindMarkers(ciona_nc.combined, ident.1 = 1, ident.2 = 3)

bodywallmuscle.markers <- bodywallmuscle.markers[bodywallmuscle.markers$p_val_adj < 0.05, ]

write.table(bodywallmuscle.markers, paste0(result_dir, 'bodywallmuscle.markers.txt'),quote = FALSE,sep = ',', col.names = TRUE, row.names = FALSE)

p1 <- SpatialFeaturePlot(nc_merge, features = 'KH2013:KH.C12.521', images = 'slide2', crop = FALSE, pt.size.factor = 1.3,image.alpha = 0.8)
p2 <- SpatialFeaturePlot(nc_merge, features = 'KH2013:KH.C7.598', images = 'slide2', crop = FALSE, pt.size.factor = 1.3,image.alpha = 0.8)
ggsave(
  filename = paste0(result_dir, "bodywallMuscle_markers.pdf"),
  plot = p1 | p2,
  width = 10, height =5
)


# ciona_nc.merge <- AddMetaData(ciona_nc.merge,ciona_nc.combined@meta.data)
# Idents(ciona_nc.merge) <- ciona_nc.merge$integrated_snn_res.0.2
# SpatialDimPlot(ciona_nc.merge, images = 'sample1', crop = FALSE, pt.size.factor = 1.3,image.alpha = 0.8) +  guides(fill=guide_legend(title="Clusters",override.aes = list(size = 5)))
# SpatialDimPlot(ciona_nc.merge, images = 'sample2', crop = FALSE, pt.size.factor = 1.3,image.alpha = 0.8) +  guides(fill=guide_legend(title="Clusters",override.aes = list(size = 5)))
# SpatialDimPlot(ciona_nc.merge, images = 'sample3', crop = FALSE, pt.size.factor = 1.3,image.alpha = 0.8) +  guides(fill=guide_legend(title="Clusters",override.aes = list(size = 5)))
# SpatialDimPlot(ciona_nc.merge, images = 'sample4', crop = FALSE, pt.size.factor = 1.3,image.alpha = 0.8) +  guides(fill=guide_legend(title="Clusters",override.aes = list(size = 5)))

# Idents(ciona_nc.merge) <- ciona_nc.merge$integrated_snn_res.0.4
# SpatialDimPlot(ciona_nc.merge, images = 'sample1', crop = FALSE, pt.size.factor = 1.3,image.alpha = 0.8) +  guides(fill=guide_legend(title="Clusters",override.aes = list(size = 5)))
# SpatialDimPlot(ciona_nc.merge, images = 'sample2', crop = FALSE, pt.size.factor = 1.3,image.alpha = 0.8) +  guides(fill=guide_legend(title="Clusters",override.aes = list(size = 5)))
# SpatialDimPlot(ciona_nc.merge, images = 'sample3', crop = FALSE, pt.size.factor = 1.3,image.alpha = 0.8) +  guides(fill=guide_legend(title="Clusters",override.aes = list(size = 5)))
# SpatialDimPlot(ciona_nc.merge, images = 'sample4', crop = FALSE, pt.size.factor = 1.3,image.alpha = 0.8) +  guides(fill=guide_legend(title="Clusters",override.aes = list(size = 5)))

# Idents(ciona_nc.combined) <- ciona_nc.combined$integrated_snn_res.0.2
# ciona_nc.0.2.markers <- FindAllMarkers(ciona_nc.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# ciona_nc.0.2.markers <- filter(ciona_nc.0.2.markers,p_val_adj < 0.05)
# write.table(ciona_nc.0.2.markers,'result/clustering/ciona_nc.0.2.markers.csv',quote = FALSE,sep = ',', col.names = TRUE, row.names = FALSE)
# ciona_nc.0.2.markers.top20 <- ciona_nc.0.2.markers %>% group_by(cluster) %>% slice_max(n=20, order_by = avg_log2FC)
# DoHeatmap(ciona_nc.combined, features = ciona_nc.0.2.markers.top20$gene) + NoLegend()

# Idents(ciona_nc.combined) <- ciona_nc.combined$integrated_snn_res.0.4
# ciona_nc.0.4.markers <- FindAllMarkers(ciona_nc.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# ciona_nc.0.4.markers <- filter(ciona_nc.0.4.markers,p_val_adj < 0.05)
# write.table(ciona_nc.0.4.markers,'result/clustering/ciona_nc.0.4.markers.csv',quote = FALSE,sep = ',', col.names = TRUE, row.names = FALSE)

# ciona_nc.0.4.markers.top20 <- ciona_nc.0.4.markers %>% group_by(cluster) %>% slice_max(n=20, order_by = avg_log2FC)
# DoHeatmap(ciona_nc.combined, features = ciona_nc.0.4.markers.top50$gene) + NoLegend()

# ciona_nc.0.8.markers <- FindAllMarkers(ciona_nc.combined, logfc.threshold = 0.25, min.pct = 0.25, only.pos = TRUE)



# ciona_human_genes <- read.table('../TF/pythonCodes/codes/results/ciona_human_homo.tsv', sep = '\t', header = TRUE)
# ciona_human_genes$gene_model <- paste0('KH2013:',ciona_human_genes$gene_model)
# rownames(ciona_human_genes) <- ciona_human_genes$gene_model

# ciona_nc.0.2.markers$human_homologs <- ciona_human_genes[ciona_nc.0.2.markers$gene,]$gene_name
# write.table(ciona_nc.0.2.markers,'result/clustering/ciona_nc.0.2.markers.csv',quote = FALSE,sep = ',', col.names = TRUE, row.names = FALSE)

# ciona_nc.0.4.markers <- read.table('result/clustering/ciona_nc.0.4.markers.csv', sep = ',', header = TRUE)
# ciona_nc.0.4.markers$human_homologs <- ciona_human_genes[ciona_nc.0.4.markers$gene,]$gene_name
# write.table(ciona_nc.0.4.markers,'result/clustering/ciona_nc.0.4.markers.csv',quote = FALSE,sep = ',', col.names = TRUE, row.names = FALSE)

# ciona_nc.0.8.markers <- read.table('result/clustering/ciona_nc.0.8.markers.csv', sep = ',', header = TRUE)
# ciona_nc.0.8.markers$human_homologs <- ciona_human_genes[ciona_nc.0.8.markers$gene,]$gene_name
# write.table(ciona_nc.0.8.markers,'result/clustering/ciona_nc.0.8.markers.csv',quote = FALSE,sep = ',', col.names = TRUE, row.names = FALSE)



# mySpatialFeatureaPlot <- function(sp,feature, save_image=TRUE){

#     p1 <- SpatialFeaturePlot(sp,features = feature,images = 'sample1', crop = FALSE, pt.size.factor = 1.3,image.alpha = 0.8)
#     p2 <- SpatialFeaturePlot(sp,features = feature,images = 'sample2', crop = FALSE, pt.size.factor = 1.3,image.alpha = 0.8)
#     p3 <- SpatialFeaturePlot(sp,features = feature,images = 'sample3', crop = FALSE, pt.size.factor = 1.3,image.alpha = 0.8) + theme(legend.position = 'bottom')
#     p4 <- SpatialFeaturePlot(sp,features = feature,images = 'sample4', crop = FALSE, pt.size.factor = 1.3,image.alpha = 0.8) + theme(legend.position = 'bottom')

#     p12 <- p1 | p2
#     p34 <- p3 | p4

#     if (save_image == TRUE){
#         ggsave(paste0('result/spatialPlot/',substring(feature, first = 8),'_sample12','.pdf'),
#                 p12,
#                 dpi = 200,
#                 units = 'cm',
#                 width = 30,
#                 height = 16)

#         ggsave(paste0('result/spatialPlot/',substring(feature, first = 8),'_sample34','.pdf'),
#                 p34,
#                 dpi = 200,
#                 units = 'cm',
#                 width = 30,
#                 height = 16)
#     }

# }

# nc_merge_sct <- AddMetaData(nc_merge_sct, nc_combined@meta.data)
# saveRDS(nc_merge_sct, file = 'result/rds/nc_merge_sct.rds')

# ##################################################

# opsin_genes <- data.frame(gene_name= c('Ci-opsin1', 'Ci-opsin2', 'Ci-opsin3', 'Ci-opsin5', 'Ci-opsin6', 'Ci-Nut2', 'Ci-Nut1'),
#                              kh_model = c('KH.L171.13', 'KH.L38.6','KH.L57.28', 'KH.C9.770', 'KH.C7.385', 'KH.C14.516', 'KH.C14.4'))
# opsin_genes$kh2013 <- paste0('KH2013:',opsin_genes$kh_model)
# opsin_genes_counts <- nc_merge_sct@assays$Spatial@data[opsin_genes$kh2013,]

# opsin_genes$num_spot_expressed <- rowSums(opsin_genes_counts > 0)
# opsin_genes$mean_in_expressed_spot <- rowSums(opsin_genes_counts)/rowSums(opsin_genes_counts > 0)
# opsin_genes$mean_in_expressed_spot <- round(opsin_genes$mean_in_expressed_spot,digits = 2)
# opsin_genes$max_expressed <- apply(opsin_genes_counts,MARGIN = 1, max)

# # check the expression level of Ci-opsin2 in cerebral ganglion.
# cluster0 <- subset(larva_ner, subset = integrated_snn_res.0.2 == '0')
# Idents(cluster0, WhichCells(object = cluster0, expression = `KH2013:KH.L38.6` > 0, slot = 'counts')) <- 'Ci-opsin2_positive'
# Idents(cluster0, WhichCells(object = cluster0, expression = `KH2013:KH.L38.6` <= 0, slot = 'counts')) <- 'Ci-opsin2_negative'

# cluster0 <- PrepSCTFindMarkers(cluster0)
# genes <- FindMarkers(cluster0, ident.1 = 'Ci-opsin2_positive', ident.2 = 'Ci-opsin2_negative', slot = 'data')


# gene_for_vln <- subset(opsin_genes,  num_spot_expressed > 4)

# VlnPlot(nc_merge_sct,gene_for_vln$kh2013,pt.size = FALSE, stack = TRUE, flip = TRUE, group.by = 'integrated_snn_res.0.2') +
#     geom_boxplot(width=0.1,fill="white",outlier.size = 0) +
#     NoLegend() +
#     theme(axis.title.x=element_blank())


# myFeaturePlot <- function(sce,feature, save_image =TRUE){

#     p1 <- FeaturePlot(sce, features = feature) + 
#             theme_bw() +
#             theme(panel.grid.major=element_blank(),
#                   panel.grid.minor=element_blank(),
#                   axis.title = element_text(face = "bold",size = rel(1)),
#                   plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
#                   legend.margin=margin(t = 0, unit='cm')) +
#             ggtitle(feature)

#     if (save_image == TRUE){
#         ggsave(paste0('result/opsin_genes/',substring(feature, first = 8),'_featurePlot','.pdf'),
#                 p1,
#                 dpi = 200,
#                 units = 'cm',
#                 width = 17,
#                 height = 15)
#     }
# }

# for (x in opsin_genes$kh2013){

#     myFeaturePlot(nc_combined, feature = x)

# }


# nc_combined <- readRDS('result/rds/nc_combined.rds')
# cluster0 <- subset(nc_combined, integrated_snn_res.0.8 == c(0,5,9,10))
# markers_sep <- FindAllMarkers(cluster0,only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# markers_sep <- subset(markers_sep, p_val_adj < 0.05)

# markers_sep$human_homologs <- ciona_human_genes[markers_sep$gene,]$gene_name


# ##get highly variable genes from samples
# ciona_merge <- readRDS('../result/rds/nc_merge.rds')
# ciona_merge <- SCTransform(ciona_merge, assay = 'Spatial', verbose = FALSE, variable.features.n = 1000)
# hvgs <- VariableFeatures(ciona_merge)

# opsin_genes <- c('KH.L171.13', 'KH.L38.6','KH.L57.28', 'KH.C9.770', 'KH.C7.385', 'KH.C14.516', 'KH.C14.4')
# hvgs_opsin <- c(hvgs, paste0("KH2013:",opsin_genes))

# nc_regulon_info <- readRDS('../../SCENIC/results/R_downstream/ciona/nc_regulon_info.rds')
# regulon_tf <- nc_regulon_info$tf_info$tf_name

# hvgs_opsin_tf <- c(hvgs_opsin, regulon_tf)
# hvgs_opsin_tf <- unique(hvgs_opsin_tf)
# hvgs_opsin_tf_merge <- paste(hvgs_opsin_tf, collapse = "|")