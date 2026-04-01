suppressMessages(library(Seurat))
suppressMessages(library(ggplot2))
suppressMessages(library(patchwork))
suppressMessages(library(cowplot))
library(tidyverse)
set.seed(123)
setwd('~/Desktop/project/Ciona_ST/')

result_dir <- 'result/result1/featurePlot/'

ciona_nc.combined <- readRDS('result/result1/clustering/ciona_nc.combined.rds')


p1 <- SpatialFeaturePlot(ciona_nc.combined, features = 'KH2013:KH.C1.215', images = 'slide2.1', crop = FALSE, pt.size.factor = 1.3,image.alpha = 0.8) +
	  scale_fill_viridis_c(option = "inferno")

ggsave(paste0(file.path(result_dir, "KH.C1.215"), ".pdf"), p1, dpi = 300, width = 6, height = 6)

p1 <- SpatialFeaturePlot(ciona_nc.combined, features = 'KH2013:KH.C1.14', images = 'slide2.1', crop = FALSE, pt.size.factor = 1.3,image.alpha = 0.8) +
    scale_fill_viridis_c(option = "inferno")

ggsave(paste0(file.path(result_dir, "KH.C1.14"), ".pdf"), p1, dpi = 300, width = 6, height = 6)





#CG 108
#BWM 243
#NGD 251
#NG 157
#CF 224

gene_num <- data.frame(cell_type = unique(ciona_nc.combined$cell_type))
gene_num$gene_num <- c(108,157, 251, 243, 224)

p1 <- ggplot(gene_num, aes(x = cell_type, y = gene_num, fill = cell_type)) +
    geom_col() +
    theme_cowplot() +
    ylab('Gene number') +
    labs(fill = NULL) + 
    theme(
        axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size = 18),
        axis.text  = element_text(size = 14),
        strip.text = element_text(size = 14),
        axis.line  = element_line(size = 0.3),
        axis.ticks = element_line(size = 0.3),
        legend.text = element_text(size = 14),
        legend.margin = margin(t = 0, unit = "cm"),
        legend.key.height = unit(0.6, "cm"), 
        legend.spacing.y = unit(1, "cm"),
    )

ggsave(paste0(file.path(result_dir, "gene_num"), ".pdf"), p1, dpi = 300, width = 6, height = 4)