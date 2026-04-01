suppressMessages(library(stringr))
suppressMessages(library(clusterProfiler))
suppressMessages(library(ggplot2))
suppressMessages(library(cowplot))

setwd('~/Desktop/project/Ciona_ST/')
result_dir <- 'result/result1/GO/'

ciona_gaf <- read.delim('result/GO_Cirobu_Aniseedv2019.gaf',skip = 1, quote = "")
ciona_gaf_list <- split(ciona_gaf, f= ciona_gaf$Aspect)

go_term2gene_func <- function(x){
     go_term2gene <- data.frame(x$GO.ID,x$DB.Reference...DB.Reference.)
     names(go_term2gene) <- c("go_term", "gene")
     go_term2gene$gene <- str_split(go_term2gene$gene,pattern = ':', simplify = TRUE)[,2]
     
     return(go_term2gene)
}
go_term2gene <- lapply(ciona_gaf_list, FUN = go_term2gene_func)

go_term2gene_all <- data.frame(ciona_gaf$GO.ID, ciona_gaf$DB.Reference...DB.Reference.)
names(go_term2gene_all) <- c("go_term", "gene")
go_term2gene_all$gene <- str_split(go_term2gene_all$gene,pattern = ':', simplify = TRUE)[,2]

go_term2name_func <- function(x){
     go_term2name <- data.frame(x$GO.ID,x$GO.name)
     names(go_term2name) <- c('go_term', 'name')

     return(go_term2name)
}

go_term2name <- lapply(ciona_gaf_list, FUN =go_term2name_func )

go_term2name_all <- data.frame(ciona_gaf$GO.ID,ciona_gaf$GO.name)
names(go_term2name_all) <- c('go_term', 'name')


################################################################################
# Read the differentially expressed genes from Spatial transcriptomic data (resolution 0.4)
DEG_0.2 <- read.table('result/result1/clustering/res2_cluster.markers.txt', sep = ',', header = TRUE)
DEG_0.2$kh_gene <- str_split(DEG_0.2$gene, pattern = ':', simplify = TRUE)[,2]
DEG_0.2_cluster_list <- split(DEG_0.2,f = DEG_0.2$cluster)
DEG_0.2_gene_list <- lapply(DEG_0.2_cluster_list,FUN = function(x){x$kh_gene})


## Go enrichment in Biological process
DEG_0.2_go_bp <- compareCluster(DEG_0.2_gene_list,
                                fun=enricher,
                                pvalueCutoff = 0.05, 
                                qvalueCutoff = 0.2,
                                pAdjustMethod = "BH", 
                                TERM2GENE = go_term2gene$P, 
                                TERM2NAME = go_term2name$P)

p1 <- dotplot(DEG_0.2_go_bp, showCategory=3) + 
    scale_y_discrete(labels=function(x) str_wrap(x, width=30)) + 
    theme(axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)) + 
    theme_cowplot() +
    theme(axis.title.x = element_blank(),
     axis.text.x = element_text(angle = 60, hjust = 1))

ggsave(paste0(file.path(result_dir, "res0.2_GO"), ".pdf"), p5, dpi = 300, width = 6, height = 6)

library(DOSE)

df <- DEG_0.2_go_bp@compareClusterResult
keep_clusters <- c("Body wall muscle", "Ciliated funnel", "Neural gland duct +   Dorsal strand")
keep_terms <- c(
    "muscle contraction",
    "skeletal muscle tissue development",
    "cardiac muscle contraction",
    "cilium movement",
    "cilium assembly",
    "outer dynein arm assembly",
    "cellular defense response",
    "macrophage activation involved in immune response")


term_labels <- str_wrap(keep_terms, width = 24)

df2 <- df %>%
    filter(Cluster %in% keep_clusters, Description %in% keep_terms) %>%
    mutate(
        Cluster = factor(Cluster, levels = keep_clusters),
        Description = str_wrap(Description, width = 24),
        Description = factor(Description, levels = rev(term_labels)),
        GeneRatio = parse_ratio(GeneRatio)
    )


p1 <- ggplot(df2, aes(x = Cluster, y = Description)) +
    geom_point(aes(size = GeneRatio, fill = -log10(p.adjust)),
               shape = 21, color = "black", stroke = 0.4) +
    scale_size(range = c(2.5, 8)) +
    scale_fill_gradient(low = "#F6B4AE", high = "#D85C57") +
    theme_cowplot(font_size = 10) +
    theme(
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 15, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 8)
    )
ggsave(paste0(file.path(result_dir, "res0.2_selected_GO"), ".pdf"), p1, dpi = 300, width = 5, height = 5)

