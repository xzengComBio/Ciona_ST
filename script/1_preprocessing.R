suppressMessages(library(Seurat))
suppressMessages(library(ggplot2))
suppressMessages(library(patchwork))

setwd('~/Desktop/project/Ciona_ST/')
set.seed(123)
result_dir <- 'result/result1/preprocessing/'

### function
filter_blank_spots <- function(
  obj,
  slice,
  imagerow_min = NULL,
  imagecol_left_max = NULL,
  imagecol_right_min = NULL
) {
  # Filter out blank/background spots using image coordinates.
  # Thresholds set to NULL will be ignored.

  coords <- obj@images[[slice]]@coordinates

  keep <- rep(TRUE, nrow(coords))

  if (!is.null(imagerow_min)) {
    keep <- keep & (coords$imagerow > imagerow_min)
  }

  if (!is.null(imagecol_left_max)) {
    keep <- keep & (coords$imagecol < imagecol_left_max)
  }

  if (!is.null(imagecol_right_min)) {
    keep <- keep & (coords$imagecol > imagecol_right_min)
  }

  keep_cells <- rownames(coords[keep, ])

  subset(obj, cells = keep_cells)
}

assign_donor_by_imagecol <- function(
  obj,
  slice,
  cut1 = 4000,
  cut2 = 7000,
  donor_colname = "donor"
) {
  # Assign donor IDs based on imagecol thresholds.
  # donor 1: imagecol < cut1
  # donor 2: cut1 <= imagecol <= cut2
  # donor 3: imagecol > cut2

  coords <- obj@images[[slice]]@coordinates

  donor <- ifelse(coords$imagecol < cut1, "D1",
                  ifelse(coords$imagecol > cut2, "D3", "D2"))

  md <- data.frame(donor = donor, row.names = rownames(coords))
  colnames(md) <- donor_colname

  AddMetaData(obj, metadata = md)
}


plot_spatial_qc <- function(
  obj,
  feature,
  prefix,
  result_dir,
  spatial_width = 4,
  spatial_height = 4,
  vln_width = 3,
  vln_height = 4
) {

  # Spatial feature plot
  p_spatial <- SpatialFeaturePlot(
    obj,
    features = feature,
    pt.size.factor = 1,
    image.alpha = 0.5,
    crop = FALSE
  ) &
    theme(
      legend.position = "top",
      legend.direction = "horizontal",
      legend.title.position = "top",
      legend.title.align = 0.5,
      legend.text = element_text(angle = 45, hjust = 1)
    )

  ggsave(
    filename = paste0(result_dir, prefix, "_spatial_", feature, ".pdf"),
    plot = p_spatial,
    width = spatial_width,
    height = spatial_height
  )

  # Violin plot
  p_vln <- VlnPlot(
    obj,
    features = feature,
    pt.size = FALSE
  ) +
    geom_boxplot(width = 0.1, fill = "white", outlier.size = 0) +
    NoLegend() +
    theme(
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank()
    ) +
    scale_x_discrete(labels = "")

  ggsave(
    filename = paste0(result_dir, prefix, "_vln_", feature, ".pdf"),
    plot = p_vln,
    width = vln_width,
    height = vln_height
  )
}
# --------------------------------------------------------------------------------- #





# pre-process of the first sample
nc1 <- Load10X_Spatial('data/ciona_brain_nc1/', slice = 'slide1')
nc1$orig.ident <- 'Slide1'

## remove spots that do not have any cells.
nc1 <- filter_blank_spots(nc1, 
						 slice = "slide1", 
						 imagerow_min = 5000)

nc1 <- assign_donor_by_imagecol(nc1, slice = "slide1", cut1 = 4000, cut2 = 7000)

plot_spatial_qc(
  obj = nc1,
  feature = "nCount_Spatial",
  prefix = "nc1",
  result_dir = result_dir
)

plot_spatial_qc(
  obj = nc1,
  feature = "nFeature_Spatial",
  prefix = "nc1",
  result_dir = result_dir
)
# --------------------- #
# pre-process of the second sample
nc2 <- Load10X_Spatial('data/ciona_brain_nc2/', slice = 'slide2')
nc2$orig.ident <- 'Slide2'
nc2 <- filter_blank_spots(nc2, 
                          slice = "slide2", 
                          imagerow_min = 4600)

nc2 <- assign_donor_by_imagecol(nc2, slice = "slide2", cut1 = 5000, cut2 = 8000)

plot_spatial_qc(
  obj = nc2,
  feature = "nCount_Spatial",
  prefix = "nc2",
  result_dir = result_dir
)

plot_spatial_qc(
  obj = nc2,
  feature = "nFeature_Spatial",
  prefix = "nc2",
  result_dir = result_dir
)

# --------------------- #
# pre-process of the third sample
nc3 <- Load10X_Spatial('data/ciona_brain_nc3/', slice = 'slide3')
nc3$orig.ident <- 'Slide3'

nc3 <- assign_donor_by_imagecol(nc3, slice = "slide3", cut1 = 4000, cut2 = 7000)

plot_spatial_qc(
  obj = nc3,
  feature = "nCount_Spatial",
  prefix = "nc3",
  result_dir = result_dir
)

plot_spatial_qc(
  obj = nc3,
  feature = "nFeature_Spatial",
  prefix = "nc3",
  result_dir = result_dir
)

# --------------------- #
# pre-process of the fourth sample
nc4 <- Load10X_Spatial('data/ciona_brain_nc4/', slice = 'slide4')
nc4$orig.ident <- 'Slide4'
nc4 <- filter_blank_spots(nc4, 
                          slice = "slide4", 
                          imagerow_min = 4400)
nc4 <- assign_donor_by_imagecol(nc4, slice = "slide4", cut1 = 4000 , cut2 = 7000)

plot_spatial_qc(
  obj = nc4,
  feature = "nCount_Spatial",
  prefix = "nc4",
  result_dir = result_dir
)

plot_spatial_qc(
  obj = nc4,
  feature = "nFeature_Spatial",
  prefix = "nc4",
  result_dir = result_dir
)


###
nc_merge <- merge(nc1,y = c(nc2,nc3,nc4), add.cell.ids = c('Slide1','Slide2', 'Slide3', 'Slide4'), project = "CIONA_NC")
saveRDS(nc_merge,paste0(result_dir, "nc_merge.rds"))


## THE END

