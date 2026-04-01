# Spatial transcriptomic landscape of the *Ciona* adult brain: functional zonalisation and cellular composition in a sessile chordate brain and a novel insight into neural gland function

## Overview

This repository contains the analysis code for a 10x Visium spatial transcriptomics study of the adult *Ciona intestinalis* (ascidian) brain. The goal is to map the spatial organization of transcriptionally distinct cell populations in the neural complex—a compact brain region comprising the cerebral ganglion, neural gland, ciliated funnel, and associated structures.

Four Visium capture areas were used, each containing sections from three individual animals placed side by side on the same slide (12 biological replicates total). The pipeline covers spot-level quality control, donor demultiplexing, cross-sample integration, unsupervised clustering, cell type annotation, Gene Ontology enrichment, and preparation of data for XFuse-based super-resolution reconstruction.

---

## Data

The processed datasets used in this study can be downloaded from Figshare: <LINK_HERE>.
  
### Input

Raw input is expected as **Space Ranger output directories**, one per slide:

```
data/
  ciona_brain_nc1/
  ciona_brain_nc2/
  ciona_brain_nc3/
  ciona_brain_nc4/
```

Each directory must contain the standard Space Ranger output structure (i.e., files loadable by `Seurat::Load10X_Spatial()`), including `spatial/tissue_positions_list.csv`, `spatial/scalefactors_json.json`, and the feature-barcode matrix.

### Auxiliary files

The following additional files are required and should be placed under `result/`:

| File | Description |
|------|-------------|
| `result/ciona_human_homo.tsv` | Tab-separated table of Ciona KH2013 gene model to human homolog mappings |
| `result/GO_Cirobu_Aniseedv2019.gaf` | *Ciona robusta* GO annotation file from ANISEED (2019) |

### Notes on slide layout

Each slide captures tissue sections from three donors placed side by side along the x-axis (imagecol). Donor boundaries differ between slides:

| Slide | D1 / D2 boundary | D2 / D3 boundary |
|-------|-------------------|-------------------|
| nc1   | 4000              | 7000              |
| nc2   | 5000              | 8000              |
| nc3   | 4000              | 7000              |
| nc4   | 4000              | 7000              |

---

## Workflow

### Step 1 — Preprocessing (`script/1_preprocessing.R`)

Loads each of the four Visium samples as a Seurat object. Background/empty spots are removed from slides 1, 2, and 4 using imagerow thresholds; slide 3 is unfiltered. Donor identity (D1, D2, D3) is assigned to each spot based on imagecol position. QC metrics (`nCount_Spatial`, `nFeature_Spatial`) are visualized as spatial feature plots and violin plots for each slide. All four objects are merged into a single Seurat object and saved.

### Step 2 — Integration and clustering (`script/2_clustering.r`)

The merged object is normalized using SCTransform. Samples are split by donor and integrated with CCA (`FindIntegrationAnchors` / `IntegrateData`). PCA is run to 50 components; UMAP and graph-based clustering are performed on the top 18 PCs. Clustering is explored at resolutions 0.1, 0.2, 0.4, and 0.8. At resolution 0.2, five clusters are manually annotated:

| Cluster | Cell type |
|---------|-----------|
| 0 | Body wall muscle |
| 1 | Cerebral ganglion |
| 2 | Neural gland duct + Dorsal strand |
| 3 | Neural gland |
| 4 | Ciliated funnel |

Differentially expressed markers (adj. p < 0.05) are identified with `FindAllMarkers` and annotated with human homologs. A heatmap of the top 10 markers per cluster and spatial dimension plots for all slides are generated. An additional pairwise DEG analysis is performed between body wall muscle sub-clusters at resolution 0.4.

### Step 3 — Gene Ontology enrichment (`script/3_Go.r`)

Cluster marker genes (from Step 2) are mapped to base KH2013 gene IDs and tested for GO Biological Process enrichment using `clusterProfiler::compareCluster` against the ANISEED 2019 GAF annotation (p < 0.05, q < 0.2, BH correction). Two figures are produced: a general dotplot showing the top 3 enriched terms per cluster, and a focused bubble plot for selected terms in three clusters (Body wall muscle, Ciliated funnel, Neural gland duct + Dorsal strand).

### Step 4 — Spatial feature visualization (`script/4_spatialFeaturePlot.r`)

Spatial expression maps for two selected genes (`KH2013:KH.C1.215`, `KH2013:KH.C1.14`) are rendered on slide 2 with the inferno colormap. A bar chart summarizes the number of differentially expressed marker genes per annotated cell type.

### Step 5 — XFuse data preparation (`script/Xfuse_split_sildes.r`)

Each slide's `tissue_positions_list.csv` is split into three donor-specific position files by masking out-of-region spots (setting the `in_tissue` flag to 0). The resulting files and corresponding H&E images are converted to XFuse format using `xfuse convert visium` at a scale factor of 0.2, producing 12 individual sample directories (`nc1_1` through `nc4_3`).

---

## Notes

- Each Visium slide captures three individual *Ciona* brains simultaneously. Donor demultiplexing relies on the spatial position of tissue sections along the slide x-axis, not on barcode-based demultiplexing.
- Gene identifiers use the KH2013 nomenclature (e.g., `KH2013:KH.C1.215`). The human homolog mapping table (`ciona_human_homo.tsv`) must be generated separately prior to running Step 2.
- The GO annotation (`GO_Cirobu_Aniseedv2019.gaf`) is specific to *Ciona robusta* and is available from the ANISEED database.
- Script `3_Go.r` contains a minor variable name inconsistency on line 63 (`p5` should be `p1`) that will cause an error if `p5` is not defined in the session. Rename accordingly before running.
- The XFuse commands in `Xfuse_split_sildes.r` reference absolute paths on a specific server (`/home/xzeng/...`). Update these paths to match your environment.

---

## Citation

> [Authors]. *Spatial transcriptomic landscape of the adult brain in Ciona intestinalis*. [Journal], [Year]. [DOI]
