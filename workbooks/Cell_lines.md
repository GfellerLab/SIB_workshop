-   [Suggested course structute:](#suggested-course-structute)
-   [TO DO:](#to-do)
-   [Single-cell level](#single-cell-level)
    -   [Standard downstream analysis](#standard-downstream-analysis)
        -   [Pre-processing](#pre-processing)
        -   [UMAP (non-linear dimensionality
            reduction)](#umap-non-linear-dimensionality-reduction)
        -   [Clustering](#clustering)
        -   [Differential expression
            analysis](#differential-expression-analysis)
-   [Data simplification (coarse-graining) – Construction of
    *metacells*](#data-simplification-coarse-graining-construction-of-metacells)
    -   [Downstream analysis of
        metacells](#downstream-analysis-of-metacells)
    -   [Sample-weighted downstream analysis of
        *metacells*](#sample-weighted-downstream-analysis-of-metacells)
        -   [Pre-processing](#pre-processing-1)
        -   [Differential expression analysis in
            metacells](#differential-expression-analysis-in-metacells)
    -   [Standard downstream analysis](#standard-downstream-analysis-1)

Suggested course structute:
===========================

We first run a **standard scRNA-seq data analysis pipeline** (i.e., data
normalization, feature selection, dimensionality reduction,
visualization, clustering and differential expression analysis) using
the [Seurat](https://satijalab.org/seurat/index.html) framework. Then,
we **simplify** the same dataset by computing *metacells* (i.e.,
grouping transcriptionally highly similar single cells into metacells).
For this, we will use a method developed in our group called
[SuperCell](https://github.com/GfellerLab/SuperCell). We will also
provide some hints on how *metacells* can be computed using alternative
approaches including [MetaCell](https://github.com/tanaylab/metacell),
[Metacell2.0](https://metacells.readthedocs.io/en/latest/readme.html),
and [SEACell](https://github.com/dpeerlab/SEACells). We will then run
**‘a standard scRNA-seq data analysis pipeline’** adjusted to the
metacell data and compare the results obtained at the single-cell and
the metacell levels.

Then, the participant can use the code we provide to build metacells and
run some basic analysis for their datasets

------------------------------------------------------------------------

**Some options:**

-   Alternatively, we can run an adjusted (ie., sample-weighted)
    pipeline for metacells or/and a standard (ie, Seurat).

-   Subsampling at the same graining level -&gt; compare the results to
    those obtained at the single-cell level

-   We could provide an example of how metacell can be used for
    RNA-velocity (on another dataset)

-   How metacells can be used for data integration (if we have nice
    examples)

-   when and if analyzing their own data, make sure to introduce
    `cell.annotation` and `cell.split.condition` arguments of
    `SCimplify()` to avoid mixing of annotated cell types / conditions
    within metacells.

------------------------------------------------------------------------

TO DO:
======

-   <s>make `supercell_DimPlot()`</s> push `supercell_DimPlot()`
-   match colors for cell lines and clusters
-   write ‘standard (not sample-weighted)’ analysis on metacells

-   code for MC2 and SEAcells computation (to provide a link), save and
    provide results of metcell membership computed with MC2 and SEACell

-   dataset for demo of gene-gene correlation/RNA velocity/ Data
    integratoin for more advanced steps of the analysis
-   ATAC-seq example ?

<!-- -->

    # make a data library (cell lines or Zilionis)
    library(SuperCell)
    library(Seurat)
    library(dplyr)

    proj.name    <- 'cell_lines'
    data.folder  <- file.path("..", "data", proj.name)

    # load single-cell (sc) count matrix and cell metadata 
    sc.counts <- readRDS(file.path(data.folder, "sc_counts_filtered.Rds"))
    sc.meta   <- readRDS(file.path(data.folder, "sc_meta_filtered.Rds"))

    # Make sure metadata and count matrix have the same cells in the same order
    if(!identical(rownames(sc.meta), colnames(sc.counts))){
      stop("Metadata (`sc.meta`) does not correspond to the count matrix (`sc.counts`)")
    }

Single-cell level
=================

Standard downstream analysis
----------------------------

Run a brief analysis at the single-cell level, lets’s use the common
[Seurat](https://satijalab.org/seurat/index.html) pipeline

<span style="color: green;"> depending on the participants’ experience
with Seurat, we can either separately run these steps of the
pre-processing or ‘skip it’ with ‘one-line’ command. For the moment, it
is more step-by -step pre-processing, that can be replaced with :
</span>
`sc <-  NormalizeData(sc, verbose=FALSE) %>% FindVariableFeatures(selection.method = "disp", nfeatures = 1000, verbose=FALSE) %>% ScaleData(verbose=FALSE) %>% RunPCA(verbose=FALSE)`

    set.seed(12345)
    sc <- CreateSeuratObject(counts = sc.counts, project = proj.name, meta.data = sc.meta)
    sc

    ## An object of class Seurat 
    ## 11786 features across 3822 samples within 1 assay 
    ## Active assay: RNA (11786 features, 0 variable features)

### Pre-processing

#### Data normalization

    sc <- NormalizeData(sc, verbose=FALSE)

    sc <- FindVariableFeatures(
      sc, 
      selection.method = "disp", # "vst" is default
      nfeatures = 1000,
      verbose=FALSE
      )

    hvg <- VariableFeatures(sc, verbose=FALSE)

    # Plot variable features 
    plot1 <- VariableFeaturePlot(sc)
    LabelPoints(plot = plot1, points = hvg[1:20], repel = TRUE)

![](Cell_lines_files/figure-markdown_strict/Normalize%20data%20and%20compute%20a%20set%20of%20highly%20variable%20genes-1.png)

#### Scaling and dimensionality redution

    sc <- ScaleData(sc, verbose=FALSE)
    sc <- RunPCA(sc, verbose=FALSE)

    # Plot PCA (2D representation of scRNA-seq data) colored by cell line
    DimPlot(sc, reduction = "pca", group.by = "cell_line")

![](Cell_lines_files/figure-markdown_strict/unnamed-chunk-2-1.png)

### UMAP (non-linear dimensionality reduction)

    sc <- RunUMAP(sc,  dims = 1:10)

    # Plot UMAP (2D representation of scRNA-seq data) colored by cell line
    DimPlot(sc, reduction = "umap", group.by = "cell_line")

![](Cell_lines_files/figure-markdown_strict/UMAP-1.png)

### Clustering

    sc <- FindNeighbors(sc, dims = 1:10)
    sc <- FindClusters(sc, resolution = 0.05)

    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 3822
    ## Number of edges: 121361
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.9890
    ## Number of communities: 5
    ## Elapsed time: 0 seconds

    # As it is a toy example with well defined cell types (i.e., cell lines), unsupervised clustering fully recapitulates cell line annotation 
    table(sc@active.ident, sc$cell_line)

    ##    
    ##     A549 H1975 H2228 H838 HCC827
    ##   0 1237     0     0    0      0
    ##   1    0     0     0  841      0
    ##   2    0     0   744    0      0
    ##   3    0     0     0    0    571
    ##   4    0   429     0    0      0

    DimPlot(sc, reduction = "umap", group.by = "ident")

![](Cell_lines_files/figure-markdown_strict/unnamed-chunk-3-1.png)

### Differential expression analysis

#### Find Markers of cell lines

    # Set idents to cell lines (as clusters are the same as cell lines)
    Idents(sc) <- "cell_line"

    # Compute upregulated genes in each cell line (versus other cells)
    sc.all.markers <-  FindAllMarkers(sc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "t")
    saveRDS(sc.all.markers, file = file.path(data.folder, "output", "sc_all_markers.Rds"))

    # Load markers (backup)
    #sc.all.markers <- readRDS(file = file.path(data.folder, "output", "sc_all_markers.Rds"))

    # Top markers (select top markers of each cell line)
    sc.top.markers <- sc.all.markers %>%
       group_by(cluster) %>%
        slice_max(n = 2, order_by = avg_log2FC)

    sc.top.markers

    ## # A tibble: 10 × 7
    ## # Groups:   cluster [5]
    ##        p_val avg_log2FC pct.1 pct.2 p_val_adj cluster gene   
    ##        <dbl>      <dbl> <dbl> <dbl>     <dbl> <fct>   <chr>  
    ##  1 0               4.17 1     0.996 0         HCC827  SEC61G 
    ##  2 0               3.05 0.998 0.941 0         HCC827  CDK4   
    ##  3 0               4.20 1     0.245 0         H838    GAGE12D
    ##  4 0               4.01 1     0.176 0         H838    GAGE12E
    ##  5 1.05e-180       3.26 0.991 0.323 1.24e-176 H1975   MT1E   
    ##  6 1.42e-142       3.07 0.963 0.083 1.68e-138 H1975   DHRS2  
    ##  7 0               4.08 0.997 0.371 0         H2228   SAA1   
    ##  8 0               3.94 0.999 0.625 0         H2228   LCN2   
    ##  9 0               5.56 0.998 0.775 0         A549    KRT81  
    ## 10 0               5.15 1     0.675 0         A549    AKR1B10

#### Plot the expression of some found markers

    VlnPlot(sc, features = sc.top.markers$gene[c(seq(1, 9, 2), seq(2, 10, 2))], ncol = 5, pt.size = 0.0)

![](Cell_lines_files/figure-markdown_strict/unnamed-chunk-5-1.png)

Data simplification (coarse-graining) – Construction of *metacells*
===================================================================

Here we compute metacells using our method called
[SuperCell](https://github.com/GfellerLab/SuperCell), but equally,
metacells can be computed with
[Metacell](https://github.com/tanaylab/metacell),
[Metacell2.0](https://metacells.readthedocs.io/en/latest/readme.html) or
[SEACell](https://github.com/dpeerlab/SEACells) algorithms and we will
see some examples below.

    gamma <- 20 # Graining level

    # Compute metacells using SuperCell package
    MC <- SCimplify(
      X = GetAssayData(sc), # single-cell log-normalized gene expression data
      genes.use = hvg, 
      gamma = gamma,
      n.pc = 10
    )

    # Compute gene expression of metacells by simply averaging gene expression within each metacell
    MC.ge <- supercell_GE(
      ge = GetAssayData(sc),
      groups = MC$membership
    )

    # Alternatively, counts can be averaged (summed up) followed by a lognormalization step (this approach is used in the MetaCell and SEACell algorithms)
    if(0){
      MC.counts <- supercell_GE(
        ge = GetAssayData(sc, slot = "counts"),
        mode = "sum", # summing counts instead of the default averaging
        groups = MC$membership
      )
      
      MC.ge <- Seurat::LogNormalize(MC.counts, verbose = FALSE)
    }

Downstream analysis of metacells
--------------------------------

There are two options to perform the downstream analysis:

-   **sample-weighted** when we account for a metacell size at each etep
    of the analysis
-   **standard** when we treat metacells as single-cell and apply a
    standard pipeline

Sample-weighted downstream analysis of *metacells*
--------------------------------------------------

For the sample-weighted analysis, we use a pipeline avalable with our
[SuperCell]() package.

### Pre-processing

#### Transfer metadata (annotate metacell to a certain cell line)

Since the cell line information is available in this dataset, we can
annotate metacells to a certain cell line. Each metacell is annotated to
the most abundant cell type in it. This also allows us to compute
metacell *purity*, that is defined as a proportion of the most abundant
cell type (use `method = "max_proportion"`) or as Shannon entropy (use
`method = "entropy"`).

    MC$cell_line <- supercell_assign(
      cluster = sc.meta$cell_line,          # single-cell assigment to cell lines 
      supercell_membership = MC$membership  # single-cell assignment to metacells
    )

    method_purity <- c("max_proportion", "entropy")[1]
    MC.purity <- supercell_purity(
      clusters = sc.meta$cell_line,
      supercell_membership = MC$membership, 
      method = method_purity
    )

    # Metacell purity distribution
    summary(MC.purity)

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##       1       1       1       1       1       1

    #hist(MC.purity, main = paste0("Purity of metacells \nin terms of cell line composition (", method_purity,")"))

#### Visualize data using metcall network

    supercell_plot(
      MC$graph.supercells, 
      group = MC$cell_line, 
      seed = 1, 
      alpha = -pi/2,
      main  = "Metacells colored by cell line assignment"
    )

![](Cell_lines_files/figure-markdown_strict/metacell%20plot%20metacell%20network%20color%20cell%20line-1.png)

#### Dimensionality reduction (sample-weighted PCA)

    MC$PCA <- supercell_prcomp(
      Matrix::t(MC.ge),
      genes.use = MC$genes.use,  # or a new set of HVG can be computed
      supercell_size = MC$supercell_size, # provide this parameter to run sample-weighted version of PCA,
      k = 20
    )

    supercell_DimPlot(
      MC, 
      groups = MC$cell_line,
      dim.name = "PCA", 
      title = paste0("PCA of metacells colored by cell line assignment")
    )

![](Cell_lines_files/figure-markdown_strict/metacell%20PCA-1.png)

#### Clustering (sample-weighted hclust)

Sample-weighted clustering computed with the hierarchical clustering,
that may accounts for sample weights

    ## compute distance
    D                <- dist(MC$PCA$x)

    ## cluster metacells
    MC$clustering    <- supercell_cluster(D = D, k = 5, supercell_size = MC$supercell_size)

    # Plot clustering result
    supercell_DimPlot(
      MC, 
      groups = factor(MC$clustering$clustering),
      dim.name = "PCA", 
      title = paste0("PCA of metacells colored by metacell clustering")
    )

![](Cell_lines_files/figure-markdown_strict/metacell%20clustering-1.png)

    table(MC$cell_line, MC$clustering$clustering)

    ##         
    ##           1  2  3  4  5
    ##   A549   51  0  0  0  0
    ##   H1975   0  0  0  0 27
    ##   H2228   0  0 40  0  0
    ##   H838    0  0  0 45  0
    ##   HCC827  0 28  0  0  0

### Differential expression analysis in metacells

    # Compute upregulated genes in each cell line (versus other cells)
    MC.all.markers <- supercell_FindAllMarkers(
      ge = MC.ge, 
      clusters = MC$cell_line, 
      supercell_size = MC$supercell_size,
      only.pos = TRUE, 
      min.pct = 0.25, 
      logfc.threshold = 0.25
    )

    saveRDS(MC.all.markers, file = file.path(data.folder, "output",  paste0("MC_gamma_", gamma, "_all_markers.Rds")))

    # Load markers (backup)
    #MC.all.markers <- readRDS(file = file.path(data.folder, "output", "paste0("MC_gamma_", gamma, "_all_markers.Rds")))

    MC.all.markers.df <- data.frame()
    for(cl in names(MC.all.markers)){
      cur <- MC.all.markers[[cl]]
      cur$cluster <- cl
      cur$gene <- rownames(cur)
      cur$avg_log2FC <- cur$logFC
      MC.all.markers.df <- rbind(MC.all.markers.df, cur)
    }


    # Top markers (select top markers of each cell line)
    MC.top.markers <- MC.all.markers.df %>%
       group_by(cluster) %>%
        slice_max(n = 2, order_by = avg_log2FC)

#### Plot the expression of some found markers (in metacells)

    supercell_VlnPlot(
      ge = MC.ge, 
      supercell_size =MC$supercell_size, 
      clusters = MC$cell_line,
      features = MC.top.markers$gene[c(seq(1, 9, 2), seq(2, 10, 2))],
      ncol = 5)

![](Cell_lines_files/figure-markdown_strict/unnamed-chunk-6-1.png)

Standard downstream analysis
----------------------------

For the standard downstream analysis, we can use the well-stablished
[Seurat]() pipeline
