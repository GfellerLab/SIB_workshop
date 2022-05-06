Simplifying large and complex single-cell RNA-Seq data with metacells
================

## Metacell concept

The exponential scaling of scRNA-seq data represents an important hurdle
for downstream analyses. One of the solutions to simplify large-scale
and noisy scRNA-seq data is to merge transcriptionally similar cells
into *metacells*. The concept was firsly introduced in [*Baran et
al., 2019*](https://doi.org/10.1186/s13059-019-1812-2) and in [*Iacono
et al., 2018*](https://doi:10.1101/gr.230771.117). Other methods to
build *metacells* are presented in [*Ben-Kiki et
al. 2022*](https://doi.org/10.1186/s13059-022-02667-1), [*Bilous et
al., 2022*](https://www.biorxiv.org/content/10.1101/2021.06.07.447430v2)
and [*Persad et
al., 2022*](https://www.biorxiv.org/content/10.1101/2022.04.02.486748v1).
Despite of some difference in implementation, all the methods are
network-based and can be summarized as follows:

1.  A single-cell network is computed based on cell-to-cell similarity
    (in transcriptomic space)
2.  Highly similar cells are identified as those forming dense reagions
    in the single-cell network
3.  Transcriptomic information within identified metacells is merges
    (i.e., summed up or averaged)
4.  Metacell data are used for the downtream analyses instead of
    large-scale single-cell data

<img src="plots/1.png" width="4932" />

The simplifiation is performed at a specific **graining level**
(\(\gamma\)), which we define as a ration between the number of single
cells in the inintial data and the number of metacells. Depening on the
simplification method, the grainin level can be either specified by a
user (in [bigSCale](https://github.com/iaconogi/bigSCale2),
[SuperCell](https://github.com/GfellerLab/SuperCell) and
[SEACells](https://github.com/dpeerlab/SEACells)) or provided by the
simlification algorithm (in
[Metacell](https://github.com/tanaylab/metacell) and
[Metacell-2](https://github.com/tanaylab/metacells)).

## Tutorial structure

  - An
    [axample](https://github.com/GfellerLab/SIB_workshop/blob/main/workbooks/Cell_lines.md)
    of the simplification of a small dataset and the comparison between
    single-cell and metacell downstream analysis results (basic steps
    including visualization, clusering and simensionality reduction).

  - [Construction of metacells]() using [Metacell-2]() and [SEACell]()
    algorithms.

  - An [axample on more interesting data]() (e.g., PBMC) and maybe more
    advanced steps for those who do not have their data to analyze. For
    those, who have own data, this code can be used for the analysis.

  - Using metacells for [data integration]() (would be nice to show an
    example where integration is possible at the metacell, but not at
    the single-cell level)

  - Using [metacells for RNA velosity]().

  - Using [metacells for scATAC-seq data analysis]() (if we have an
    example)
