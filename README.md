
<!-- README.md is generated from README.Rmd. Please edit that file -->

# fastBCR

### a heuristic method for fast BCR clonal family inference from large-scale AIRR-seq data

<!-- badges: start -->
<!-- badges: end -->

The goal of fastBCR is to rapid inferencing B cell clonal families from
massive BCR heavy chain sequences. To facilitate the post-clustering
analysis, we also provide a series of functional modules in fastBCR,
including multiple sequence alignment (MSA), phylogenetic tree
construction, somatic hypermutation (SHM), and class switch analysis
(CSR), which can be useful for BCR repertoire analysis and antibody
discovery.

### REFERENCE

Kaixuan Wang, Xihao Hu, Jian Zhang\*. **Fast Clonal Family Inference from Large-scale 
B Cell Repertoire Sequencing Data.** (Available at SSRN: https://ssrn.com/abstract=4431127 or http://dx.doi.org/10.2139/ssrn.4431127)

## Installation

Before installing fastBCR, you need to download the dependency packages
msa, ggtree and ggmsa using Bioconductor. To install these packages, start R
(version “4.1.0”) and enter:

``` r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("proj4", "msa","ggtree","ggmsa"))
```

Or you can install all the required R packages through the _requirement.R_ file by entering:

``` r
source("requirement.R")
```

Now you can install the development version of fastBCR like so:

``` r
install.packages("devtools")
devtools::install_github("ZhangLabTJU/fastBCR")
library(fastBCR)
```

## Sample dataset

A real example AIRR Rearrangement dataset (‘test_data’) is included in
the fastBCR package. The dataset consists of BCR sequencing data from
peripheral blood samples of a COVID-19 patient (Galson et al., 2020).
Inferencing clonal families requires the following columns to be present
in the table:

``` r
## "sequence_id"
## "v_call"
## "j_call"
## "junction_aa"
## "junction" (Optional. Needed for phylogenetic tree construction.)
## "c_call" (Optional. Needed for isotypes related analyis.)
```

## Example

fastBCR is an automatic BCR clonal family inference method, which also
incorporates multiple functional modules for downstream analyses.
fastBCR is composed of three parts: **Data Processing**, **Clonal Family Inferencing**, and
**Downstream Analysis**. In addition, fastBCR can be combined with annotation 
software (_e.g._ IgBlast) for B-cell **Clonal Family Simulation**.

<img src="man/figures/README-overview-1.png" width="100%" />

### Data Processing

Data Processing is the first step of fastBCR to process raw data. The
input of the function needs to meet the AIRR standard format. Only
productive sequences whose junction amino acid lengths between 9 and 26
are reserved. Sequences with the same ‘v_call’, ‘j_call’ and
‘junction_aa’ are considered identical and deduplicated.

``` r
#> Loading required package: ggplot2
#> Loading required package: network
#> 
#> 'network' 1.18.1 (2023-01-24), part of the Statnet Project
#> * 'news(package="network")' for changes since last version
#> * 'citation("network")' for citation information
#> * 'https://statnet.org' for help, support, and other information
#> Loading required package: ggmsa
#> Registered S3 methods overwritten by 'ggalt':
#>   method                  from   
#>   grid.draw.absoluteGrob  ggplot2
#>   grobHeight.absoluteGrob ggplot2
#>   grobWidth.absoluteGrob  ggplot2
#>   grobX.absoluteGrob      ggplot2
#>   grobY.absoluteGrob      ggplot2
#> ggmsa v1.4.0  Document: http://yulab-smu.top/ggmsa/
#> 
#> If you use ggmsa in published research, please cite:
#> L Zhou, T Feng, S Xu, F Gao, TT Lam, Q Wang, T Wu, H Huang, L Zhan, L Li, Y Guan, Z Dai*, G Yu* ggmsa: a visual exploration tool for multiple sequence alignment and associated data. Briefings in Bioinformatics. DOI:10.1093/bib/bbac222
#> Loading required package: msa
#> Loading required package: Biostrings
#> Loading required package: BiocGenerics
#> 
#> Attaching package: 'BiocGenerics'
#> The following objects are masked from 'package:stats':
#> 
#>     IQR, mad, sd, var, xtabs
#> The following objects are masked from 'package:base':
#> 
#>     anyDuplicated, aperm, append, as.data.frame, basename, cbind,
#>     colnames, dirname, do.call, duplicated, eval, evalq, Filter, Find,
#>     get, grep, grepl, intersect, is.unsorted, lapply, Map, mapply,
#>     match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
#>     Position, rank, rbind, Reduce, rownames, sapply, setdiff, sort,
#>     table, tapply, union, unique, unsplit, which.max, which.min
#> Loading required package: S4Vectors
#> Loading required package: stats4
#> 
#> Attaching package: 'S4Vectors'
#> The following objects are masked from 'package:base':
#> 
#>     expand.grid, I, unname
#> Loading required package: IRanges
#> Loading required package: XVector
#> Loading required package: GenomeInfoDb
#> 
#> Attaching package: 'Biostrings'
#> The following object is masked from 'package:base':
#> 
#>     strsplit
data("test_data")
input <- data.pro(test_data)
```

### BCR Clustering

BCR Clustering is the core step of fastBCR to infer clonal families from
processed data. BCR Clustering consists of two steps, k-mer-based
pre-clustering and optimized clustering, to infer clonal families fast
and accurately.

``` r
data("input")
bcr_clusters <- BCR.cluster(input) # A few parameters are optional. See the function introduction (?BCR.cluster) for details.
#> use default substitution matrix
#> use default substitution matrix
#> ...
#> use default substitution matrix
```

### Downstream Analysis

Downstream Analysis contains a variety of functional modules, such as
multiple sequence alignment (MSA), phylogenetic tree construction, and
SHM/CSR analysis. You can choose the plot function for downstream
analysis according to your needs.

#### Pieplot of V/J usage frequency

``` r
data("input")
pie.freq(input, "v_call")
```
<img src="man/figures/README-pie.freq-1.png" width="50%" />

#### Histogram of junction length

``` r
data("input")
junc.len(input, "AA")
```
<img src="man/figures/README-junc.len-1.png" width="50%" />

#### Plot of phylogenetic tree and multiple sequences alignment

``` r
data(bcr_clusters)
msa.tree(bcr_clusters, 20)
#> use default substitution matrix
#> Parameter setting = useAbundance:  True ; revision:  True ; trim: True
#> done
```
<img src="man/figures/README-msa.tree-1.png" width="80%" />

#### Reconstructing B cell lineage trees with minimum spanning tree and genotype abundances (Optional)

ClonalTree is a new algorithm to reconstruct BCR lineage trees that incorporates genotype abundance into a minimum spanning tree to infer maximum parsimony trees.(URL: https://github.com/julibinho/ClonalTree)

##### Requirements 

  * numpy :
      ```
      $ conda install numpy
      ```
      or 
      ```
      $ pip install numpy
      ```

  * Biopython
      ```
      $ pip install biopython
      ```

  * ete3 :
      ```
      $ pip install ete3
      ```

  * networkx (for the evaluation) :
      ```
      $ pip install networkx
      ```

##### Using ClonalTree      

``` r
data(bcr_clusters)
Clonal.tree(bcr_clusters, 1)
#> use default substitution matrix
#> Parameter setting = useAbundance:  True ; revision:  True ; trim: True
#> done
```
<img src="man/figures/README-Clonal.tree-1.png" width="70%" />

#### Visualization of junction_aa sequences within a cluster

``` r
data(bcr_clusters)
msa.logo(bcr_clusters, 20)
#> use default substitution matrix
#> Coordinate system already present. Adding new coordinate system, which will
#> replace the existing one.
```

<img src="man/figures/README-msa.logo-1.png" width="30%" />

#### Boxplot of SHM ratios in four isotypes in different samples

``` r
data("SHM_ratio_1")
data("SHM_ratio_2")
data("SHM_ratio_3")
data("SHM_ratio_4")
data("SHM_ratio_5")
df <- data.frame(
  Isotypes = factor(c("IGHD", "IGHM", "IGHA", "IGHG"),
    levels = c("IGHD", "IGHM", "IGHA", "IGHG")
  ),
  Sample1 = SHM_ratio_1,
  Sample2 = SHM_ratio_2,
  Sample3 = SHM_ratio_3,
  Sample4 = SHM_ratio_4,
  Sample5 = SHM_ratio_5
)
SHM.sample(df)
#> Using Isotypes as id variables
```

<img src="man/figures/README-SHM.sample-1.png" width="50%" />

#### Scatter diagram of cluster size (x axis), cluster number (y axis) and SHM ratio (point size) of all clusters in a sample

``` r
data("bcr_clusters")
SHM.cluster(bcr_clusters)
#> use default substitution matrix
#> use default substitution matrix
#> ...
#> use default substitution matrix
```

<img src="man/figures/README-SHM.cluster-1.png" width="50%" />

#### Network of CSR within a sample

``` r
data("bcr_clusters")
CSR.sample(bcr_clusters)
#> Registered S3 method overwritten by 'GGally':
#>   method from   
#>   +.gg   ggplot2
```

<img src="man/figures/README-CSR.sample-1.png" width="50%" />

#### Network of CSR within a cluster

``` r
data("bcr_clusters")
CSR.cluster(bcr_clusters, 1)
```

<img src="man/figures/README-CSR.cluster-1.png" width="50%" />

### Clonal Family Simulation

fastBCR can simulate the generation of B cell clonal families to evaluate the 
performance of different clonal family inference methods. Specifically, it 
begins by generating an ancestor cell through V(D)J random recombinant and 
simulated the process of antigen activation that led to multiple rounds of 
expansion, mutation in the junction region and elimination, ultimately 
resulting in the formation of a B cell clonal family. This process needs to be 
combined with sequence annotation software, and here we recommend using IgBlast.
(URL: https://changeo.readthedocs.io/en/stable/examples/igblast.html)

<img src="man/figures/README-Simulation-1.png" width="100%" />

``` r
germline2fas(100, filename = "Simulation/Germline.fasta")
# Annotation
germline_data = read.table('Simulation/Germline_igblast_db-pass_parse-select.tsv', header = T, sep = "\t")
CF2fas(germline_data, CF_n = 10, mut_ratio = 0.001, filename = 'Simulation/10_0.001.fasta')
# Annotation
```
