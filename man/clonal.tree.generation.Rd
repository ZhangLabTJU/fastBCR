% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Downstream_Analysis.R
\name{clonal.tree.generation}
\alias{clonal.tree.generation}
\title{Function: Reconstructing a B cell lineage tree with minimum spanning tree and genotype abundances using ClonalTree}
\usage{
clonal.tree.generation(bcr_clusters, index, raw_data, python_path)
}
\arguments{
\item{bcr_clusters}{Clonal families inferred by fastBCR}

\item{index}{Index of cluster}

\item{python_path}{The absolute path of the Python interpreter
#' @param raw_data The raw data which the clonal families inferred from
fastBCR will retrieve all the DNA sequences, which can be multiple sequences due to the degeneracy of codons, that correspond to the amino acid sequence of each clonotype from the raw data}
}
\value{
ClonalTree returns two files in the 'ClonalTree/Examples/output' folder
ClonalFamily_index.nk: the reconstructed BCR lineage tree in newick format
ClonalFamily_index.nk.csv: a table in csv format, containing the parent relationship and cost
}
\description{
The junction of BCR sequences within a cluster will be written as 'ClonalFamily_index.fasta' in 'ClonalTree/Examples/input' folder
}
