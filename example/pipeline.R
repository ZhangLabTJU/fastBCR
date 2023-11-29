library(fastBCR)

# You should set "fastBCR" folder as working directory before running pipeline.

### Fast BCR clonal family inference
## 1. Data loading
# Path to the "COVID"/"HC" folder in the "example" folder in fastBCR package.
COVID_folder_path <- "example/COVID"
HC_folder_path <- "example/HC"
# Load files from "COVID_folder_path"/"HC_folder_path" into list.
# The storage format of data can be in "csv", "tsv", or "Rdata" format.
# The compressed files in the above storage format in ‘7zip’, ‘cab’, ‘cpio’, ‘iso9660’, ‘lha’, ‘mtree’, ‘shar’, ‘rar’, ‘raw’, ‘tar’, ‘xar’, ‘zip’, ‘warc’ format can also be read in.
COVID_raw_data_list <- data.load(folder_path = COVID_folder_path, storage_format = "csv")
HC_raw_data_list <- data.load(folder_path = HC_folder_path, storage_format = "csv")
# # Or you can load one file at a time.
# COVID_file_path <- "~/Documents/Rpackage/fastBCR/example/COVID/COVID_01.zip"
# COVID_load_file <- data.load(folder_path = COVID_file_path, storage_format = "csv")

## 2. Data preprocessing
# Preprocessing of raw data to meet the input requirements for clonal family inference.
# The input of the function needs to meet the AIRR standard format (containing at least "sequence_id", "v_call", "j_call", and "junction_aa" information).
# Only productive sequences whose junction amino acid lengths between 9 and 26 are reserved.
# Sequences with the same "v_call", "j_call" and "junction_aa" are considered identical and deduplicated.
COVID_pro_data_list <- data.preprocess(data_list = COVID_raw_data_list)
HC_pro_data_list <- data.preprocess(data_list = HC_raw_data_list)

## 3. BCR clonal family inference
# Fast clonal family inference from preprocessed data.
# The "cluster_thre" parameter represents minimal clustering criteria (minimum number of sequences needed to form a cluster) and defaults to 3.
# For high efficiency, the "cluster_thre" is increased by 1 for every 100,000 entries of input data.
# The "overlap_thre" parameter represents overlap coefficient threshold for merging two clusters, selectable range (0,1) and defaults to 0.1.
# Lower "overlap_thre" may lead to overclustering while higher thresholds may lead to the split of clonal families.
# The "consensus_thre" parameter represents the consensus score threshold for filtering candidates and defaults to 0.8.
# A higher "consensus_thre" means stricter inference of the cluster.
COVID_clusters_list <- data.BCR.clusters(pro_data_list = COVID_pro_data_list, cluster_thre = 3, overlap_thre = 0.1, consensus_thre = 0.8)
HC_clusters_list <- data.BCR.clusters(pro_data_list = HC_pro_data_list, cluster_thre = 3, overlap_thre = 0.1, consensus_thre = 0.8)

## 4. Classification of clustered and unclustered sequences
# Merge all the clustered sequences in each sample into "clustered_seqs".
# Merge all the unclustered sequences in each sample into "unclustered_seqs".
COVID_seqs_list <- Clustered.seqs(pro_data_list = COVID_pro_data_list, clusters_list = COVID_clusters_list)
HC_seqs_list <- Clustered.seqs(pro_data_list = HC_pro_data_list, clusters_list = HC_clusters_list)



### Downstream analysis
## 1. Single cluster analysis
## 1.1 Conserved motifs distribution
# Visualization of multiple sequence alignment (MSA) of junction sequences within a cluster.
# The parameter "index" allows you to choose a cluster for visualization.
# The parameter "type" can be "DNA" for deoxyribonucleic acid or "AA" for amino acid.
COVID_01_clusters = COVID_clusters_list$COVID_01
HC_01_clusters = HC_clusters_list$HC_01
msa.plot(bcr_clusters = COVID_01_clusters, index = 200, type = "DNA")
msa.plot(bcr_clusters = COVID_01_clusters, index = 200, type = "AA")

## 1.2 Sequence logo
# Visualization of sequence logo of junction sequences within a cluster.
seqlogo.plot(bcr_clusters = COVID_01_clusters, index = 200, type = "DNA")
seqlogo.plot(bcr_clusters = COVID_01_clusters, index = 200, type = "AA")

## 1.3 Evolutionary tree
# Reconstruct a B cell lineage tree with minimum spanning tree and genotype abundances using ClonalTree.
# The junction of BCR sequences within a cluster will be written as "ClonalFamily_index.fasta" in "ClonalTree/Examples/input" folder.
# ClonalTree uses Python for compilation, so it needs the absolute path of the Python interpreter in the "python_path" parameter.
# ClonalTree returns two files in the "ClonalTree/Examples/output" folder.
# ClonalFamily_index.nk: the reconstructed BCR lineage tree in newick format.
# ClonalFamily_index.nk.csv: a table in csv format, containing the parent relationship and cost.
clonal.tree.generation(bcr_clusters = COVID_01_clusters, index = 200, python_path = "/Users/wangkaixuan/anaconda3/bin/python")
# You can use the clonal.tree.plot() function to visualize the evolutionary tree in R.
# Orange points represent nodes and blue points represent tips.
# Abbreviation of "sequence_id" (last 5 characters) are shown.
# The x-axis shows the absolute genetic distance.
clonal.tree.plot(nk_path = "ClonalTree/Examples/output/ClonalFamily_200.abRT.nk")

## 1.4 Class switch recombination (CSR) network
# Visualization of isotype co-occurrence within clusters within a cluster
# Circle size represents the number of sequences carrying a given isotype.
# Lines connecting two circles indicate the enrichment level of observing switches in the two corresponding immunoglobulin subclasses.
# The enrichment level is the ratio of observed and expected switches if immunoglobulin isotypes are assumed to be independently distributed among cluster.
# Matrix of values of connected edges between clustered sequences in different isotypes is printed.
CSR.cluster.plot(bcr_clusters = COVID_01_clusters, index = 50)


## 2. Single sample/group analysis
## 2.1 Summary of clusters from a sample
# Summarize the number of clusters, the average size of clusters and the proportion of clustered sequences.
COVID_clusters_summary = Clusters.summary(pro_data_list = COVID_pro_data_list, clusters_list = COVID_clusters_list)
HC_clusters_summary = Clusters.summary(pro_data_list = HC_pro_data_list, clusters_list = HC_clusters_list)
COVID_01_summary = COVID_clusters_summary$COVID_01
print(COVID_01_summary)
HC_01_summary = HC_clusters_summary$HC_01
print(HC_01_summary)

## 2.2 Visualization of clusters from a sample
# Point diagram showing clusters from a sample where a circle represents a cluster.
# The size and color of the circle represents the size of the cluster.
# The parameter "index" represents the index of the sample for visualization.
Clusters.visualization(pro_data_list = COVID_pro_data_list, clusters_list = COVID_clusters_list, index = 1)
Clusters.visualization(pro_data_list = HC_pro_data_list, clusters_list = HC_clusters_list, index = 1)

## 2.3 V/J gene usage
## 2.3.1 Pie chart: V/J gene
# Pie chart showing the gene usage frequency of clustered sequences.
# The top ten most frequent genes are shown, and the rest are represented by "others".
# The parameter "colname" can be "v_call" for V gene or "j_call" for J gene.
# (1) single sample
COVID_01_clustered_seqs = COVID_seqs_list$clustered_seqs$COVID_01
COVID_01_unclustered_seqs = COVID_seqs_list$unclustered_seqs$COVID_01
HC_01_clustered_seqs = HC_seqs_list$clustered_seqs$HC_01
pie.freq.plot(clustered_seqs = COVID_01_clustered_seqs, colname = "v_call")
pie.freq.plot(clustered_seqs = HC_01_clustered_seqs, colname = "v_call")
# (2) multiple samples in a group.
# All the clustered/unclustered sequences from multiple samples in a group can be merged.
COVID_all_clustered_seqs = NULL
for(i in 1:length(COVID_seqs_list$clustered_seqs)){
  COVID_all_clustered_seqs = rbind(COVID_all_clustered_seqs, COVID_seqs_list$clustered_seqs[[i]])
}
COVID_all_unclustered_seqs = NULL
for(i in 1:length(COVID_seqs_list$unclustered_seqs)){
  COVID_all_unclustered_seqs = rbind(COVID_all_unclustered_seqs, COVID_seqs_list$unclustered_seqs[[i]])
}
HC_all_clustered_seqs = NULL
for(i in 1:length(HC_seqs_list$clustered_seqs)){
  HC_all_clustered_seqs = rbind(HC_all_clustered_seqs, HC_seqs_list$clustered_seqs[[i]])
}
HC_all_unclustered_seqs = NULL
for(i in 1:length(HC_seqs_list$unclustered_seqs)){
  HC_all_unclustered_seqs = rbind(HC_all_unclustered_seqs, HC_seqs_list$unclustered_seqs[[i]])
}
pie.freq.plot(clustered_seqs = COVID_all_clustered_seqs, colname = "v_call")
pie.freq.plot(clustered_seqs = HC_all_clustered_seqs, colname = "v_call")
## 2.3.2 Heatmap: V-J gene pair
# Heatmap showing the V-J gene pair frequency of clustered sequences
# (1) single sample.
vjpair.sample.plot(clustered_seqs = COVID_01_clustered_seqs)
# (2) multiple samples in a group.
vjpair.sample.plot(clustered_seqs = COVID_all_clustered_seqs)
vjpair.sample.plot(clustered_seqs = HC_all_clustered_seqs)

## 2.4 Junction length distribution
# Histogram and density plot showing the length distribution of junction amino acid of clustered sequences.
# (1) single sample
len.sample.plot(clustered_seqs = COVID_01_clustered_seqs)
# (2) multiple samples in a group
len.sample.plot(clustered_seqs = COVID_all_clustered_seqs)
# Density ridges showing the length distributions of junction amino acid of clustered sequences and unclustered sequences.
# Statistical comparisons are carried out by the two-sided Wilcoxon rank-sum test.
# (1) single sample
len.clustered.plot(clustered_seqs = COVID_01_clustered_seqs,
                   unclustered_seqs = COVID_01_unclustered_seqs)
# (2) multiple samples in a group
len.clustered.plot(clustered_seqs = COVID_all_clustered_seqs,
                   unclustered_seqs = COVID_all_unclustered_seqs)

## 2.5 Class switch recombination (CSR) network
# Visualization of isotype co-occurrence within clusters from a sample.
# Circle size represents the number of sequences carrying a given isotype.
# Lines connecting two circles indicate the enrichment level of observing switches in the two corresponding immunoglobulin subclasses.
# The enrichment level is the ratio of observed and expected switches if immunoglobulin isotypes are assumed to be independently distributed among cluster.
# Matrix of values of connected edges between clustered sequences in different isotypes is printed.
CSR.sample.plot(bcr_clusters = COVID_01_clusters)
CSR.sample.plot(bcr_clusters = HC_01_clusters)

### 2.6 Neutralizing antibody (NAb) sequence query
## 2.6.1 Load the public antibody database to get the information of neutralizing antibody (NAb) sequences.
# Here is an example of the Coronavirus Antibody Database (CoV-AbDab).
CoV_AbDab = read.csv("example/CoV-AbDab_130623.csv")
## 2.6.2 NAb query
# Query the corresponding sequence from the public antibody database.
# The parameter "method" represents the CDRH3 matching method.
# It can be "NA" for perfect match, "hamming" for hamming distance or "lv" for Levenshtein distance. Defaults to "NA".
# The parameter "species" can be "Mouse" or "Human". Defaults to "Human".
# The parameter "maxDist" represents the maximum distance allowed for matching when the argument "method" is "hamming" or "lv". Defaults to "NA".
# example to perfect match
human_perfect_match <- NAb.query(bcr_clusters = COVID_01_clusters, AbDab = CoV_AbDab, method = NA, maxDist = NA, species = "Human")
head(human_perfect_match)
# example with perfect match in "Mouse" species
mouse_perfect_match <- NAb.query(bcr_clusters = COVID_01_clusters, AbDab = CoV_AbDab, method = NA, maxDist = NA, species = "Mouse")
# example with fuzzy matching with "hamming" method and max distance of 1
human_hamming_1_match <- NAb.query(bcr_clusters = COVID_01_clusters, AbDab = CoV_AbDab, method = "hamming", maxDist = 1, species = "Human")


## 3. Inter group analysis
## 3.1 V/J gene usage
## 3.1.1 Boxplot: V/J gene
# Boxplot showing the V/J gene usage of the clustered sequences between two groups.
# The parameter "colname" can be "v_call" for V gene or "j_call" for J gene.
# Statistical comparisons are carried out by the two-sided Wilcoxon rank-sum test.
gene.fre.plot(group1_seqs_list = COVID_seqs_list,
              group1_all_clustered_seqs = COVID_all_clustered_seqs,
              group1_label = "COVID-19",
              group2_seqs_list = HC_seqs_list,
              group2_all_clustered_seqs = HC_all_clustered_seqs,
              group2_label = "HC",
              colname = "v_call")
gene.fre.plot(group1_seqs_list = COVID_seqs_list,
              group1_all_clustered_seqs = COVID_all_clustered_seqs,
              group1_label = "COVID-19",
              group2_seqs_list = HC_seqs_list,
              group2_all_clustered_seqs = HC_all_clustered_seqs,
              group2_label = "HC",
              colname = "j_call")
## 3.1.2 Heatmap: V-J gene pair
# Heatmap showing the fold change of V-J gene pair frequency of clustered sequences between two groups.
# Log fold change (log FC) is calculated as the log2 ratio of the average values between group1 and group2 samples.
# Statistical comparisons are carried out by the two-sided Wilcoxon rank-sum test.
# FDR correction was performed with the Benjamini–Hochberg procedure.
# The V-J gene pair frequency of samples with a frequency less than the minimum value (1‰) is set to the minimum value.
vjpair.group.plot(group1_seqs_list = COVID_seqs_list,
                  group1_all_clustered_seqs = COVID_all_clustered_seqs,
                  group1_label = "COVID-19",
                  group2_seqs_list = HC_seqs_list,
                  group2_all_clustered_seqs = HC_all_clustered_seqs,
                  group2_label = "HC")

## 3.2 Junction length
# Density ridges showing the length distribution of junction amino acid of clustered sequences between two groups.
# Statistical comparisons are carried out by the two-sided Wilcoxon rank-sum test.
len.group.plot(group1_all_clustered_seqs = COVID_all_clustered_seqs, group1_label = "COVID-19",
               group2_all_clustered_seqs = HC_all_clustered_seqs, group2_label = "HC")

## 3.3 Diversity analysis
## 3.3.1 Number and size of clusters
# Bubble plot showing the size and number of clusters between two groups.
clu.size.plot(clusters_list1 = COVID_clusters_list, group1_label = "COVID-19",
              clusters_list2 = HC_clusters_list, group2_label = "HC")
## 3.3.2 Tcf20 score
# Tcf20 score represents the proportion of sequences attributed to the top 20 clonal families out of the total number of BCR sequences.
# Boxplot showing the Tcf20 scores between two groups.
# Statistical comparisons are carried out by the two-sided Wilcoxon rank-sum test.
Tcf20.plot(clusters_list1 = COVID_clusters_list, group1_label = "COVID-19",
           clusters_list2 = HC_clusters_list, group2_label = "HC")

## 3.4 Somatic hypermutation (SHM) analysis
# Calculate the average SHM ratio of all clusters in each sample.
# The calculation of SHM ratios may take a while.
SHM_df = SHM.calculation(clusters_list1 = COVID_clusters_list, group1_label = "COVID-19",
                         clusters_list2 = HC_clusters_list, group2_label = "HC")
# Boxplot showing the SHM ratios between two groups.
# Statistical comparisons are carried out by the two-sided Wilcoxon rank-sum test.
SHM.plot(SHM_df = SHM_df)
# Calculate the SHM ratios of clustered sequences in different isotypes in each sample.
# The calculation of SHM ratios may take a while.
SHM_iso_df = SHM.iso.calculation(clusters_list1 = COVID_clusters_list, group1_label = "COVID-19",
                                 clusters_list2 = HC_clusters_list, group2_label = "HC")
# Boxplot showing the SHM ratios of clustered sequences in different isotypes between two groups.
# Statistical comparisons are carried out by the two-sided Wilcoxon rank-sum test.
SHM.iso.plot(SHM_iso_df = SHM_iso_df)

### 3.5 NAb ratio calculation
## 3.5.1 Load the public antibody database to get the information of neutralizing antibody (NAb) sequences.
# Here is an example of the Coronavirus Antibody Database (CoV-AbDab).
CoV_AbDab = read.csv("example/CoV-AbDab_130623.csv")
# Obtain the IGHV gene, IGHJ gene and CDRH3 of all antibody sequences
v = sapply(strsplit(CoV_AbDab$Heavy.V.Gene, " "), function(x) x[1])
j = sapply(strsplit(CoV_AbDab$Heavy.J.Gene, " "), function(x) x[1])
cdr3 = sapply(strsplit(CoV_AbDab$CDRH3, " "), function(x) x[1])
vjcdr3 = unique(paste(v, j, cdr3))
## 3.5.2 NAb ratio calculation
# NAb ratio is established as an indicator of the proportional prevalence of neutralizing antibody sequences within expanded clonotypes in each sample.
# It is defined as the fraction of the number of NAb sequences within clonal families to the total number of NAb sequences present in each sample.
NAb_ratio_df = NAb.ratio.calculation(pro_data_list1 = COVID_pro_data_list, clusters_list1 = COVID_clusters_list, group1_label = "COVID-19",
                                     pro_data_list2 = HC_pro_data_list, clusters_list2 = HC_clusters_list, group2_label = "HC", NAb_vjcdr3 = vjcdr3)
# Boxplot showing the NAb ratios between the two groups.
# Statistical comparisons are carried out by the two-sided Wilcoxon rank-sum test.
NAb.ratio.plot(NAb_ratio_df)
