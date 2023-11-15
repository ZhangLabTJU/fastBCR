#### Fast BCR clonal family inference
### 1. Data loading
# Absolute path to the 'COVID'/'HC' folder in the 'example' folder in fastBCR package.
COVID_folder_path <- '~/Documents/Rpackage/fastBCR/example/COVID'
HC_folder_path <- '~/Documents/Rpackage/fastBCR/example/HC'
# Load files from 'COVID_folder_path'/'HC_folder_path' into list.
# The type of files can be 'csv', 'tsv', or 'Rdata'.
COVID_raw_data_list <- data.load(folder_path = COVID_folder_path, type = 'csv')
HC_raw_data_list <- data.load(folder_path = HC_folder_path, type = 'csv')
# Or you can load one file at a time.
# COVID_file_path <- '~/Documents/Rpackage/fastBCR/example/COVID/COVID_01.csv'
# COVID_load_file <- data.load(folder_path = COVID_file_path, type = 'csv')

### 2. Data preprocessing
# Preprocessing of raw data to meet fastBCR requirements for input data.
# The input of the function needs to meet the AIRR standard format (containing at least 'sequence_id', 'v_call', 'j_call', and 'junction_aa' information).
# Only productive sequences whose junction amino acid lengths between 9 and 26 are reserved.
# Sequences with the same 'v_call', 'j_call' and 'junction_aa' are considered identical and deduplicated.
COVID_pro_data_list <- data.preprocess(pro_data_list = COVID_raw_data_list)
HC_pro_data_list <- data.preprocess(pro_data_list = HC_raw_data_list)

### 3. BCR clonal family inferring
# Fast clonal family inference from preprocessed data.
# The 'cluster_thre' parameter represents minimal clustering criteria (minimum number of sequences needed to form a cluster) and defaults to 3.
# For high efficiency, the 'cluster_thre' is increased by 1 for every 100,000 entries of input data.
# The 'overlap_thre' parameter represents overlap coefficient threshold for merging two clusters, selectable range (0,1) and defaults to 0.1.
# Lower 'overlap_thre' may lead to overclustering while higher thresholds may lead to the split of clonal families.
# The 'consensus_thre' parameter represents the consensus score threshold for filtering candidates and defaults to 0.8.
# A higher 'consensus_thre' means stricter inference of the cluster.
COVID_clusters_list <- data.BCR.clusters(pro_data_list = COVID_pro_data_list, cluster_thre = 3, overlap_thre = 0.1, consensus_thre = 0.8)
HC_clusters_list <- data.BCR.clusters(pro_data_list = HC_pro_data_list, cluster_thre = 3, overlap_thre = 0.1, consensus_thre = 0.8)

### 4. Classification of clustered and unclustered sequences
# Merge all the clustered sequences in each sample into 'clustered_seqs'.
# Merge all the unclustered sequences in each sample into 'unclustered_seqs'.
COVID_seqs_list <- Clustered.seqs(pro_data_list = COVID_pro_data_list, clusters_list = COVID_clusters_list)
HC_seqs_list <- Clustered.seqs(pro_data_list = HC_pro_data_list, clusters_list = HC_clusters_list)
# All the clustered/unclustered sequences from different samples in a group can also be merged.
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


#### Downstream analysis of clonal families
### 1. Summary of clusters from a sample
# Summarize the number of clusters, the average size of clusters and the proportion of clustered sequences.
COVID_clusters_summary = Clusters.summary(pro_data_list = COVID_pro_data_list, clusters_list = COVID_clusters_list)
HC_clusters_summary = Clusters.summary(pro_data_list = HC_pro_data_list, clusters_list = HC_clusters_list)
COVID_01_summary = COVID_clusters_summary$COVID_01
print(COVID_01_summary)
HC_01_summary = HC_clusters_summary$HC_01
print(HC_01_summary)

### 2. Diversity analysis
# (1) Number and size of clusters
# Bubble plot showing the size and number of clusters between two groups
clu.size.plot(clusters_list1 = COVID_clusters_list, group1_label = 'COVID',
              clusters_list2 = HC_clusters_list, group2_label = 'HC')

# (2) Tcf20 score
# Tcf20 score represents the proportion of sequences attributed to the top 20 clonal families out of the total number of BCR sequences.
# Boxplot showing the Tcf20 scores between two groups.
# Statistical comparisons are carried out by the two-sided Wilcoxon rank-sum test.
Tcf20.plot(clusters_list1 = COVID_clusters_list, group1_label = 'COVID',
           clusters_list2 = HC_clusters_list, group2_label = 'HC')

### 3. Visualization of clusters from a sample
# Point diagram showing clusters from a sample where a circle represents a cluster.
# The size and color of the circle represents the size of the cluster.
Clusters.visualization(pro_data_list = COVID_pro_data_list, clusters_list = COVID_clusters_list, index = 1)
Clusters.visualization(pro_data_list = HC_pro_data_list, clusters_list = HC_clusters_list, index = 1)

### 4. V/J gene usage
## (1) Pie chart: V/J gene
# Pie chart showing the frequency of gene usage.
# The top ten most frequent genes are shown, and the rest are represented by 'others'.
# The parameter 'colname' can be 'v_call' for V gene or 'j_call' for J gene.
# (single sample)
COVID_01_clustered_seqs = COVID_seqs_list$clustered_seqs$COVID_01
HC_01_clustered_seqs = HC_seqs_list$clustered_seqs$COVID_01
pie.freq.plot(clustered_seqs = COVID_01_clustered_seqs, colname = 'v_call')
pie.freq.plot(clustered_seqs = HC_01_clustered_seqs, colname = 'v_call')
# (single group)
pie.freq.plot(clustered_seqs = COVID_all_clustered_seqs, colname = 'j_call')
pie.freq.plot(clustered_seqs = HC_all_clustered_seqs, colname = 'j_call')

## (2) Boxplot: V/J gene (between groups)
# Boxplot showing the V/J gene usage of the clustered sequences between two groups.
# The parameter 'colname' can be 'v_call' for V gene or 'j_call' for J gene.
# Statistical comparisons are carried out by the two-sided Wilcoxon rank-sum test.
gene.fre.plot(group1_seqs_list = COVID_seqs_list,
              group1_all_clustered_seqs = COVID_all_clustered_seqs,
              group1_label = 'COVID',
              group2_seqs_list = HC_seqs_list,
              group2_all_clustered_seqs = HC_all_clustered_seqs,
              group2_label = 'HC',
              colname = 'v_call')
gene.fre.plot(group1_seqs_list = COVID_seqs_list,
              group1_all_clustered_seqs = COVID_all_clustered_seqs,
              group1_label = 'COVID',
              group2_seqs_list = HC_seqs_list,
              group2_all_clustered_seqs = HC_all_clustered_seqs,
              group2_label = 'HC',
              colname = 'j_call')

## (3) Heatmap: V-J gene pair
# Heatmap showing the frequency of V-J gene pair.
# (single sample)
vjpair.sample.plot(clustered_seqs = COVID_01_clustered_seqs)
# (single group)
vjpair.sample.plot(clustered_seqs = COVID_all_clustered_seqs)
vjpair.sample.plot(clustered_seqs = HC_all_clustered_seqs)
# # (between groups) # TODO
# vjpair.group.plot(group1_all_clustered_seqs = COVID_all_clustered_seqs, group1_label = 'COVID',
#                   group2_all_clustered_seqs = HC_all_clustered_seqs, group2_label = 'HC')

### 5. Junction length
# Histogram and density plot showing junction amino acid length of clustered sequences.
# (single sample)
len.sample.plot(clustered_seqs = COVID_01_clustered_seqs)
# (single group)
len.sample.plot(clustered_seqs = COVID_all_clustered_seqs)
# Density ridges showing junction amino acid length of clustered sequences between groups.
# Statistical comparisons are carried out by the two-sided Wilcoxon rank-sum test.
len.group.plot(group1_all_clustered_seqs = COVID_all_clustered_seqs, group1_label = 'COVID',
               group2_all_clustered_seqs = HC_all_clustered_seqs, group2_label = 'HC')
# Or between clustered sequences and unclustered sequences in a group.
len.group.plot(group1_all_clustered_seqs = COVID_all_clustered_seqs, group1_label = 'Clustered',
               group2_all_clustered_seqs = COVID_all_unclustered_seqs, group2_label = 'Unclustered')

### 6. Multiple sequence alignment (MSA)
# Visualization of multiple sequence alignment (MSA) of junction sequences within a cluster.
# The parameter 'index' allows you to choose a cluster for visualization.
# The parameter 'type' can be 'DNA' for deoxyribonucleic acid or 'AA' for amino acid.
COVID_01_clusters = COVID_clusters_list$COVID_01
msa.plot(bcr_clusters = COVID_01_clusters, index = 200, type = 'DNA')
msa.plot(bcr_clusters = COVID_01_clusters, index = 200, type = 'AA')

### 7. Sequence logo
# Visualization of sequence logo of junction sequences within a cluster
seqlogo.plot(bcr_clusters = COVID_01_clusters, index = 200, type = 'DNA')
seqlogo.plot(bcr_clusters = COVID_01_clusters, index = 200, type = 'AA')

### 8. Clonal tree
# Reconstructing B cell lineage trees with minimum spanning tree and genotype abundances using ClonalTree
clonal.tree.plot(bcr_clusters = COVID_01_clusters, index = 200, file_path = '~/Documents/Rpackage/fastBCR/ClonalTree')

### 9. Affinity maturation analysis
## (1) Somatic hypermutation (SHM)
# Boxplot showing the SHM ratios between two groups.
# The calculation of SHM ratios may take a while.
# Statistical comparisons are carried out by the two-sided Wilcoxon rank-sum test.
SHM.plot(clusters_list1 = COVID_clusters_list, group1_label = 'COVID',
         clusters_list2 = HC_clusters_list, group2_label = 'HC')
# Boxplot showing the SHM ratios of clustered sequences in different isotypes between two groups.
# Statistical comparisons are carried out by the two-sided Wilcoxon rank-sum test.
SHM.iso.plot(clusters_list1 = COVID_clusters_list, group1_label = 'COVID',
             clusters_list2 = HC_clusters_list, group2_label = 'HC')

## (2) Class switch recombination (CSR)
# Visualization of isotype co-occurrence within a BCR cluster.
# Circle size represents the number of sequences carrying a given isotype.
# Lines connecting two circles indicate the enrichment level of observing switches in the two corresponding immunoglobulin subclasses.
# The enrichment level is the ratio of observed and expected switches if immunoglobulin isotypes are assumed to be independently distributed among cluster.
# (single cluster)
CSR.cluster.plot(bcr_clusters = COVID_01_clusters, index = 50)
# (single sample)
CSR.sample.plot(bcr_clusters = COVID_01_clusters)
CSR.sample.plot(bcr_clusters = HC_01_clusters)

### 10. NAb query and NAb ratio calculation
## (1) Load the public antibody database to get the information of neutralizing antibody (NAb) sequences.
# Here is an example of the Coronavirus Antibody Database (CoV-AbDab).
CoV_AbDab = read.csv("example/CoV-AbDab_130623.csv")
# Obtain the IGHV gene, IGHJ gene and CDRH3 of all antibody sequences
v = sapply(strsplit(CoV_AbDab$Heavy.V.Gene, ' '), function(x) x[1])
j = sapply(strsplit(CoV_AbDab$Heavy.J.Gene, ' '), function(x) x[1])
cdr3 = sapply(strsplit(CoV_AbDab$CDRH3, ' '), function(x) x[1])
vjcdr3 = unique(paste(v, j, cdr3))

## (2) NAb query
# Query the corresponding sequence from the public antibody database.
# The parameter 'method' represents the CDRH3 matching method.
# It can be 'NA' for perfect match, 'hamming' for hamming distance or 'lv' for Levenshtein distance. Defaults to 'NA'.
# The parameter 'species' can be 'Mouse' or 'Human'. Defaults to 'Human'.
# The parameter 'maxDist' represents the maximum distance allowed for matching when the argument 'method' is 'hamming' or 'lv'. Defaults to 'NA'.
# example to perfect match
human_perfect_match <- NAb.query(bcr_clusters = COVID_01_clusters, AbDab = CoV_AbDab, method = NA, maxDist = NA, species = 'Human')
head(human_perfect_match)
# example with perfect match in 'Mouse' species
mouse_perfect_match <- NAb.query(bcr_clusters = COVID_01_clusters, AbDab = CoV_AbDab, method = NA, maxDist = NA, species = 'Mouse')
# example with fuzzy matching with 'hamming' method and max distance of 1
human_hamming_1_match <- NAb.query(bcr_clusters = COVID_01_clusters, AbDab = CoV_AbDab, method = 'hamming', maxDist = 1, species = 'Human')

## (3) NAb ratio calculation
# NAb ratio is established as an indicator of the proportional prevalence of neutralizing antibody sequences within expanded clonotypes in each sample.
# It is defined as the fraction of the number of NAb sequences within clonal families to the total number of NAb sequences present in each sample.
# Boxplot showing the NAb ratios between the two groups.
# Statistical comparisons are carried out by the two-sided Wilcoxon rank-sum test.
NAb.ratio.plot(pro_data_list1 = COVID_pro_data_list, clusters_list1 = COVID_clusters_list, group1_label = 'COVID',
               pro_data_list2 = HC_pro_data_list, clusters_list2 = HC_clusters_list, group2_label = 'HC', NAb_vjcdr3 = vjcdr3)

