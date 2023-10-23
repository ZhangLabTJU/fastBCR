library(fastBCR)

# Load data
load('example/COVID_01.Rdata')
load('example/COVID_02.Rdata')
load('example/COVID_03.Rdata')
load('example/COVID_04.Rdata')
load('example/COVID_05.Rdata')
load('example/HC_01.Rdata')
load('example/HC_02.Rdata')
load('example/HC_03.Rdata')
load('example/HC_04.Rdata')
load('example/HC_05.Rdata')

# Data processing
COVID_01 = data.pro(COVID_01)
COVID_02 = data.pro(COVID_02)
COVID_03 = data.pro(COVID_03)
COVID_04 = data.pro(COVID_04)
COVID_05 = data.pro(COVID_05)
HC_01 = data.pro(HC_01)
HC_02 = data.pro(HC_02)
HC_03 = data.pro(HC_03)
HC_04 = data.pro(HC_04)
HC_05 = data.pro(HC_05)

# BCR clustering
COVID_01_clusters = BCR.cluster(COVID_01)
COVID_02_clusters = BCR.cluster(COVID_02)
COVID_03_clusters = BCR.cluster(COVID_03)
COVID_04_clusters = BCR.cluster(COVID_04)
COVID_05_clusters = BCR.cluster(COVID_05)
HC_01_clusters = BCR.cluster(HC_01)
HC_02_clusters = BCR.cluster(HC_02)
HC_03_clusters = BCR.cluster(HC_03)
HC_04_clusters = BCR.cluster(HC_04)
HC_05_clusters = BCR.cluster(HC_05)

## treemap visualization (single sample)
plot.treemap(COVID_01, COVID_01_clusters)

# Clusters metrics
COVID_01_metrics = Metrics.df(COVID_01, COVID_01_clusters, "COVID")
COVID_02_metrics = Metrics.df(COVID_02, COVID_02_clusters, "COVID")
COVID_03_metrics = Metrics.df(COVID_03, COVID_03_clusters, "COVID")
COVID_04_metrics = Metrics.df(COVID_04, COVID_04_clusters, "COVID")
COVID_05_metrics = Metrics.df(COVID_05, COVID_05_clusters, "COVID")
HC_01_metrics = Metrics.df(HC_01, HC_01_clusters, "HC")
HC_02_metrics = Metrics.df(HC_02, HC_02_clusters, "HC")
HC_03_metrics = Metrics.df(HC_03, HC_03_clusters, "HC")
HC_04_metrics = Metrics.df(HC_04, HC_04_clusters, "HC")
HC_05_metrics = Metrics.df(HC_05, HC_05_clusters, "HC")

## Number&Size (between groups)
Covid_clu.size = sapply(c(COVID_01_clusters, COVID_02_clusters, COVID_03_clusters, COVID_04_clusters, COVID_05_clusters),
                        function(x) nrow(x))
Covid_clu_size = as.numeric(names(table(Covid_clu.size)))
Covid_clu_n = as.numeric(table(Covid_clu.size))
Healthy_clu.size = sapply(c(HC_01_clusters, HC_02_clusters, HC_03_clusters, HC_04_clusters, HC_05_clusters),
                          function(x) nrow(x))
Healthy_clu_size = as.numeric(names(table(Healthy_clu.size)))
Healthy_clu_n = as.numeric(table(Healthy_clu.size))
clu_size_df = data.frame(group = c(rep('COVID', length(Covid_clu_size)),
                          rep('HC', length(Healthy_clu_size))),
                         size = c(Covid_clu_size, Healthy_clu_size),
                         num = c(Covid_clu_n, Healthy_clu_n))
plot.clu.size(clu_size_df)

# Downstream analysis
## Diversity:Tcf20 score (between groups)
Tcf20_df = data.frame(group = c(COVID_01_metrics$Group, COVID_02_metrics$Group, COVID_03_metrics$Group, COVID_04_metrics$Group, COVID_05_metrics$Group,
                                HC_01_metrics$Group, HC_02_metrics$Group, HC_03_metrics$Group, HC_04_metrics$Group, HC_05_metrics$Group),
                      value = c(COVID_01_metrics$Tcf20, COVID_02_metrics$Tcf20, COVID_03_metrics$Tcf20, COVID_04_metrics$Tcf20, COVID_05_metrics$Tcf20,
                                HC_01_metrics$Tcf20, HC_02_metrics$Tcf20, HC_03_metrics$Tcf20, HC_04_metrics$Tcf20, HC_05_metrics$Tcf20))
plot.Tcf20(Tcf20_df)

## Clustered seqs
COVID_01_clustered_seqs = clu2df(COVID_01_clusters)
COVID_02_clustered_seqs = clu2df(COVID_02_clusters)
COVID_03_clustered_seqs = clu2df(COVID_03_clusters)
COVID_04_clustered_seqs = clu2df(COVID_04_clusters)
COVID_05_clustered_seqs = clu2df(COVID_05_clusters)
HC_01_clustered_seqs = clu2df(HC_01_clusters)
HC_02_clustered_seqs = clu2df(HC_02_clusters)
HC_03_clustered_seqs = clu2df(HC_03_clusters)
HC_04_clustered_seqs = clu2df(HC_04_clusters)
HC_05_clustered_seqs = clu2df(HC_05_clusters)
COVID_clustered_seqs = rbind(COVID_01_clustered_seqs, COVID_02_clustered_seqs, COVID_03_clustered_seqs, COVID_04_clustered_seqs, COVID_05_clustered_seqs)
HC_clustered_seqs = rbind(HC_01_clustered_seqs, HC_02_clustered_seqs, HC_03_clustered_seqs, HC_04_clustered_seqs, HC_05_clustered_seqs)

## Gene usage
### Pie chart:V/J gene (single sample)
plot.pie.freq(COVID_01_clustered_seqs, "v_call")
plot.pie.freq(COVID_01_clustered_seqs, "j_call")
### Heatmap:VJ-pair (single group)
plot.vjpair.sample(COVID_clustered_seqs)
plot.vjpair.sample(HC_clustered_seqs)
### Heatmap:VJ-pair (between groups)
plot.vjpair.group(COVID_clustered_seqs, 'COVID',
                  HC_clustered_seqs, 'HC')
### Boxplot (between groups)
COVID_vcall = table(COVID_clustered_seqs$v_call)
HC_vcall = table(HC_clustered_seqs$v_call)
uni_gene = union(names(COVID_vcall), names(HC_vcall))
v_fre = rbind(gene.fre.df(COVID_01_clustered_seqs, 'v_call', uni_gene, 'COVID'),
              gene.fre.df(COVID_02_clustered_seqs, 'v_call', uni_gene, 'COVID'),
              gene.fre.df(COVID_03_clustered_seqs, 'v_call', uni_gene, 'COVID'),
              gene.fre.df(COVID_04_clustered_seqs, 'v_call', uni_gene, 'COVID'),
              gene.fre.df(COVID_05_clustered_seqs, 'v_call', uni_gene, 'COVID'),
              gene.fre.df(HC_01_clustered_seqs, 'v_call', uni_gene, 'HC'),
              gene.fre.df(HC_02_clustered_seqs, 'v_call', uni_gene, 'HC'),
              gene.fre.df(HC_03_clustered_seqs, 'v_call', uni_gene, 'HC'),
              gene.fre.df(HC_04_clustered_seqs, 'v_call', uni_gene, 'HC'),
              gene.fre.df(HC_05_clustered_seqs, 'v_call', uni_gene, 'HC'))
plot.gene.fre(v_fre)

COVID_jcall = table(COVID_clustered_seqs$j_call)
HC_jcall = table(HC_clustered_seqs$j_call)
uni_gene = union(names(COVID_jcall), names(HC_jcall))
j_fre = rbind(gene.fre.df(COVID_01_clustered_seqs, 'j_call', uni_gene, 'COVID'),
              gene.fre.df(COVID_02_clustered_seqs, 'j_call', uni_gene, 'COVID'),
              gene.fre.df(COVID_03_clustered_seqs, 'j_call', uni_gene, 'COVID'),
              gene.fre.df(COVID_04_clustered_seqs, 'j_call', uni_gene, 'COVID'),
              gene.fre.df(COVID_05_clustered_seqs, 'j_call', uni_gene, 'COVID'),
              gene.fre.df(HC_01_clustered_seqs, 'j_call', uni_gene, 'HC'),
              gene.fre.df(HC_02_clustered_seqs, 'j_call', uni_gene, 'HC'),
              gene.fre.df(HC_03_clustered_seqs, 'j_call', uni_gene, 'HC'),
              gene.fre.df(HC_04_clustered_seqs, 'j_call', uni_gene, 'HC'),
              gene.fre.df(HC_05_clustered_seqs, 'j_call', uni_gene, 'HC'))
plot.gene.fre(j_fre)

## Junction length
### single sample
plot.len.sample(COVID_01_clustered_seqs)

### between groups
length_df = data.frame(group = c(rep(COVID_01_metrics$Group, nrow(COVID_01_clustered_seqs)),
                                 rep(COVID_02_metrics$Group, nrow(COVID_02_clustered_seqs)),
                                 rep(COVID_03_metrics$Group, nrow(COVID_03_clustered_seqs)),
                                 rep(COVID_04_metrics$Group, nrow(COVID_04_clustered_seqs)),
                                 rep(COVID_05_metrics$Group, nrow(COVID_05_clustered_seqs)),
                                 rep(HC_01_metrics$Group, nrow(HC_01_clustered_seqs)),
                                 rep(HC_02_metrics$Group, nrow(HC_02_clustered_seqs)),
                                 rep(HC_03_metrics$Group, nrow(HC_03_clustered_seqs)),
                                 rep(HC_04_metrics$Group, nrow(HC_04_clustered_seqs)),
                                 rep(HC_05_metrics$Group, nrow(HC_05_clustered_seqs))),
                       value = c(nchar(COVID_01_clustered_seqs$junction_aa),
                                 nchar(COVID_02_clustered_seqs$junction_aa),
                                 nchar(COVID_03_clustered_seqs$junction_aa),
                                 nchar(COVID_04_clustered_seqs$junction_aa),
                                 nchar(COVID_05_clustered_seqs$junction_aa),
                                 nchar(HC_01_clustered_seqs$junction_aa),
                                 nchar(HC_02_clustered_seqs$junction_aa),
                                 nchar(HC_03_clustered_seqs$junction_aa),
                                 nchar(HC_04_clustered_seqs$junction_aa),
                                 nchar(HC_05_clustered_seqs$junction_aa)))
plot.len.group(length_df)

## MSA
### DNA (single cluster)
plot.msa(COVID_01_clusters, 200, 'DNA')
### AA (single cluster)
plot.msa(COVID_01_clusters, 200, 'AA')

## Seqlogo:DNA/AA (single cluster)
### DNA (single cluster)
plot.seqlogo(COVID_01_clusters, 200, 'DNA')
### AA (single cluster)
plot.seqlogo(COVID_01_clusters, 200, 'AA')

## Clonal tree (single cluster)
plot.clonal.tree(COVID_01_clusters, 200)

## SHM
### cluster SHM (between groups)
SHM_df = data.frame(group = c(COVID_01_metrics$Group, COVID_02_metrics$Group, COVID_03_metrics$Group, COVID_04_metrics$Group, COVID_05_metrics$Group,
                              HC_01_metrics$Group, HC_02_metrics$Group, HC_03_metrics$Group, HC_04_metrics$Group, HC_05_metrics$Group),
                    value = c(COVID_01_metrics$SHM_ratio, COVID_02_metrics$SHM_ratio, COVID_03_metrics$SHM_ratio, COVID_04_metrics$SHM_ratio, COVID_05_metrics$SHM_ratio,
                              HC_01_metrics$SHM_ratio, HC_02_metrics$SHM_ratio, HC_03_metrics$SHM_ratio, HC_04_metrics$SHM_ratio, HC_05_metrics$SHM_ratio))
plot.SHM(SHM_df)

### isotype SHM (within groups)
SHM_iso_df = data.frame(group = c(rep(COVID_01_metrics$Group, 4),
                                  rep(COVID_02_metrics$Group, 4),
                                  rep(COVID_03_metrics$Group, 4),
                                  rep(COVID_04_metrics$Group, 4),
                                  rep(COVID_05_metrics$Group, 4)),
                        Isotypes = factor(rep(c("IGHD", "IGHM", "IGHA", "IGHG"), 5)), # Total number of samples
                        value = c(COVID_01_metrics$SHM_iso, COVID_02_metrics$SHM_iso, COVID_03_metrics$SHM_iso, COVID_04_metrics$SHM_iso, COVID_05_metrics$SHM_iso,
                                  HC_01_metrics$SHM_iso, HC_02_metrics$SHM_iso, HC_03_metrics$SHM_iso, HC_04_metrics$SHM_iso, HC_05_metrics$SHM_iso))
plot.SHM.iso(SHM_iso_df)

## CSR
### single cluster
plot.CSR.cluster(COVID_01_clusters, 50)

### single sample
plot.CSR.sample(COVID_01_clusters)
plot.CSR.sample(HC_01_clusters)

## NAb ratio
### boxplot (between groups)
CoV_AbDab = read.csv("example/CoV-AbDab_130623.csv")
CoV_AbDab_vjcdr3 = vector(length = nrow(CoV_AbDab))
v = sapply(strsplit(CoV_AbDab$Heavy.V.Gene, ' '), function(x) x[1])
j = sapply(strsplit(CoV_AbDab$Heavy.J.Gene, ' '), function(x) x[1])
cdr3 = sapply(strsplit(CoV_AbDab$CDRH3, ' '), function(x) x[1])
CoV_AbDab_vjcdr3 = paste(v, j, cdr3)
CoV_AbDab_vjcdr3 = unique(CoV_AbDab_vjcdr3)
NAb_v = unlist(strsplit(CoV_AbDab_vjcdr3, ' '))[seq(1,3*length(CoV_AbDab_vjcdr3),3)]
NAb_j = unlist(strsplit(CoV_AbDab_vjcdr3, ' '))[seq(2,3*length(CoV_AbDab_vjcdr3),3)]
NAb_cdr3 = unlist(strsplit(CoV_AbDab_vjcdr3, ' '))[seq(3,3*length(CoV_AbDab_vjcdr3),3)]
NAb_Ratio = data.frame(group = c(COVID_01_metrics$Group, COVID_02_metrics$Group, COVID_03_metrics$Group, COVID_04_metrics$Group, COVID_05_metrics$Group,
                                 HC_01_metrics$Group, HC_02_metrics$Group, HC_03_metrics$Group, HC_04_metrics$Group, HC_05_metrics$Group),
                value = c(NAb.ratio(COVID_01, COVID_01_clusters, NAb_v, NAb_j, NAb_cdr3),
                          NAb.ratio(COVID_02, COVID_02_clusters, NAb_v, NAb_j, NAb_cdr3),
                          NAb.ratio(COVID_03, COVID_03_clusters, NAb_v, NAb_j, NAb_cdr3),
                          NAb.ratio(COVID_04, COVID_04_clusters, NAb_v, NAb_j, NAb_cdr3),
                          NAb.ratio(COVID_05, COVID_05_clusters, NAb_v, NAb_j, NAb_cdr3),
                          NAb.ratio(HC_01, HC_01_clusters, NAb_v, NAb_j, NAb_cdr3),
                          NAb.ratio(HC_02, HC_02_clusters, NAb_v, NAb_j, NAb_cdr3),
                          NAb.ratio(HC_03, HC_03_clusters, NAb_v, NAb_j, NAb_cdr3),
                          NAb.ratio(HC_04, HC_04_clusters, NAb_v, NAb_j, NAb_cdr3),
                          NAb.ratio(HC_05, HC_05_clusters, NAb_v, NAb_j, NAb_cdr3)))
plot.NAb.ratio(NAb_Ratio)

### ROC (between groups)
roc = pROC::roc(NAb_Ratio[,1], NAb_Ratio[,2],
          levels = c('HC','COVID'),
          direction = "<",
          auc = TRUE,
          ci = TRUE,
          smooth = F)
plot.NAb.roc(roc)

## BCR data simulation
