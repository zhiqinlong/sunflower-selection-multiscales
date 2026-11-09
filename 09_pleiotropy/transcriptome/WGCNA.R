options(repos = c(CRAN = "https://cloud.r-project.org"))
if (!("BiocManager" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  install.packages("BiocManager")
}
if (!("WGCNA" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install('WGCNA')
}
library(WGCNA)
enableWGCNAThreads(as.integer(Sys.getenv('SLURM_CPUS_PER_TASK','1')))
if (!("tidyverse" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  install.packages("tidyverse")
}

if (!('DESeq2' %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("DESeq2")
}
if (!('igraph' %in% installed.packages())) {
  # Install this package if it isn't installed yet
  install.packages("igraph")
}
library(igraph)
library(parallel)
ncores=as.integer(Sys.getenv('SLURM_CPUS_PER_TASK','1'))
library(DESeq2)
library(tidyverse)
###run from the analysis folder (r.sh cd's into it): fir /home/zlong/scratch/co_expression_135 or rorqual /home/zlong/links/scratch/transcriptome_annuus/co_expression_135
###expression_values = rsem_out/*.genes.results from rorqual transcriptome_annuus (transcriptome_analysis.sh):
###135 libraries; the 45 top-up runs (-b) were merged with their library at the fastq level before RSEM.
###The old 180-sample run (run1 and -b as separate samples) is in /home/zlong/scratch/co_expression
files=list.files('./expression_values',pattern = 'genes.results',full.names = F)
print(length(files))
stopifnot(length(files)==as.integer(Sys.getenv('N_SAMPLES','135')))
if (Sys.getenv('N_SAMPLES','135')=='135') stopifnot(!any(grepl('-b',files)))
###which edge sets to keep: Rscript WGCNA.R 5 (top5 only) or Rscript WGCNA.R 1 5 (default)
tops=as.numeric(commandArgs(trailingOnly = TRUE))
if (length(tops)==0) tops=c(1,5)
pdf(paste0('WGCNA_plots_top',paste(tops,collapse='_'),'.pdf'))
###expression data should be quantified and normalized gene expression data.
combine_tissues<- function (files) {
  data_ex=data.frame()
  for (file in files) {
    sample=gsub('.genes.results','',file)
    data=read.csv(paste('./expression_values/',file,sep=''),header=T,sep='\t') %>% 
      select(gene_id,expected_count)%>%
      dplyr::rename(!!sample :=expected_count)
    if (nrow(data_ex)==0){data_ex=data}else{data_ex=merge(data_ex,data,by='gene_id',all=T)}
  }
  return(data_ex)
}
data_ex=combine_tissues(files) 
sum(is.na(data_ex))
##filter out the genes have not been expressed or low expressed 
###change the dataformat to dataframe with column=gene and row=sample 
data_ex_filt=data_ex%>%
  mutate(gene_id=gsub('mRNA:','',gene_id))%>%
  column_to_rownames(.,var = "gene_id") %>%
  dplyr::filter(apply(.,1, function(row) any(row>10)))%>%
  mutate(across(everything(),~as.integer(round(.x))))
#data_ex_filt=data_ex_filt[sample(1500),]
######## DESeq2 package is used to normalize the counts data. It will normalize the Library size and compositinal 
###We can not use TPM/to build the network.
#####
colData <- data.frame(row.names = colnames(data_ex_filt),
                      group = factor(rep(1, ncol(data_ex_filt))))
dds <- DESeqDataSetFromMatrix(
  countData = data_ex_filt, # Our prepped data frame with counts
  colData = colData, # Data frame with annotation for our samples
  design = ~1 # Here we are not specifying a model
)
# Normalize and transform the data in the `DESeqDataSet` object using the `vst()`
# function from the `DESEq2` R package
dds_norm <- vst(dds)
# Retrieve the normalized data from the `DESeqDataSet`
normalized_counts <- assay(dds_norm) %>%
  t() # Transpose this data
#Identifying Outlier Genes and remove them 
#################################
##########1000 genes to run test
#################################
#set.seed(123);n_test_genes <- 1000

#test_genes <- sample( colnames(normalized_counts),size = min(n_test_genes, ncol(normalized_counts)),replace = FALSE)

#normalized_counts_test <- normalized_counts[, test_genes, drop = FALSE]; dim(normalized_counts_test)
#gsg <-goodSamplesGenes(normalized_counts_test)
####################################
gsg <-goodSamplesGenes(normalized_counts)
summary(gsg)
gsg$allOK ##this is Ture indicating there is no outlier genes
if (!gsg$allOK){
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(colnames(normalized_counts)[!gsg$goodGenes], collapse = ", "))); #Identifies and prints outlier genes
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(normalized_counts)[!gsg$goodSamples], collapse = ", "))); #Identifies and prints oulier samples
  normalized_counts <- normalized_counts[gsg$goodSamples == TRUE, gsg$goodGenes == TRUE] # Removes the offending genes and samples from the data
}

#Identifying outlier samples
sampleTree <- hclust(dist(normalized_counts), method = "average") #Clustering samples based on distance 

#Plotting the cluster dendrogram
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
###if you have outliers you should remove them
#cut.sampleTree <- cutreeStatic(sampleTree, cutHeight = 15, minSize = 10) #returns numeric vector
#Remove outlier
#normalized_counts <- normalized_counts[cut.sampleTree==1, ]
###Pairwise Gene Co-expression similarity

spt <- pickSoftThreshold(normalized_counts,networkType = "unsigned") 
par(mar=c(1,1,1,1))
plot(spt$fitIndices[,1],spt$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(spt$fitIndices[,1],spt$fitIndices[,2],col="red",labels=spt$fitIndices[,1])
abline(h=0.80,col="red")

plot(spt$fitIndices[,1], spt$fitIndices[,6],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(spt$fitIndices[,1], spt$fitIndices[,6], labels= spt$fitIndices[,1],col="red")
###the beta should be maximize the R2 (>0.8) and minimize the loss of mean connectivity
#### for my data, the recommanded beta=spt$powerEstimate:14. becuase the R2>0.9 and the mean connectivity is similar to when power=12
softPower <- spt$powerEstimate
###the 180-sample run used beta=14; with 135 samples check the plots in WGCNA_plots.pdf
print(spt$fitIndices)
print(paste('softPower:',softPower))
stopifnot(!is.na(softPower))

#calculating the similarity measurement and transforming the similarity by the adjacency function and generates a weighted network adjacency matrix.
adjacency <- adjacency(normalized_counts, power = softPower,type='unsigned')
### writing notes about the weight and mode
###This adjacency should show 1 represent high connection and 0 represent no corelation
calculate_interaction<-function(ad_matrix){
graph=graph_from_adjacency_matrix(ad_matrix,mode = "undirected", weighted = TRUE, diag =FALSE)
###same values as igraph::betweenness(graph,weights = NA,normalized = T) and harmonic_centrality(graph,weights = NA,normalized = TRUE),
###only split into chunks of genes run on ncores CPUs (single-core igraph took ~20 h for top5).
###betweenness = sum over source genes, so the chunks (sources) are summed and normalized as igraph::betweenness does;
###harmonic centrality is computed per gene. Checked against the single-core functions in check_parallel() below.
chunks=split(seq_len(vcount(graph)),cut(seq_len(vcount(graph)),max(ncores,1),labels=FALSE))
###stop if any chunk failed (mclapply returns an error object instead of stopping)
check_chunks<-function(res,len){stopifnot(all(sapply(res,is.numeric)),all(lengths(res)==len)); res}
betweenness<-Reduce('+',check_chunks(mclapply(chunks,function(s) igraph:::betweenness_subset_impl(graph,vids=V(graph),sources=s,targets=V(graph),weights=NA),mc.cores=ncores),vcount(graph)))
vc=as.numeric(vcount(graph)); betweenness<-2*betweenness/(vc*vc-3*vc+2)
########
###Ordinary closeness can be misleading in a disconnected network because it only considers reachable nodes. A gene in a very small component may appear highly central simply because it is close to the few genes it can reach, even though it is disconnected from most of the network.

#Use harmonic centrality instead because unreachable genes contribute zero, appropriately penalizing genes in small or isolated components.
#######
harmonic <- unlist(check_chunks(mclapply(chunks,function(s) igraph::harmonic_centrality(graph,vids=s,weights = NA,normalized = TRUE),mc.cores=ncores),lengths(chunks)))
node_degree=degree(graph,loops = FALSE)
node_strength=strength(graph,loops = FALSE)
data_inter=data.frame(gene_id=colnames(ad_matrix),
                      betweenness=betweenness,
                      harmonic_centrality=harmonic,
                      node_degree=node_degree,
                      node_strength=node_strength)
return(data_inter)
}
###check the parallel version gives the igraph single-core values on a 3000-gene part of the network before the full run
check_parallel<-function(ad_matrix){
graph=graph_from_adjacency_matrix(ad_matrix,mode = "undirected", weighted = TRUE, diag =FALSE)
par=calculate_interaction(ad_matrix)
stopifnot(isTRUE(all.equal(par$betweenness,unname(igraph::betweenness(graph,weights = NA,normalized = T)),tolerance=1e-12)),
          identical(unname(par$harmonic_centrality),unname(igraph::harmonic_centrality(graph,weights = NA,normalized = TRUE))))
print('parallel betweenness/harmonic centrality = igraph single-core values')
}

## too many nodes and it took over than 12 hours
# Exclude diagonal and duplicate matrix entries
edge_values <- adjacency[upper.tri(adjacency)]
dev.off()
# Retain approximately the strongest 1% and 5% of edges (top1 first so it is saved even if top5 runs out of time)
for (top in tops){
adj_cutoff <- quantile(edge_values, 1-top/100, na.rm = TRUE)

adjacency_filtered <- adjacency
adjacency_filtered[adjacency_filtered < adj_cutoff] <- 0
diag(adjacency_filtered) <- 0
print(paste0('top',top))
print('this is the adj_cutoff')
print(adj_cutoff)
print('this is the absolute correlation')
print(adj_cutoff^(1 / softPower))  # corresponding absolute correlation
set.seed(1); g3000=sample(ncol(adjacency_filtered),3000)
check_parallel(adjacency_filtered[g3000,g3000])
adjacency_filtered%>%calculate_interaction(.)%>%
  write.table(.,file=paste0('adjacency_interactino_top',top,'.txt'),quote=FALSE,col.names = F,row.names = F)
rm(adjacency_filtered); gc()
}
