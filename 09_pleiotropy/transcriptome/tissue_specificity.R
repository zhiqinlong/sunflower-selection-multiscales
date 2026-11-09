library(tidyverse)
library(tibble)
bedtools="/Users/zhiqinlong/Desktop/computer/work/Bio_tools/bedtools2/bin/bedtools"
setwd('/Users/zhiqinlong/Desktop/computer/work/PhD/transcriptome/new_results_135/')
###expression_values = rsem_out/*.genes.results from rorqual transcriptome_annuus (pipeline_rorqual/transcriptome_analysis_135.sh):
###135 libraries; the 45 top-up runs (-b) were merged with their library at the fastq level before RSEM.
###The old 180-sample version is in ../old_results
files=list.files('./expression_values',pattern = 'genes.results',full.names = F)
print(length(files))
stopifnot(length(files)==135, !any(grepl('-b',files)))

combine_tissues<- function (files) {
  data_ex=data.frame()
  for (file in files) {
  sample=gsub('.genes.results','',file)
  if  (nrow(data_ex) == 0){
    data=read.csv(paste('./expression_values/',file,sep=''),header=T,sep='\t') %>%  dplyr::select(gene_id,TPM)%>%
      dplyr::rename(!!sample :=TPM)
    data_ex=data
  }
  else{
    data_ex=read.csv(paste('./expression_values/',file,sep=''),header=T,sep='\t') %>% dplyr::select(gene_id,TPM)%>%
    dplyr::rename(!!sample :=TPM)%>%
    merge(.,data_ex,by='gene_id',all=T)
  }
  }
  return(data_ex)
}
###the gtf is the reference annotation used to build the RSEM index (kept in old_results, not duplicated here)
gene_pos=read.csv('/Users/zhiqinlong/Desktop/computer/work/PhD/transcriptome/old_results/2023Aug11_FINAL_INTERPRO.AED_.sorted.interpro.new.gtf',sep='\t',header=FALSE)%>%
  filter(V3=="transcript")%>%
  separate(V9, into = c("gene_id", "ts"), sep = "; transcript_")%>%
  mutate(gene_id=gsub('gene_id mRNA:',"",gene_id))%>%
  select(c(gene_id,V1,V4,V5))

data_ex=combine_tissues(files)
data_ex_filt=data_ex%>%
  mutate(gene_id=gsub('mRNA:','',gene_id),gene_id=gsub('ncRNA:','',gene_id))%>%
  column_to_rownames(.,var = "gene_id") %>%
 # filter(apply(.,1,max)>1)
  filter(apply(.,1, function(row) any(row>10)))
###sample name = <plant>_<tissue>_<treatment>_<rep>; match the tissue field exactly
sample_tissue=unlist(lapply(colnames(data_ex_filt),function(x){strsplit(x,"_")[[1]][2]}))
tissues=unique(sample_tissue)
number_tissues=table(sample_tissue)
###135 libraries: anthere 23, leaf 39, nectary 23, ovary 36, pistil 14
print(number_tissues)
stopifnot(length(tissues)==5)

calculate_mean_ex<-function(data_ex_filt){
  data_ex_mean=data.frame(gene_id=row.names(data_ex_filt))%>%column_to_rownames(.,var = "gene_id")
  for (t in tissues) {
    data_ex_mean[,t]=apply(data_ex_filt[,sample_tissue==t],1,mean)
  }
  return(data_ex_mean)
}

N=length(tissues)
data_ex_mean=calculate_mean_ex(data_ex_filt)
data_ex_mean=1-data_ex_mean/apply(data_ex_mean,1,max)
data_ex_mean$tau=apply(data_ex_mean,1,sum)/(N-1)
###output
data_ex_mean%>% rownames_to_column("gene_id")%>%
merge(.,gene_pos,by='gene_id')%>%write.csv(.,file='./tissue_specificity.csv',quote = F,row.names = FALSE)

###adjacency_interactino_top*.txt = WGCNA.R output from fir:/home/zlong/scratch/co_expression_135 (copy here first)
read.csv('./adjacency_interactino_top5.txt',
                           col.names=c('gene_id','betweenness','harmonic_centrality','node_degree','node_strength'),header=F,sep=' ')%>%
  merge(.,gene_pos,by="gene_id")%>%
  write.csv(.,file='./adjacency_interactino_top5_right_id.txt',quote = F,row.names = FALSE)


read.csv('./adjacency_interactino_top1.txt',
         col.names=c('gene_id','betweenness','harmonic_centrality','node_degree','node_strength'),header=F,sep=' ')%>%
  merge(.,gene_pos,by="gene_id")%>%
  write.csv(.,file='./adjacency_interactino_top1_right_id.txt',quote = F,row.names = FALSE)
