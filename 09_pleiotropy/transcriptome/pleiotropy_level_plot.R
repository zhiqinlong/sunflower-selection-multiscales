###135-library version: inputs (tissue_specificity.csv, adjacency_interactino_top5_right_id.txt) come from tissue_specificity.R in new_results_135;
###outputs are written to new_results_135. The old 180-sample version is in ../old_results
library(tidyverse)
library(ggpubr)
library(data.table)
genes_pos=read.csv('/Users/zhiqinlong/Desktop/computer/work/PhD/annotation/FromJose2023Aug11_FINAL_INTERPRO.AED_genes.bed',sep='\t',header=F)%>%
  mutate(gene_id=gsub('mRNA:|ncRNA:','',V4))%>%filter(!grepl('Ha412HOChr00c',gene_id))%>%as_tibble()%>%
  select(c(gene_id,V1,V2,V3))
colnames(genes_pos)=c('gene_id','chr','s','e')
########################################################################
co_ex=read.csv("/Users/zhiqinlong/Desktop/computer/work/PhD/transcriptome/new_results_135/adjacency_interactino_top5_right_id.txt",
               col.names=c("gene_id","betweenness","harmonic_centrality","node_degree",
                           "node_strength","chr","s","e"))%>%
  select(-gene_id)%>%
  filter(node_degree>0)%>%
  left_join(genes_pos,by=c("chr","s","e"))

########################################################################
tau=read.csv('/Users/zhiqinlong/Desktop/computer/work/PhD/transcriptome/new_results_135/tissue_specificity.csv')%>%
  filter(!is.na(tau))%>%
  select(tau,V1,V4,V5)%>%
  rename(chr=V1,s=V4,e=V5)%>%
  left_join(genes_pos,by=c('chr','s','e'))
########################################################################
og_ann<-read.table("/Users/zhiqinlong/Desktop/computer/work/PhD/selective_sweep/analysis/duplication_level/Orthogroups.txt",sep=':',header=FALSE)%>%
  separate_rows(V2,sep=' ')%>%
  mutate(V2=trimws(V2))%>%
  filter(grepl('Ha412HOChr',V2))%>%
  distinct()%>%
  mutate(gene_id=V2,Orthogroup=V1)%>%
  dplyr::select(c(gene_id,Orthogroup))

dup <- fread("/Users/zhiqinlong/Desktop/computer/work/PhD/selective_sweep/analysis/duplication_level/Duplications.tsv") %>%
  filter(Support==1)
####Species-tree nodes and gene-tree nodes in orthogroups
####A species-tree node represents the common ancestor of two or more species and therefore corresponds to a speciation event.
####A gene-tree node represents the common ancestor of gene copies. It can correspond to either:
####  a speciation event, producing orthologs; or
####a gene-duplication event, producing paralogs.
####OrthoFinder reconciles each orthogroup’s gene tree with the species tree to identify gene duplication, loss, and orthology relationships.
####The species tree is shared by all orthogroups, whereas each orthogroup has its own gene tree.
####An orthogroup contains all genes descended from a single ancestral gene present at a particular species-tree node.
wrong_id_gene_pos=read.csv('/Users/zhiqinlong/Desktop/computer/work/PhD/transcriptome/old_results/2023Aug11_FINAL_INTERPRO.AED_.sorted.interpro.new.gtf',sep='\t',header=FALSE)%>%
  filter(V3=="transcript")%>%
  separate(V9, into = c("gene_id", "ts"), sep = "; transcript_")%>%
  mutate(gene_id=gsub('gene_id mRNA:',"",gene_id))%>%
  select(c(gene_id,V1,V4,V5))
ann_dup <- read.csv('/Users/zhiqinlong/Desktop/computer/work/PhD/selective_sweep/analysis/duplication_level/Duplications_per_Orthogroup.tsv',sep='\t')%>%
  merge(.,og_ann,by='Orthogroup',all=T)%>%
  merge(.,wrong_id_gene_pos,by='gene_id',all=T)%>%
  rename(n_dup=Duplications..50..support.)%>%
  select(-gene_id)%>%
  rename(chr=V1,s=V4,e=V5)%>%
  left_join(genes_pos,by=c('chr','s','e'))%>%
  filter(!is.na(n_dup)&!is.na(chr))
############
############ Combine all measurements
all_pleiotropy <- co_ex %>%
  distinct(gene_id, .keep_all=TRUE) %>%
  select(gene_id, betweenness, harmonic_centrality,
         node_degree, node_strength) %>%
  full_join(tau %>%
              distinct(gene_id, .keep_all=TRUE) %>%
              select(gene_id, tau),
            by="gene_id") %>%
  full_join(ann_dup %>%
              distinct(gene_id, .keep_all=TRUE) %>%
              select(gene_id, n_dup),
            by="gene_id") 

all_pleiotropy_log=all_pleiotropy%>%
  mutate(across(c(n_dup,betweenness,harmonic_centrality,node_degree,node_strength),~log10(.x+1)))
pairs=list(c("tau","betweenness"),c("tau","harmonic_centrality"),c("tau","node_degree"),c("tau","node_strength"),c("n_dup","betweenness"),c("n_dup","harmonic_centrality"),c("n_dup","node_degree"),c("n_dup","node_strength"),c("n_dup","tau"))
library(patchwork)
lab=list(
  tau=expression(tau),
  n_dup=expression(log[10](n_dup+1)),
  betweenness=expression(log[10](betweenness+1)),
  harmonic_centrality=expression(log[10](harmonic~centrality+1)),
  node_degree=expression(log[10](node~"degree"+1)),  # quoted: bare degree is plotmath for the ° symbol
  node_strength=expression(log[10](node~strength+1))
)

p=map(pairs,\(v)
      ggplot(all_pleiotropy_log,
             aes(x=.data[[v[1]]],y=.data[[v[2]]]))+
        geom_point(size=.5,alpha=.15,color="grey")+
        geom_smooth(method="lm",se=T,
                    color="#E94E74",linewidth=.6)+
        stat_cor(method="spearman",
                 aes(label=paste("rho~`=`~",after_stat(r),
                                 "~`,`~",after_stat(p.label))),
                 label.x.npc=.05,label.y.npc=.95)+
        labs(x=lab[[v[1]]],y=lab[[v[2]]])+
        theme_classic()
)
p_cor_all=wrap_plots(p,ncol=4)
ggsave(p_cor_all,file="/Users/zhiqinlong/Desktop/computer/work/PhD/transcriptome/new_results_135/cors_among_pleiotropic_parameters_top5.png",width=10,height=8)

ggsave(p_cor_all,file="/Users/zhiqinlong/Desktop/computer/work/PhD/transcriptome/new_results_135/cors_among_pleiotropic_parameters_top5.pdf",width=10,height=8)
########################################################################
########wavelet
########################################################################
wavelet_dir="/Users/zhiqinlong/Desktop/computer/work/PhD/selective_sweep/analysis/wavelet/separate_scale/top0001/"
sps=c("ann","pet","fal","arg")
wavelet_repeat <- map_dfr(sps, \(sp)
                          read.delim(paste0(wavelet_dir, sp, "_top0001_adap_gene.bed"),
                                     header=FALSE) %>%
                            transmute(gene_id=V4, sp=toupper(sp)) %>% distinct()
) %>%
  summarise(repeat_n_w=n_distinct(sp),
            species_w=paste(sort(sp), collapse=","),
            .by=gene_id)


########################################################################
########global adaptation
########################################################################
global_dir="/Users/zhiqinlong/Desktop/computer/work/PhD/selective_sweep/analysis/old_methods_sweep/window_based_99/"
sps=c("ann","pet","fal","arg")

global_repeat <- map_dfr(sps, \(sp)
                         read.delim(paste0(global_dir, sp, "_global_0.99_geneid.txt"),
                                    header=FALSE) %>%
                           transmute(gene_id=V1, sp=toupper(sp)) %>% distinct()
) %>%
  summarise(repeat_n_g=n_distinct(sp),
            species_g=paste(sort(sp), collapse=","),
            .by=gene_id)

df <- all_pleiotropy %>%
  left_join(wavelet_repeat,by="gene_id") %>%
  left_join(global_repeat,by="gene_id")

pars <- c("tau","n_dup","betweenness","harmonic_centrality",
          "node_degree","node_strength")


set.seed(123)

result <- crossing(
  sp=c("PET","FAL","ARG","ANN"),
  method=c("Wavelet","Global"),
  parameter=pars
) %>% pmap_dfr(\(sp,method,parameter){
  
  col <- ifelse(method=="Wavelet","species_w","species_g")
  sel <- str_detect(coalesce(df[[col]],""),
                    paste0("(^|,)",sp,"($|,)"))
  x <- df[[parameter]]
  if(parameter=="n_dup") x <- log10(x+1)
  
  candidate <- x[sel & is.finite(x)]
  bg <- x[!sel & is.finite(x)]
  
  stat <- if(parameter=="n_dup") mean else median

  n <- length(candidate)
  obs <- stat(candidate)
  null <- replicate(1000,stat(sample(bg,n)))
  
  p <- (sum(abs(null-median(null)) >=
              abs(obs-median(null)))+1)/1001
  
  tibble(sp,method,parameter,n,observed=obs,null,p)
})

###colours as in the manuscript figures (Fig. 4, Fig. S17): blue = species-wide (Global), pink = local (Wavelet)
box_cols <- c(Global="#A4DAF5",Wavelet="#F2A1BC")   # Fig. 4 null boxes (light)
point_cols <- c(Global="#68C2EE",Wavelet="#EA6290") # observed diamonds; Fig. S17 boxes
text_cols <- c(Global="#3290BF",Wavelet="#C03061")  # significance labels
method_labs <- c(Global="Species-wide adaptation",Wavelet="Local adaptation")
###panel shapes (height/width of the axis frame) measured in the manuscript figures; legend below as in the figures.
###Box widths already match (Fig. 4: width .55; Fig. S17: default .75).
fig4_shape <- theme(aspect.ratio=217/354,legend.position="bottom")      # Fig. 4 A-C: 354 x 217 px
s17_shape <- theme(aspect.ratio=276/434,legend.position="bottom")       # Fig. S17 tau, duplication
s17_cent_shape <- theme(aspect.ratio=262/450,legend.position="bottom")  # Fig. S17 centrality panels (2 x 2)
cent_labs <- c(betweenness="Betweenness",harmonic_centrality="Harmonic centrality",
               node_degree="Node degree",node_strength="Node strength")

plot_rand <- function(x,y){
  obs <- distinct(x,sp,method,parameter,observed,p) %>%
    mutate(lab=case_when(
      p<.001 ~ "***",
      p<.01  ~ "**",
      p<.05  ~ "*",
      TRUE   ~  'ns'
    ))
  pd <- position_dodge(.7)
  
  ggplot(x,aes(sp,null,fill=method))+
    geom_boxplot(width=.55,outlier.shape=NA,position=pd)+
    ###diamond colour taken from each row (a fill vector outside aes() is misassigned in faceted plots)
    geom_point(data=mutate(obs,dfill=point_cols[as.character(method)]),aes(sp,observed,group=method,fill=I(dfill)),
               inherit.aes=FALSE,shape=23,color="black",
               size=2.7,position=pd)+
    geom_text(data=obs,aes(sp,Inf,label=lab,color=method,
                           group=method),
              inherit.aes=FALSE,vjust=1.1,size=3.5,
              position=pd)+
    scale_fill_manual(values=box_cols,labels=method_labs)+
    scale_color_manual(values=text_cols,guide="none")+
    scale_y_continuous(expand=expansion(mult=c(.05,.22)))+
    labs(x=NULL,y=y,fill=NULL,color=NULL)+
    theme_classic()
}

p_tau <- plot_rand(
  filter(result,parameter=="tau"),"Tau"
)+fig4_shape
ggsave(p_tau,file="/Users/zhiqinlong/Desktop/computer/work/PhD/transcriptome/new_results_135/p_tau_singlesps.pdf",width=4,height=3)

p_dup <- plot_rand(
  filter(result,parameter=="n_dup"),
  expression("Mean "*log[10]*"(duplication events + 1)")
)+fig4_shape
ggsave(p_dup,file="/Users/zhiqinlong/Desktop/computer/work/PhD/transcriptome/new_results_135/p_dup_singlesps.pdf",width=4,height=3)

###Fig. 4C: harmonic centrality alone, same panel size as Fig. 4A and 4B
p_harmonic <- plot_rand(
  filter(result,parameter=="harmonic_centrality"),"Harmonic centrality"
)+fig4_shape
ggsave(p_harmonic,file="/Users/zhiqinlong/Desktop/computer/work/PhD/transcriptome/new_results_135/p_harmonic_singlesps.pdf",width=4,height=3)

p_cent <- plot_rand(
  filter(result,parameter %in% pars[3:6]),"Centrality"
)+facet_wrap(~parameter,nrow=1,scales="free_y",labeller=as_labeller(cent_labs))+fig4_shape
ggsave(p_cent,file="/Users/zhiqinlong/Desktop/computer/work/PhD/transcriptome/new_results_135/p_cent_singlesps.pdf",width=13,height=3)

##################

pars <- c("tau","n_dup","betweenness",
          "harmonic_centrality","node_degree","node_strength")
groups <- c("Unique","2 species","3 species","4 species")

set.seed(123)

result_g <- crossing(
  method=c("Wavelet","Global"),
  group=groups,parameter=pars
) %>% pmap_dfr(\(method,group,parameter){
  
  r <- df[[ifelse(method=="Wavelet","repeat_n_w","repeat_n_g")]]
  sel <- r==match(group,groups)
  
  x <- df[[parameter]]
  if(parameter=="n_dup") x <- log10(x+1)
  
  candidate <- x[sel %in% TRUE & is.finite(x)]
  bg <- x[is.na(r) & is.finite(x)]
  stat <- if(parameter=="n_dup") mean else median
  
  obs <- stat(candidate)
  null <- replicate(1000,stat(sample(bg,length(candidate))))
  p <- (sum(abs(null-median(null)) >=
              abs(obs-median(null)))+1)/1001
  
  tibble(group,method,parameter,n=length(candidate),
         observed=obs,null,p)
}) %>%
  mutate(group=factor(group,levels=groups))

plot_group <- function(x,y){
  obs <- distinct(x,group,method,parameter,observed,p) %>%
    mutate(lab=case_when(
      p<.001~"***",p<.01~"**",p<.05~"*",
      TRUE~'ns'
    ))
  
  pd <- position_dodge(.7)
  
  ggplot(x,aes(group,null,fill=method))+
    geom_boxplot(outlier.shape=NA,position=pd)+
    ###diamond colour taken from each row (a fill vector outside aes() is misassigned in faceted plots)
    geom_point(data=mutate(obs,dfill=point_cols[as.character(method)]),aes(group,observed,group=method,fill=I(dfill)),
               inherit.aes=FALSE,shape=23,color="black",
               size=2.5,position=pd)+
    geom_text(data=obs,aes(group,Inf,label=lab,color=method,
                           group=method),
              inherit.aes=FALSE,vjust=1.1,size=3,position=pd)+
    scale_fill_manual(values=point_cols,labels=method_labs)+
    scale_color_manual(values=text_cols,guide="none")+
    scale_y_continuous(expand=expansion(mult=c(.05,.22)))+
    labs(x=NULL,y=y,fill=NULL,color=NULL)+
    theme_classic()
}

p_tau_g <- plot_group(
  filter(result_g,parameter=="tau"),"Tau"
)+s17_shape

p_dup_g <- plot_group(
  filter(result_g,parameter=="n_dup"),
  expression("Mean "*log[10]*"(duplication events + 1)")
)+s17_shape

###Fig. S17 centrality: 2 x 2 panels, each with its own y axis title, one shared legend
p_cent_g <- wrap_plots(lapply(pars[3:6],function(v)
  plot_group(filter(result_g,parameter==v),cent_labs[[v]])+s17_cent_shape),ncol=2)+
  plot_layout(guides="collect")&theme(legend.position="bottom")


ggsave(p_tau_g,file="/Users/zhiqinlong/Desktop/computer/work/PhD/transcriptome/new_results_135/p_tau_multiple_sps.pdf",width=4.5,height=3.3)

ggsave(p_dup_g,file="/Users/zhiqinlong/Desktop/computer/work/PhD/transcriptome/new_results_135/p_dup_multiple_sps.pdf",width=4.5,height=3.3)

ggsave(p_cent_g,file="/Users/zhiqinlong/Desktop/computer/work/PhD/transcriptome/new_results_135/p_cent_multiple_sps.pdf",width=9,height=6.2)

