# Per-scale wavelet outliers and locally adaptive genes, and Fig. S1 (run locally after 02_wavelet/wavelet_scans.sh).
# For each taxon:
#  1. Reads the GV P values of the six analysed scales (snp_var_pvalue.R output: <sp>.pvalue.maf0.1.Rdata, ...).
#     The six scales are scale indices (of 15 log-spaced scales) ann 4,6,7,8,9,10; pet 8-13; fal 6,7,8,9,10,12;
#     arg 6,8,9,10,11,12 - the same as in 02_wavelet/snp_var_pvalue.R and Fig. 1B.
#  2. Outliers: P in the lowest 0.1% over the six scales pooled -> top0001/<sp>_top0001_adap_loci.bed (with scale).
#  3. Fig. S1: Manhattan plot of P values per scale, with haploblocks (<sp>top0001.pvalue.png).
#  4. Genes whose body +/- 10 kb contains an outlier -> top0001/<sp>_top0001_adap_gene.bed
#     (chr, start and end of gene +/- 10 kb, gene ID, outlier start, outlier end, P, scale in km), used by every
#     later analysis.
library(data.table)
library(ggplot2)
library(tidyverse)
library(RColorBrewer)
sps=c('pet','fal','arg','ann')

setwd('/Users/zhiqinlong/Desktop/computer/work/PhD/selective_sweep/analysis/wavelet/separate_scale/wavelet_output')
pops=read.csv('../../sunflower_4_sps_pop_info.csv')

for (sp in sps){
  load(paste(sp,'.pvalue.maf0.1.Rdata',sep=''))
  load(paste(sp,'.wavvar.maf0.1.Rdata',sep=''))
  load(paste(sp,'.rm_snp.Rdata',sep=''))
  snp_pos=fread(paste(sp,'.maf0.1.snp.pos',sep=''),header=F)%>%dplyr::select(c(V1,V4))%>%as_tibble()
  
  tmp2samp <- read.csv(paste('../../',sp,'.sample',sep=''),sep='\t',header=F)%>%
    dplyr::select(V1)%>%
    dplyr::rename("ind_ID"='V1')
  #calculate a distance matrix between sampling locations
  pop_coordinates<- merge(pops,tmp2samp,by='ind_ID')%>%
    dplyr::select(c(Population.ID,Latitude, Longitude))%>%
    distinct()
  DistMat<- as.matrix(dist(pop_coordinates[,c('Latitude', 'Longitude')]))*111.3195
  minscal <- quantile(DistMat[upper.tri(DistMat)], 0.01,na.rm=T)
  # minimum scale = 1st percentile of pairwise distances (15 km if that is 0); maximum = half the 99th percentile
  if (minscal==0) {minscal=minscal+15}
  maxscal <- quantile(DistMat[upper.tri(DistMat)], 0.99,na.rm=T)/2 
  if (sp == 'ann') {
    scales <- exp(seq(log(minscal), log(maxscal), length.out = 15))[c(4,6,7,8,9,10)]
  } else if (sp == 'pet') {
    scales <- exp(seq(log(minscal), log(maxscal), length.out = 15))[c(8,9,10,11,12,13)]
  } else if (sp == 'fal') {
    scales <- exp(seq(log(minscal), log(maxscal), length.out = 15))[c(6,7,8,9,10,12)]
  } else if (sp == 'arg') {
    scales <- exp(seq(log(minscal), log(maxscal), length.out = 15))[c(6,8,9,10,11,12)]
  }
  
  pvalue_var=data.frame(var=c(scalWavVar[1, ],scalWavVar[2, ],scalWavVar[3, ],scalWavVar[4, ],scalWavVar[5, ],scalWavVar[6, ]),
                        pvalue=c(jl_out[[1]][[3]],jl_out[[2]][[3]],jl_out[[3]][[3]],jl_out[[4]][[3]],jl_out[[5]][[3]],jl_out[[6]][[3]]),
                        qvalue=c(jl_out[[1]][[4]],jl_out[[2]][[4]],jl_out[[3]][[4]],jl_out[[4]][[4]],jl_out[[5]][[4]],jl_out[[6]][[4]]),
                        scale=rep(round(scales[1:6],2),each=length(jl_out[[5]][[3]])))
  if (length(rm_snp)!=0){
    snp_pos=snp_pos[-rm_snp,]
  }
  pos_new=do.call("rbind", replicate(6,snp_pos , simplify = FALSE))

  plot_data=cbind(pos_new,pvalue_var)%>%dplyr::rename(chr=V1,pos=V4)
  
  top0001=quantile(plot_data$pvalue,0.001)
  dir.create("../top0001", showWarnings = FALSE)

  local_adap=plot_data[plot_data$pvalue<=top0001,]%>%mutate(end=pos)%>%mutate(pos=(pos-1))
  write.table(local_adap[!duplicated(local_adap[,c('chr','pos','end','pvalue','scale')]),c('chr','pos','end','pvalue','scale')],file=paste('../top0001/',sp,'_top0001_adap_loci.bed',sep=''),sep='\t',col.names=F,quote=F,row.names=F)
  chr_len <- plot_data%>%
    group_by(chr) %>% 
    summarise(chrlen=max(pos))
  chr_pos <- chr_len  %>%
    mutate(total = cumsum(as.numeric(chrlen)) - chrlen) %>%
    dplyr::select(-chrlen)
  Snp_pos <- chr_pos %>%
    left_join(plot_data, ., by="chr") %>%
    arrange(chr, pos) %>%
    mutate( poscum = pos + total) %>%
    mutate(chr_num = as.numeric(gsub("Ha412HOChr", "", chr)),
      chr_group = ifelse(chr_num %% 2 == 0, "#040000", "#D5D5D5"))
  
  X_axis <-  Snp_pos %>% group_by(chr) %>% summarize(center=( max(poscum) + min(poscum) ) / 2 )
  ######inversion
  inversions=read.csv('/Users/zhiqinlong/Desktop/computer/work/PhD/selective_sweep/analysis/data/Marco_nature/all_inversions_location.csv')
  if (sp=='fal'){sp_inv='pet'}else if(sp=='ann'){sp_inv='annuus'}else {sp_inv=sp}
  
  inversion_sp=inversions%>%
    filter(sps==sp_inv)%>%
    mutate(ystart=0.7,yend=1,START=START*1000000,END=END*1000000)%>%
    mutate(chr=case_when(CHR<10~paste('Ha412HOChr0',CHR,sep=''),
                         CHR>=10~paste('Ha412HOChr',CHR,sep='')))
  
  inversion_sp_new<-chr_pos %>%
    left_join(inversion_sp, ., by="chr") %>%
    arrange(chr, START) %>%
    mutate( start_cum = START + total,
            end_cum = END + total)
  
  ############
  
  #pdf(paste(sp,'.pvalue.pdf',sep=''),width=14,height=4)
  #png(paste(sp,'.maf0.1.pvalue.png',sep=''),width=1400,height=400)
  line_df <- data.frame(threshold=-log10(top0001),type="threshold")
  p_label <- bquote(italic(P-value) == .(format(as.numeric(top0001),
                                          scientific=TRUE,digits=2)))
  p_text <- formatC(as.numeric(top0001),format="e",digits=1)
  x_label <- paste0(
    toupper(sp),
    "&nbsp;&nbsp;&nbsp;&nbsp;",
    "       ",
    "<span style='color:#3366CC;'>— —</span>&nbsp;",
    "<i>P</i> = ",p_text
  )
  p1 <- ggplot() +
    geom_rect(data=inversion_sp_new,
              aes(xmin=start_cum,xmax=end_cum,ymin=ystart,ymax=yend),
              fill="darkblue",colour="white",linewidth=.1) +
    geom_point(data=subset(Snp_pos,pvalue<.1 & pvalue>top0001),
               aes(poscum,-log10(pvalue)),colour="grey65",size=.25,alpha=.4) +
    geom_point(data=subset(Snp_pos,pvalue<=top0001),
               aes(poscum,-log10(pvalue)),colour="#D33662",size=.7,alpha=.9) +
    geom_hline(yintercept=-log10(top0001),
               colour="#3366CC",linewidth=.5,linetype="longdash") +
    facet_wrap(~scale,strip.position="right",scales="free_y",ncol=1) +
    scale_linetype_manual(name=NULL,values="longdash",label=p_label) +
    scale_x_continuous(
      breaks=X_axis$center,
      labels=gsub("Ha412HO","",X_axis$chr),
      expand=c(.01,0),
      guide=guide_axis(angle=45,check.overlap=FALSE)
    ) +
    scale_y_continuous(expand=c(0,0)) +
    labs(x=x_label,y=expression(-log[10](italic(P)))) +
    guides(linetype=guide_legend(
      override.aes=list(colour="#3366CC",linewidth=.6)
    )) +
    theme_classic(base_size=11) +
    
    theme(
      legend.position="none",
      axis.title.x=ggtext::element_markdown(size=12),
      panel.border=element_rect(colour="black",fill=NA,linewidth=.35),
      panel.spacing.y=grid::unit(.12,"cm"),
      strip.background=element_blank(),
      strip.text=element_text(angle=0,size=9),
      axis.text.x=element_text(size=8)
    )
  
  ggsave(paste0("../",sp,"top0001.pvalue.png"),
         p1,width=5,height=4,dpi=800,bg="white")
}

bedtools="/Users/zhiqinlong/Desktop/computer/work/Bio_tools/bedtools2/bin/bedtools"
library(dplyr)
library(data.table)
library(ggpubr)
genes_pos=read.csv('/Users/zhiqinlong/Desktop/computer/work/PhD/annotation/FromJose2023Aug11_FINAL_INTERPRO.AED_genes.bed',sep='\t',header=F)%>%
  mutate(gene_id=gsub('mRNA:|ncRNA:','',V4))%>%as_tibble()
for (sp in sps){
  adap_file=paste('/Users/zhiqinlong/Desktop/computer/work/PhD/selective_sweep/analysis/wavelet/separate_scale/top0001/',sp,'_top0001_adap_loci.bed',sep='')
  gene_file='/Users/zhiqinlong/Desktop/computer/work/PhD/annotation/FromJose2023Aug11_FINAL_INTERPRO.AED_genes_updown10kb.bed'
  out_gene=paste('/Users/zhiqinlong/Desktop/computer/work/PhD/selective_sweep/analysis/wavelet/separate_scale/top0001/',sp,'_top0001_adap_gene.bed',sep='')
  lines=system(paste('wc -l', adap_file,sep=' '),intern = T)
  print(lines)
  command=paste(bedtools, 'intersect -a', adap_file,'-b',gene_file,'-wa -wb ', sep=' ')
 system(command,intern = T)%>%
    read.table(text=.,header=F,sep = "\t")%>%
    mutate(gene_id=V9)%>%
    merge(.,genes_pos,by='gene_id',all=F)%>%
    dplyr::select(c(2,8,9,1,3,4,5,6))%>%
    write.table(.,out_gene,quote=F, row.names=F,col.names = F,sep='\t')
}

for (sp in sps){
 adap_file=paste('/Users/zhiqinlong/Desktop/computer/work/PhD/selective_sweep/analysis/wavelet/separate_scale/top0001/',sp,'_top0001_adap_loci.bed',sep='')
  df=read.csv(adap_file,sep='\t',header=F)%>%
    select(c(V1,V2,V3))%>%
    distinct()%>%
    nrow()
  print(df)
 }
