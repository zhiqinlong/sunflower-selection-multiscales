# Figs. S3 and S4: enrichment of wavelet candidate genes among (A) GEA candidate genes and (B) the haploblocks of
# Todesco et al. (2020), pooled across scales (Fig. S3) and per scale (Fig. S4). Run locally after gea_genes.R.
#  - GEA: for each wavelet gene set, the same number of genes is sampled without replacement from all annotated genes
#    on the 17 chromosomes, 1,000 times; overlap with the GEA genes (union over retained variables) gives the null.
#  - Haploblocks: the haploblocks are shuffled 1,000 times (bedtools shuffle -noOverlapping, same number and lengths,
#    seeded); the number of wavelet genes overlapping them gives the null.
#  - Fold enrichment = observed / mean of the null; grey bar = 2.5-97.5% of null / mean; P = (#null >= observed + 1)/1001.
#  - A scale without candidate genes is shown as "No candidate genes".
# Output: FigS3_enrichment_all_scales.pdf/.png, FigS4_enrichment_by_scale.pdf/.png, FigS3_S4_enrichment_results.csv
library(data.table); library(tidyverse); library(patchwork)

base_dir   = '/Users/zhiqinlong/Desktop/computer/work/PhD/selective_sweep/analysis/wavelet'
out_dir    = file.path(base_dir,'separate_scale')
bedtools   = '/Users/zhiqinlong/Desktop/computer/work/Bio_tools/bedtools2/bin/bedtools'
genome     = file.path(base_dir,'gea_loci/genome.txt')
sps        = c('pet','fal','arg','ann')
sp_levels  = c('ANN','ARG','PET','FAL')
nperm      = 1000
diamond_col= '#C0392B'   # red diamonds, as in the figure captions

## genomic background: all annotated genes on assembled chromosomes
all_genes = fread('/Users/zhiqinlong/Desktop/computer/work/PhD/annotation/FromJose2023Aug11_FINAL_INTERPRO.AED_genes.bed') %>%
  mutate(gene_id=gsub('mRNA:|ncRNA:','',V4)) %>% as_tibble() %>%
  filter(!grepl('Ha412HOChr00c',V1)) %>% distinct(gene_id) %>% pull(gene_id)

## the six analysed scales per taxon (same definition as the wavelet analysis)
pops = read.csv(file.path(base_dir,'sunflower_4_sps_pop_info.csv'))
scales_list = list()
for(sp in sps){
  tmp2samp = read.csv(file.path(base_dir,paste0(sp,'.sample')),sep='\t',header=F) %>% select(V1) %>% rename(ind_ID=V1)
  pop_coordinates = merge(pops,tmp2samp,by='ind_ID') %>% select(Population.ID,Latitude,Longitude) %>% distinct()
  DistMat = as.matrix(dist(pop_coordinates[,c('Latitude','Longitude')]))*111.3195
  minscal = quantile(DistMat[upper.tri(DistMat)],.01,na.rm=T); if(minscal==0) minscal=minscal+15
  maxscal = quantile(DistMat[upper.tri(DistMat)],.99,na.rm=T)/2
  all_scales = exp(seq(log(minscal),log(maxscal),length.out=15))
  ids = switch(sp,ann=c(4,6,7,8,9,10),pet=c(8,9,10,11,12,13),
               fal=c(6,7,8,9,10,12),arg=c(6,8,9,10,11,12))
  scales_list[[sp]] = all_scales[ids]
}

read_wavelet = function(sp){
  w = read.table(file.path(out_dir,'top0001',paste0(sp,'_top0001_adap_gene.bed')),sep='\t')
  ## exact scale matching (file stores scales rounded to 2 decimals)
  w$scale_id = sapply(w$V8, function(v){
    hit = which(abs(round(scales_list[[sp]],2) - v) < 0.006)
    if(length(hit)!=1) stop('Scale ',v,' in ',sp,' does not match an analysed scale')
    hit })
  w
}

fmt_label = function(real,p,n) ifelse(n==0,'No candidate genes',
  paste0(real,'; ',ifelse(p<.001,'P<0.001',paste0('P=',signif(p,2)))))

summarise_perm = function(sp,scale,n,real,pers){
  if(n==0) return(tibble(sp,scale,n_genes=0,real=0,fold=NA_real_,low=NA_real_,high=NA_real_,p=NA_real_))
  e = mean(pers)
  tibble(sp,scale,n_genes=n,real,fold=real/e,
         low=unname(quantile(pers/e,.025)),high=unname(quantile(pers/e,.975)),
         p=(sum(pers>=real)+1)/(nperm+1))
}

########## GEA enrichment ##########
gea_res = list()
for(k in seq_along(sps)){
  sp = sps[k]; set.seed(123+k)
  gea = read.csv(file.path(base_dir,'gea_loci',paste0(sp,'_gea_gene.bed')),sep='\t',header=F)
  gea_genes = intersect(unique(na.omit(gea$V8)),all_genes)
  w = read_wavelet(sp)
  sets = c(list(All=unique(w$V4)),
           lapply(seq_along(scales_list[[sp]]),function(i) unique(w$V4[w$scale_id==i])))
  names(sets) = c('All',seq_along(scales_list[[sp]]))
  for(s in names(sets)){
    wg = intersect(sets[[s]],all_genes); n = length(wg)
    real = length(intersect(wg,gea_genes))
    pers = if(n>0) replicate(nperm,sum(sample(all_genes,n) %in% gea_genes)) else NULL
    scale = if(s=='All') NA_real_ else scales_list[[sp]][as.integer(s)]
    gea_res[[paste(sp,s)]] = summarise_perm(sp,scale,n,real,pers)
  }
}
gea_df = bind_rows(gea_res)

########## Haploblock (inversion) enrichment ##########
tmp = tempfile('inv_'); dir.create(tmp)
count_by_scale = function(hits_file,n_scales){
  if(file.info(hits_file)$size==0) return(c(All=0,setNames(rep(0,n_scales),1:n_scales)))
  h = read.table(hits_file,sep='\t')   # V4 gene, V5 scale_id
  c(All=n_distinct(h$V4), sapply(1:n_scales,function(i) n_distinct(h$V4[h$V5==i])) %>% setNames(1:n_scales))
}
inv_res = list()
for(k in seq_along(sps)){
  sp = sps[k]; n_sc = length(scales_list[[sp]])
  inv = file.path(base_dir,'gea_loci',paste0(sp,'.inversion.bed'))
  w = read_wavelet(sp)
  wbed = file.path(tmp,paste0(sp,'_wavelet.bed'))
  w %>% transmute(V1,V2,V3,V4,scale_id) %>% distinct() %>% arrange(V1,V2) %>%
    write.table(wbed,sep='\t',quote=F,row.names=F,col.names=F)
  hits = file.path(tmp,'hits.bed'); shuf = file.path(tmp,'shuffle.bed')
  system(paste(bedtools,'intersect -u -a',wbed,'-b',inv,'>',hits))
  real = count_by_scale(hits,n_sc)
  pers = sapply(1:nperm,function(i){
    system(paste(bedtools,'shuffle -noOverlapping -seed',10000*k+i,'-i',inv,'-g',genome,'>',shuf))
    system(paste(bedtools,'intersect -u -a',wbed,'-b',shuf,'>',hits))
    count_by_scale(hits,n_sc)
  })
  n_genes = c(All=n_distinct(w$V4),setNames(sapply(1:n_sc,function(i) n_distinct(w$V4[w$scale_id==i])),1:n_sc))
  for(s in names(real)){
    scale = if(s=='All') NA_real_ else scales_list[[sp]][as.integer(s)]
    inv_res[[paste(sp,s)]] = summarise_perm(sp,scale,unname(n_genes[s]),unname(real[s]),pers[s,])
  }
  message(sp,' inversion permutations done')
}
inv_df = bind_rows(inv_res)

res = bind_rows(mutate(gea_df,test='GEA'),mutate(inv_df,test='Inversions')) %>%
  mutate(sp=toupper(sp),scale_km=round(scale,1),label=fmt_label(real,p,n_genes))
write.csv(res,file.path(out_dir,'FigS3_S4_enrichment_results.csv'),row.names=F)

########## plotting ##########
theme_enrich = theme_classic(base_size=11)+
  theme(plot.title=element_text(hjust=.5),strip.background=element_blank(),
        strip.text=element_text(size=11),plot.tag=element_text(face='bold'))

enrich_layers = function(){ list(
  geom_vline(xintercept=1,linetype='dashed',color='grey50'),
  geom_segment(aes(x=1,xend=fold,yend=y),color='grey75',na.rm=T),
  geom_segment(aes(x=low,xend=high,yend=y),color='grey70',linewidth=3,lineend='round',na.rm=T),
  geom_point(aes(x=1),shape=21,fill='white',color='grey40',size=2),
  geom_point(aes(x=fold),shape=18,color=diamond_col,size=3.2,na.rm=T),
  geom_text(aes(x=ifelse(is.na(fold),1,pmax(fold,high)),label=label),hjust=-.12,size=3),
  scale_x_continuous(expand=expansion(mult=c(.05,.55))),
  coord_cartesian(clip='off'))
}

## Fig. S3: pooled across scales
s3_panel = function(tst){
  d = res %>% filter(test==tst,is.na(scale)) %>% mutate(y=factor(sp,levels=rev(sp_levels)))
  ggplot(d,aes(y=y))+enrich_layers()+labs(title=tst,x='Fold enrichment relative to null expectation',y=NULL)+theme_enrich
}
figS3 = (s3_panel('GEA')|s3_panel('Inversions'))+plot_annotation(tag_levels='A')
ggsave(file.path(out_dir,'FigS3_enrichment_all_scales.png'),figS3,width=9,height=3,dpi=300)
ggsave(file.path(out_dir,'FigS3_enrichment_all_scales.pdf'),figS3,width=9,height=3)

## Fig. S4: per scale
s4_panel = function(tst){
  d = res %>% filter(test==tst,!is.na(scale)) %>%
    mutate(sp=factor(sp,levels=sp_levels)) %>%
    mutate(y=factor(scale_km,levels=sort(unique(scale_km),decreasing=T)))
  ggplot(d,aes(y=y))+enrich_layers()+
    facet_wrap(~sp,nrow=1,scales='free')+
    labs(title=tst,x='Fold enrichment relative to null expectation',y='Scale (km)')+theme_enrich
}
figS4 = (s4_panel('GEA')/s4_panel('Inversions'))+plot_annotation(tag_levels='A')
ggsave(file.path(out_dir,'FigS4_enrichment_by_scale.png'),figS4,width=13,height=6,dpi=300)
ggsave(file.path(out_dir,'FigS4_enrichment_by_scale.pdf'),figS4,width=13,height=6)


