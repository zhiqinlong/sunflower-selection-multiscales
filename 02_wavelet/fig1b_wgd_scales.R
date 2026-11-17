# Fig. 1B: mean genome-wide wavelet genetic dissimilarity (WGD) across the 15 scales for each taxon, with the null
# (2.5-97.5% of 100 permutations of population locations, grey). Scales to which fewer than five locations
# contribute are open symbols; the six scales used for the locus-specific scans are circled in orange.
# Input: <sp>.GenomicWavDissim.Rdata from wavelet_scans.sh step 4. Output: <sp>_genomewide_genetic_dissimilarity.pdf
library(tidyverse)
library(data.table)
setwd('/Users/zhiqinlong/Desktop/computer/work/PhD/selective_sweep/analysis/wavelet/separate_scale')
pops=read.csv('../sunflower_4_sps_pop_info.csv')
sps=c('pet','fal','arg','ann')
for (sp in sps){
tmp2samp <- read.csv(paste(sp,'.sample',sep=''),sep='\t',header=F)%>%
  dplyr::select(V1)%>%
  dplyr::rename("ind_ID"='V1')
#calculate a distance matrix between sampling locations
pop_coordinates<- merge(pops,tmp2samp,by='ind_ID')%>%
    select(c(Population.ID,Latitude, Longitude))%>%
    distinct()
DistMat<- as.matrix(dist(pop_coordinates[,c('Latitude', 'Longitude')]))*111.3195
minscal <- quantile(DistMat[upper.tri(DistMat)], 0.01,na.rm=T)
if (minscal==0) {minscal=minscal+15}
maxscal <- quantile(DistMat[upper.tri(DistMat)], 0.99,na.rm=T) / 2
scales <- exp(seq(log(minscal), log(maxscal), length.out = 15)) 

load(paste('wavelet_output/',sp,'.GenomicWavDissim.Rdata',sep=''))

tmpylim <- c(min(gWavDiss_out[['null_95bounds']][,1],na.rm=T), max(apply(gWavDiss_out[['obs_wavelets']], 2, mean, na.rm = T), na.rm = T))

#plot the change in mean genomic wavelet dissimilarity across scales, compared to null
pdf(paste(sp,'_genomewide_genetic_dissimilarity.pdf',sep=''),width=4,height=4)
plot.new()
plot.window(xlim = range(scales), ylim =tmpylim, log = 'xy')
polygon(c(scales, rev(scales)), c(gWavDiss_out[['null_95bounds']][,1], rev(gWavDiss_out[['null_95bounds']][,2])), col = grey(0.6), border = NA)
text(scales[5], mean(gWavDiss_out[['null_95bounds']][,2]), 'Null', pos = 4, col = gray(0.6))
lines(scales, apply(gWavDiss_out[['obs_wavelets']], 2, mean, na.rm = T))
text(scales[5], mean(gWavDiss_out[['null_95bounds']][,1])*1.5, 'Observed', pos = 4)
#scales with fewer than 5 contributing locations (isolated locations are masked in gWavDiss) are drawn as open symbols
n_used <- colSums(!is.na(gWavDiss_out[["obs_wavelets"]]))
#scales used for the locus-specific scans (orange circles; same indices as snp_var_pvalue.R)
sel_idx <- switch(sp, ann=c(4,6,7,8,9,10), pet=c(8,9,10,11,12,13), fal=c(6,7,8,9,10,12), arg=c(6,8,9,10,11,12))
points(scales[sel_idx], apply(gWavDiss_out[["obs_wavelets"]], 2, mean, na.rm=T)[sel_idx],
       pch=21, col="#F28C28", lwd=2, cex=2.1)
points(scales, apply(gWavDiss_out[["obs_wavelets"]], 2, mean, na.rm=T),
       pch=21, bg=ifelse(n_used >= 5, "#17499D", "white"), col=ifelse(n_used >= 5, "grey", "#17499D"), cex=1.2)
axis(1, at=scales, labels=round(scales, 0), las=2, cex.axis=.8)
axis(2)
box()
title(xlab = paste0(sp,':s (spatial scale km)'))
title(ylab = expression(Mean~of~wavelet~dissimilarities~D[list(a,b)]^{wav} ~(s)), line = 2)
dev.off()
}
