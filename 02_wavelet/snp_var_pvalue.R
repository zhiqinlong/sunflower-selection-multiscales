# Scale-specific genetic variance (GV) of every SNP at the six analysed scales of a taxon, and upper-tail P values.
# GV = variance across locations of the wavelet-transformed allele frequency divided by the SNP's standard deviation.
# Locations whose nearest neighbour is farther than half the scale are masked at that scale. A chi-square
# distribution is fitted by maximum likelihood (jl_iterator, OutFLANK-style: 2.5% of values trimmed at each tail,
# outliers with q < 0.05 removed iteratively) and gives the P values.
# Six analysed scales (indices of the 15 scales): ann 4,6,7,8,9,10; pet 8-13; fal 6,7,8,9,10,12; arg 6,8,9,10,11,12.
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(optparse))
source('./WaveletSpatialGenetic/JL_OutflankFuns5_Fst.R')
source("./WaveletSpatialGenetic/WavFuns5.R")
source("./WaveletSpatialGenetic/Likelihood functions for OutFLANK_JRL.R")
option_list <- list(
  make_option(c("-s","--order_file"), type="character", default=NULL),
  make_option(c("-n","--spname"), type="character", default=NULL),
  make_option(c("-p","--output_pvalue"),type="character",default=NULL),
  make_option(c("-m","--define_minscale"),type="integer",default=NULL),
  make_option(c("-g","--genotype_path"),type="character",default=NULL),
  make_option(c("-r","--remove_snp"),type="character",default=NULL),
  make_option(c("-v","--output_wavvar"),type="character",default=NULL))
options(error=traceback)
parser <- OptionParser(usage = "%prog -s sp.sample -p pvalue -v wavvar  -g ./sp/ -m define_minscale",option_list=option_list)
opt = parse_args(parser)
pops=read.csv('./sunflower_4_sps_pop_info.csv')
pop_coordinates<-fread(opt$order_file)%>%as_tibble()%>%
	dplyr::select(c(1))%>%
	rename(Population.ID=V1)%>%
	merge(.,pops,by='Population.ID')%>%
	select(c(Population.ID,Latitude, Longitude))%>%
	distinct()
#calculate a distance matrix between sampling locations
DistMat<- as.matrix(dist(pop_coordinates[,c('Latitude', 'Longitude')]))*111.3195


DistMat_scale=DistMat ###used to define the scale. if data is one individual from one pop, there is no need do this

minscal <- quantile(DistMat_scale[upper.tri(DistMat_scale)], 0.01,na.rm=T)
if (minscal==0) {minscal=minscal+opt$define_minscale}
maxscal <- quantile(DistMat_scale[upper.tri(DistMat_scale)], 0.99,na.rm=T)/2 
 
if (opt$spname == 'ann') {
  scales <- exp(seq(log(minscal), log(maxscal), length.out = 15))[c(4,6,7,8,9,10)]
} else if (opt$spname == 'pet') {
  scales <- exp(seq(log(minscal), log(maxscal), length.out = 15))[c(8,9,10,11,12,13)]
} else if (opt$spname == 'fal') {
  scales <- exp(seq(log(minscal), log(maxscal), length.out = 15))[c(6,7,8,9,10,12)]
} else if (opt$spname == 'arg') {
  scales <- exp(seq(log(minscal), log(maxscal), length.out = 15))[c(6,8,9,10,11,12)]
}
#I have not yet tested this with missing allele calls. Currently I recommend using it with imputed data.
fr_fun<-function(path){
  out_data<-fread(path)%>%as_tibble()
  stopifnot(identical(out_data$V1,pop_coordinates$Population.ID))   # populations in the same order
	out_data<-out_data[,-1]
  return(out_data)
}
SNPsF<- do.call(cbind,
                lapply(list.files(path = opt$genotype_path,full.names = T,pattern = '.freq'), fr_fun))

#the SNP matrix is just a matrix with alternate allele counts for each individual. Rows are samples, columns are SNPs
SNPsF <-as.matrix(SNPsF)
#convert so major allele is always 1 
#calculate the genome-wide wavelet dissimliarity for the selected scales, using permutations to generate a null ('nperm' gives the number of permutations)

#distance matrix

locSD <- apply(SNPsF, 2, sd)# locus standard deviation
####when all individuals showed 0/1 genotype the sd will be 0
rm_snp=which(locSD==0)
save(rm_snp,file=opt$remove_snp)
if(length(rm_snp)>0){
SNPsF=SNPsF[,-c(rm_snp)]
locSD=locSD[-c(rm_snp)]
}
locSD <- matrix(locSD, byrow = T, nrow = nrow(SNPsF), ncol = ncol(SNPsF))
#empty matrix for scaled-specific genetic variances
scalWavVar <- matrix(NA, nrow = length(scales), ncol = ncol(SNPsF))

outFtail <- 0.025

jl_out <- list()

for(i in 1:length(scales)){ #run test for all scales
  
  tmpW <- twavf_allOmega_b(SNPsF, OmegaDist = DistMat, s = scales[i]) #SNPs should be matrix
  
  #mask those scale-sample combinations where samples are too isolated
  tmpdist <- DistMat
  tmpdist[tmpdist==0]<-NA
  minsclz <- apply(tmpdist, 1, min, na.rm = T) * 2 #to avoid calculating wavelets for small scales for isolated samples
  tmpW[minsclz > scales[i],] <- NA #for those that have no neighbors closer than scale / 2
  scalWavVar[i, ] <- apply(tmpW/locSD, 2, var, na.rm = T)# locus variance in wavelet
  # starting value for the chi-square fit (df guess); results were the same for other starting values
  landscape_size=maxscal*0.1
	jl_out[[i]] <- jl_iterator(scalWavVar[i, ], landscape_size/scales[i], quantile(scalWavVar[i, ], outFtail), quantile(scalWavVar[i, ], 1 - outFtail))
}

save(jl_out,file=opt$output_pvalue)
save(scalWavVar,file=opt$output_wavvar)

