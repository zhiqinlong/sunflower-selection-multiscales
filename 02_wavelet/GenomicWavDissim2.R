# Genome-wide wavelet genetic dissimilarity (WGD) at 15 log-spaced scales, with a null from 100 permutations of the
# population locations (gWavDiss, WaveletSpatialGenetic). Distances: Euclidean on latitude/longitude x 111.3 km.
# Minimum scale = 1st percentile of pairwise distances (-m if that is 0); maximum = half the 99th percentile.
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(foreach))
suppressPackageStartupMessages(library(doParallel)) 
suppressPackageStartupMessages(source("./WaveletSpatialGenetic/WavFuns5.R"))
suppressPackageStartupMessages(library(optparse))
option_list <- list(
  make_option(c("-s","--order_file"), type="character", default=NULL),
  make_option(c("-o","--output"),type="character",default=NULL),
  make_option(c("-m","--define_minscale"),type="integer",default=NULL),
  make_option(c("-p","--genotype_path"),type="character",default=NULL))
options(error=traceback)
parser <- OptionParser(usage = "%prog -s sp.sample -o sp.genomicwavdissim -p ./sp/ -m define_minscale",option_list=option_list)
opt = parse_args(parser)
pops=read.csv('./sunflower_4_sps_pop_info.csv')
###if one population includes many individulas we need transfer it to allele frequency. other wise the null distribution will be influenced 

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
maxscal <- quantile(DistMat_scale[upper.tri(DistMat_scale)], 0.99,na.rm=T)/2 #largest scale = about half the landscape extent
scales <- exp(seq(log(minscal), log(maxscal), length.out = 15)) 



#I have not yet tested this with missing allele calls. Currently I recommend using it with imputed data.

fr_fun<-function(path){
  out_data<-fread(path)%>%as_tibble()
  stopifnot(identical(out_data$V1,pop_coordinates$Population.ID))   # populations in the same order
	out_data<-out_data[,-1]
  return(out_data)
}
#the SNP matrix is just a matrix with alternate allele counts for each individual. Rows are samples, columns are SNPs
SNPsF<- do.call(cbind,
                lapply(list.files(path = opt$genotype_path,full.names = T,pattern = 'freq'), fr_fun))
SNPsF <-as.matrix(SNPsF)

set.seed(1) #for reproducibility
#calculate the genome-wide wavelet dissimliarity for the selected scales, using permutations to generate a null ('nperm' gives the number of permutations)
gWavDiss_out <- gWavDiss(SNPsF, OmegaDist = DistMat, scales = scales, nperm = 100) 
save(gWavDiss_out,file=opt$output)
