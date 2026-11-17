# PLINK allele counts per individual (--recode A) -> allele frequency per population (rows sorted by Population.ID).
# Counts are recoded so that the allele that is major across the taxon's individuals is counted.
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option(c("-s","--sample_file"), type="character", default=NULL),
  make_option(c("-o","--output"),type="character",default=NULL),
  make_option(c("-g","--genotype_path"),type="character",default=NULL))
options(error=traceback)
parser <- OptionParser(usage = "%prog -s sp.sample -o out_frq -g genotype_path",option_list=option_list)
opt = parse_args(parser)
pops=read.csv('./sunflower_4_sps_pop_info.csv')

tmp2samp <- read.csv(opt$sample_file,sep=' ',header=F)%>%
  dplyr::select(V1)%>%
  rename("ind_ID"=V1)

out_data<-fread(opt$genotype_path)%>%as_tibble()%>%dplyr::select(-c(1:6))               
nsamp <- nrow(tmp2samp)   
#convert so major allele is always 1  
out_data[,colSums(out_data) < (nsamp * 0.5 * 2)] <- 2 - out_data[,colSums(out_data) < (nsamp * 0.5 * 2)]####
out_data<- out_data %>%
        cbind(.,tmp2samp)%>%
        merge(.,pops[c('ind_ID','Population.ID')],by="ind_ID")%>%
        select(-ind_ID) %>%
        group_by(Population.ID)%>%
        summarize(across(everything(),~sum(.x)/(length(.x) * 2)))

write.table(out_data,file=opt$output,sep=' ',quote = FALSE,row.names = FALSE,col.names = FALSE)
