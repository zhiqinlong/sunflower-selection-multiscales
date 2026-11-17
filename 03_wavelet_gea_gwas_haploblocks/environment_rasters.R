# Environmental background layers for the allele-frequency maps in Figs. S5C / S6C (run once, locally).
#  - Climate variables: 1961-1990 climate normals (Normal_1961_1990_<variable>.tif, ClimateNA), cropped to the USA
#    and to each taxon's map extent.
#  - Soil variables (and elevation): values measured at the sampled populations (environmental_factors.csv),
#    interpolated by inverse-distance weighting (gstat; 7 nearest populations, power 0.5) using the populations of
#    all four taxa together.
# The list of variables and their type is data/evs_types.csv. Output: all_dataframe_envs_data.Rdata ('rasters').
library(sf)
library(raster)
library(tidyr)
library(readxl)
library(ggplot2)
library(dplyr)
library(stringr)
library(scatterpie)
library(terra)
library(gstat)

USA=st_read('/Users/zhiqinlong/Desktop/computer/work/PhD/selective_sweep/analysis/data/shp/States_shapefile_shp/States_shapefile.shp')
USA= USA %>%
  filter(State_Code != 'AK'&State_Code != 'HI')
########
####crop the tiff to smaller scale
#######
####only need to run one time
if (FALSE) {
  tif_files=list.files('/Users/zhiqinlong/Desktop/computer/work/PhD/selective_sweep/analysis/wavelet/separate_scale/Normal_1961_1990_bioclim/',pattern = 'tif$',full.names = FALSE)
  for (tif_file in tif_files){
    tif_path=paste0('/Users/zhiqinlong/Desktop/computer/work/PhD/selective_sweep/analysis/wavelet/separate_scale/Normal_1961_1990_bioclim/',tif_file)
    tif=rast(tif_path)
    tif <- project(tif, crs(USA))
    r_crop <- mask(tif,vect(USA))
    tif_out=paste0('/Users/zhiqinlong/Desktop/computer/work/PhD/selective_sweep/analysis/wavelet/separate_scale/Normal_1961_1990_bioclim_USA/',tif_file)
    writeRaster(r_crop, tif_out, overwrite = TRUE)
  }
}



pop=read_excel('/Users/zhiqinlong/Desktop/computer/work/PhD/selective_sweep/analysis/data//Marco_nature/allKMap.xlsx')[,c(2,4,5,6,7)]%>%as.data.frame()%>%
  rename(Population_ID=`Population ID`)

######empty raster for interpolation
tif_file_interpolation=paste0('/Users/zhiqinlong/Desktop/computer/work/PhD/selective_sweep/analysis/wavelet/separate_scale/Normal_1961_1990_bioclim_USA/Normal_1961_1990_MAT.tif')
tif_interpolation=rast(tif_file_interpolation)
tif_interpolation <- project(tif_interpolation, crs(USA))

#######
####create the functions for color the map
####
creategroup <- function(tiff){
  colnames(tiff) <- c("x","y","bio_value")
  min_val <- floor(min(tiff$bio_value, na.rm = TRUE))
  max_val <- max(tiff$bio_value, na.rm = TRUE)
  breaks <- seq(min_val, max_val, length.out = 9)
  tiff$levels <- cut(
    tiff$bio_value,
    breaks = breaks,
    include.lowest = TRUE,
    right = FALSE
  )
  tiff$label_order <- as.numeric(tiff$levels) - 1
  return(tiff)
}
####
####plot function
sps=c('pet','fal','arg','ann')
rasters=list()
all_evs=read.csv('/Users/zhiqinlong/Desktop/computer/work/PhD/selective_sweep/analysis/wavelet/separate_scale/evs_types.csv',header=F,na.strings = "" )   # = data/evs_types.csv
for (sp in sps){
  #####set up the range
  if (sp=='ann'){
    yrange=c(20,50)
    xrange=c(-130,-90)
  }else if (sp=='arg'){
    yrange=c(26.5,30)
    xrange=c(-98.5,-96.5)
  }else if (sp=='pet'){
    yrange=c(35,42)
    xrange=c(-110,-95)
  }else if (sp=='fal'){
    yrange=c(30,40) 
    xrange=c(-115,-100)
  }
  bbox <- st_bbox(c(xmin = xrange[1],
                    xmax =  xrange[2],
                    ymin = yrange[1],
                    ymax =  yrange[2]
  ), crs = crs(USA))  
  shapefile=vect(st_as_sfc(bbox))
  interpolation_crop <- mask(tif_interpolation,shapefile)
  all_rasters=list()
  #############
  for (i in 1:nrow(all_evs)){
    ev=all_evs[i,]$V1
    if (all_evs[i,]$V2=='climate'){
      tif_file=paste0('/Users/zhiqinlong/Desktop/computer/work/PhD/selective_sweep/analysis/wavelet/separate_scale/Normal_1961_1990_bioclim_USA/Normal_1961_1990_',ev,'.tif')
      tif=rast(tif_file)
      tif <- project(tif, crs(USA))
      r_crop<- mask(tif,shapefile)
      all_rasters[[ev]]<-as.data.frame(r_crop,xy=TRUE)%>%
        drop_na()%>%
        creategroup(.)
    }else{
       evs=read.csv('/Users/zhiqinlong/Desktop/computer/work/PhD/selective_sweep/analysis/data/Marco_nature/environmental_factors.csv',na.strings = "",check.names = FALSE)%>%
        #filter(grepl(taxon_name,Taxon))%>%
        dplyr::select(latitude,longitude,ev)%>%
        rename(x=longitude,y=latitude,value=ev)%>%
         filter(value!='NA')%>%
         mutate(value=as.numeric(value))
      model <- gstat(
        formula = value ~ 1,
        locations = ~x+y,
        data=evs,
        nmax = 7,             # number of nearest populations
        set = list(idp = .5)   # inverse-distance power
      )
      
      all_rasters[[ev]] <-predict(model,newdata =as.data.frame(interpolation_crop,xy=TRUE))%>%
        dplyr::select(c(1,2,3))%>%
        creategroup(.)
    }
    rasters[[sp]]=all_rasters
  }
}
save(rasters,file='/Users/zhiqinlong/Desktop/computer/work/PhD/selective_sweep/analysis/wavelet/separate_scale/all_dataframe_envs_data.Rdata')    

