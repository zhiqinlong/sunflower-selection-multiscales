# map_plot(): pie charts of allele frequencies per population on a map of one environmental variable
# (used by wavelet_gea_gwas.R for Figs. S5C / S6C). 'scale' is passed but not drawn.
library(tidyr)
library(readxl)
library(ggplot2)
library(dplyr)
library(stringr)
library(scatterpie)
colors=c("#4575b7","#5588bd","#a7d6e6","#f2f8dc","#FFE3A7","#EEAF9F","#E87C5F","#D82D17")
pop=read_excel('/Users/zhiqinlong/Desktop/computer/work/PhD/selective_sweep/analysis/data//Marco_nature/allKMap.xlsx')[,c(2,4,5,6,7)]%>%as.data.frame()%>%
  rename(Population_ID=`Population ID`)
####
####plot function
map_plot<-function(sp,plot_freq,ev,scale,resdf,gene_id){
  plot_data <- merge(pop, plot_freq, by = "Population_ID") %>%
      mutate(across(c(3,4,6,7), as.numeric))
  if (sp=='ann'){
    yrange=c(20,50)
    xrange=c(-125,-90)
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
  
  map_color=data.frame(label=0:7,color=colors)%>%
    filter(label %in% unique(resdf$label_order))
 # xend=round(scale/110.67,2)
  p=ggplot()+
    geom_tile(aes(x=x,y=y,fill=factor(levels)),data=resdf,stat="identity")+
    coord_sf(xlim=xrange,ylim=yrange)+
    geom_scatterpie(aes(x=Longitude, y=Latitude),data=plot_data,cols=c('major','minor'),alpha=0.8,pie_scale = 1)+
    scale_fill_manual(values=c(map_color$color,'#007ACC','#BA0062'))+
    #annotate("segment", x = xrange[2]-1, xend = xend +xrange[2]-1, y = yrange[2]-1.5, yend = yrange[2]-1.5)+
    #annotate('text',label=paste("scale", scale,'km',sep=' '),x=xrange[2]-5,y=yrange[2]-1)+
    theme_bw()+
    labs(x="Longitude",y="Latitude")+
   labs(fill = ev) +
  theme(
  legend.title = element_text(size = 7),
  legend.text = element_text(size = 6),
  legend.key.size = grid::unit(0.3, "cm"),
  legend.spacing.y = grid::unit(0.04, "cm"),
  panel.grid = element_blank()
)
return(p)
}

