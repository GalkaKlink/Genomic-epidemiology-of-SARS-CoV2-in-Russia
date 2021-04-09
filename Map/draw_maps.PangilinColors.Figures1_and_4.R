# Loading required packages
library(ggplot2)
library(sf)
library(scatterpie)
library(RColorBrewer)
library(stringi)
library(reshape2)
library(dplyr)
library(ggpubr)
library(gridExtra)

projection =3576 #EPSG code
#Map of Russian Federation was obtained from GADM:
#https://gadm.org/download_country_v3.html
rus.map.init <- readRDS("gadm36_RUS_1_sf.rds")#%>%st_transform(projection)
ukr.map <- readRDS("gadm36_UKR_1_sf.rds")#%>%st_transform(projection)
##subset Crimea
Crimea.map <- ukr.map[ukr.map$NAME_1 %in% c("Crimea","Sevastopol'"),]
Crimea.map$NL_NAME_1<-c("Р РµСЃРїСѓР±Р»РёРєР° РљСЂС‹Рј","Рі. РЎРµРІР°СЃС‚РѕРїРѕР»СЊ")
####
rus.map<-rbind(rus.map.init, Crimea.map)#%>%st_transform(projection)
rus.map$NL_NAME_1[rus.map$NL_NAME_1 == "РџРµСЂРјСЃРєР°СЏ РєСЂР°Р№"] <- "РџРµСЂРјСЃРєРёР№ РєСЂР°Р№"
rus.map$NL_NAME_1[rus.map$NL_NAME_1 == "РљР°РјС‡Р°С‚СЃРєР°СЏ РєСЂР°Р№"] <- "РљР°РјС‡Р°С‚СЃРєРёР№ РєСЂР°Р№"
rus.map$NL_NAME_1[rus.map$NL_NAME_1 == "Р РµСЃРїСѓР±Р»РёРєР° Р§РµС‡РµРЅРѕ-РРЅРіСѓС€СЃРєР°СЏ"] <- "Р§РµС‡РµРЅСЃРєР°СЏ СЂРµСЃРїСѓР±Р»РёРєР°"
rus.map$NL_NAME_1[rus.map$NL_NAME_1 == "Р РµСЃРїСѓМЃР±Р»РёРєР° РРЅРіСѓС€РµМЃС‚РёСЏ"] <- "Р РµСЃРїСѓР±Р»РёРєР° РРЅРіСѓС€РµС‚РёСЏ"
rus.map$NL_NAME_1[rus.map$NL_NAME_1 == "РЎР°РЅРєС‚-РџРµС‚РµСЂР±СѓСЂРі (РіРѕСЂСЃРѕРІРµС‚)"] <-"Рі. РЎР°РЅРєС‚-РџРµС‚РµСЂР±СѓСЂРі"
rus.map$NL_NAME_1[is.na(rus.map$NL_NAME_1)] <-"Рі. РњРѕСЃРєРІР°"

###############
#general graphical values
dimh<-18
dimw<-32
#
lwdp = 0.05
common = theme_classic()+
  theme(axis.line = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position="right",
        legend.text = element_text(size =10),
        legend.title = element_text(size =10),
        legend.key.size = unit(0.4, "cm"),
        plot.margin = unit(c(0,0,0,0), "lines"))
####
##################
###Draw number of genomes on map
GenomeData<-read.csv("FinalForMap.merged_pangolin_lineages.WithAncNodes.txt", sep = "\t", header = TRUE)

GenomeData_merge_MosSPb<-GenomeData
GenomeData_merge_MosSPb$NL_NAME[GenomeData_merge_MosSPb$NL_NAME == "РњРѕСЃРєРѕРІСЃРєР°СЏ РѕР±Р»Р°СЃС‚СЊ"]<-"Рі. РњРѕСЃРєРІР°"
GenomeData_merge_MosSPb$NL_NAME[GenomeData_merge_MosSPb$NL_NAME == "Р›РµРЅРёРЅРіСЂР°РґСЃРєР°СЏ РѕР±Р»Р°СЃС‚СЊ"]<-"Рі. РЎР°РЅРєС‚-РџРµС‚РµСЂР±СѓСЂРі"

#########
###upgrade of initial table for mapping fraction of unique and shared lineages
spl=split(GenomeData_merge_MosSPb,GenomeData_merge_MosSPb$Lineage)
lineage_regions = data.frame()
for(i in (1:length(spl))){
 sp=spl[i]
 sp=as.data.frame(sp)
 colnames(sp)=colnames(GenomeData_merge_MosSPb)
 un=unique(sp$NL_NAME)
 flag=ifelse(length(un)==1,"unique","shared")
 lineage_regions[i,1]=sp[1,11]
 lineage_regions[i,2]=flag
}
unique=subset(lineage_regions,lineage_regions$V2=="unique")
unique=unique$V1
GenomeData_merge_MosSPb$lineage_status=ifelse(GenomeData_merge_MosSPb$Lineage %in% unique,"unique","shared")

###without division by source:
seqs.by.region.aggr<-aggregate(GISAID_name~NL_NAME,FUN=length,data=GenomeData)
seqs.by.region.aggr[is.na(seqs.by.region.aggr)] <- 0
seqs.by.region.aggr_sum = seqs.by.region.aggr
seqs.by.region.aggr_sum$All = seqs.by.region.aggr_sum$GISAID_name
##################
##Map manipulations
rus.map_transf<-rus.map %>% st_transform(projection)
labels<-as.data.frame(cbind(rus.map$NL_NAME_1,st_centroid(rus.map_transf)))

labels_transformed <- labels %>%
  st_as_sf(crs = projection)
labels_transformed_with_lat_lon <- cbind(labels_transformed, st_coordinates(labels_transformed))

#####add map values to genome counts data frame
labels_seqs<-subset(labels_transformed_with_lat_lon, 
                    labels_transformed_with_lat_lon$NL_NAME_1 %in% seqs.by.region.aggr_sum$NL_NAME)
labels_seqs_df<-data.frame(t(apply(labels_seqs, 1, unlist)))
options(digits =14)
labels_seqs_df$X<-as.double(labels_seqs_df$X)
labels_seqs_df$Y<-as.double(labels_seqs_df$Y)
###Adjuct centroid coordinates for Arkhangelsk
labels_seqs_df$Y[labels_seqs_df$NL_NAME_1=="РђСЂС…Р°РЅРіРµР»СЊСЃРєР°СЏ РѕР±Р»Р°СЃС‚СЊ"]<-labels_seqs_df$Y[labels_seqs_df$NL_NAME_1=="РђСЂС…Р°РЅРіРµР»СЊСЃРєР°СЏ РѕР±Р»Р°СЃС‚СЊ"] -220000
labels_seqs_df$X[labels_seqs_df$NL_NAME_1=="РђСЂС…Р°РЅРіРµР»СЊСЃРєР°СЏ РѕР±Р»Р°СЃС‚СЊ"]<-labels_seqs_df$X[labels_seqs_df$NL_NAME_1=="РђСЂС…Р°РЅРіРµР»СЊСЃРєР°СЏ РѕР±Р»Р°СЃС‚СЊ"] -330000
##
seqs.by.region.aggr_sort<-seqs.by.region.aggr_sum[match(labels_seqs_df$NL_NAME_1, seqs.by.region.aggr_sum$NL_NAME),]
#######
labels_seqs_num<-cbind(labels_seqs_df,seqs.by.region.aggr_sort)
#Add radius
labels_seqs_num$radius<-log(labels_seqs_num$All+2)*30000
######################

sequenced_genomes<-ggplot() +
  geom_sf(data = rus.map,
          lwd = .5, fill ="#ffffff", color ='#bdbdbd') +
  coord_sf(crs=projection) +
  #guides(size =FALSE)+
  common+
  theme(legend.position = c(0.6,0.8))

###PANGOLIN

pangolin.df<-GenomeData_merge_MosSPb[,c(1,8,10)]


###FOR TWO TIME INTERVALS:
#GenomeData_merge_MosSPb$Date = as.Date(GenomeData_merge_MosSPb$Date)
#GenomeData_merge_MosSPb1 = subset(GenomeData_merge_MosSPb,as.POSIXct(GenomeData_merge_MosSPb$Date)<"2020-08-01")
#GenomeData_merge_MosSPb2 = subset(GenomeData_merge_MosSPb,as.POSIXct(GenomeData_merge_MosSPb$Date)>="2020-08-01")
#pangolin.df<-GenomeData_merge_MosSPb2[,c(1,8,10)]
#####



pangolin_aggr<-pangolin.df %>% group_by(NL_NAME,Pangolin_clade) %>% summarise_each(funs(length))
pangolin_for_graph<-dcast(pangolin_aggr, NL_NAME~Pangolin_clade)
pangolin_for_graph[is.na(pangolin_for_graph)] <- 0

pangolin_for_graph_sum<-transform(pangolin_for_graph, All=rowSums(pangolin_for_graph[,c(2:12)]))
pangolin_for_graph_sort<-pangolin_for_graph_sum[match(labels_seqs_df$NL_NAME_1, pangolin_for_graph_sum$NL_NAME),]

labels_pan_num<-cbind(labels_seqs_df,pangolin_for_graph_sort)
#Add radius
labels_pan_num = subset(labels_pan_num,labels_pan_num$NL_NAME != "NA")
labels_pan_num$radius<-log(labels_pan_num$All+2)*30000     


##########

colorsPang = c('#C4C3D0', '#708090',  '#37A1ED', '#63C76D', '#ADE053', '#F9D23C', '#FD8930', '#EB1E2C', '#CC5BE3', '#7641F2', '#0A14A6')
names(colorsPang) = c('A..','B..', 'B.1', 'B.1..', 'B.1.1..', 'B.1.1','B.1.1.31', 'B.1.1.294', 'B.1.5..', 'B.1.1.184', 'B.1.1.163')


numbers=labels_pan_num$All
numbers = numbers[!is.na(numbers)]
Fig_pangolin<-sequenced_genomes + geom_scatterpie(aes(x=X, y=Y, group = NL_NAME_1,
                                                      r=radius),
                                                  data = labels_pan_num,
                                                  cols=colnames(labels_pan_num[,c(17:27)]),
                                                  alpha=1) + 
                  
               geom_scatterpie_legend(labels_pan_num$radius,x = -2000000, y = -700000,n=3,labeller = function(x) round(exp(1)**(x/30000)-2)) +               
  
  scale_fill_manual(name = "PANGOLIN clades",
                    values=colorsPang,
                    breaks = c('A..','B..', 'B.1', 'B.1..', 'B.1.1..', 'B.1.1','B.1.1.31', 'B.1.1.294', 'B.1.5..', 'B.1.1.184', 'B.1.1.163'),
                    labels = c('A.*','B.*', 'B.1', 'B.1.*', 'B.1.1.*', 'B.1.1','B.1.1.31', 'B.1.1.294', 'B.1.5.*', 'B.1.1.184', 'B.1.1.163'))
Fig_pangolin


ggsave("Fig_pangolin_clades.png",plot= Fig_pangolin, 
        width = 30, height = 19, dpi=300, units = "cm")

#####################################################################
#######Map unique-shared lineages
###First two rows are for drawing number of LINEAGES instead of number of samples
GenomeData_merge_MosSPb$nonunique = paste(GenomeData_merge_MosSPb$NL_NAME,GenomeData_merge_MosSPb$Lineage,sep="_")
GenomeData_merge_MosSPb1 = GenomeData_merge_MosSPb[!duplicated(GenomeData_merge_MosSPb$nonunique),]

pangolin.df<-GenomeData_merge_MosSPb1[,c(1,12,10)]

pangolin_aggr<-pangolin.df %>% group_by(NL_NAME,lineage_status) %>% summarise_each(funs(length))
pangolin_for_graph<-dcast(pangolin_aggr, NL_NAME~lineage_status)
pangolin_for_graph[is.na(pangolin_for_graph)] <- 0
pangolin_for_graph_sum<-transform(pangolin_for_graph, All=rowSums(pangolin_for_graph[,c(2:3)]))
pangolin_for_graph_sort<-pangolin_for_graph_sum[match(labels_seqs_df$NL_NAME_1, pangolin_for_graph_sum$NL_NAME),]

labels_pan_num<-cbind(labels_seqs_df,pangolin_for_graph_sort)
#Add radius
labels_pan_num$radius<-log(labels_pan_num$All+2)*30000   #radius - number of lineages in the region

colorsPang<- c("#FD8930", "#37A1ED")
numbers=labels_pan_num$radius
numbers = numbers[!is.na(numbers)]
Fig_pangolin<-sequenced_genomes + geom_scatterpie(aes(x=X, y=Y, group = NL_NAME_1,
                                                      r=radius),
                                                  data = labels_pan_num,
                                                  cols=colnames(labels_pan_num[,c(17:18)]),
                                                  alpha=.8) + 
              scale_fill_manual(name = "lineages",values=colorsPang,labels=colnames(labels_pan_num[,c(17:18)])) +
              geom_scatterpie_legend(numbers,x = -2000000, y = -700000,n=3,labeller = function(x) round(exp(1)**(x/30000)-2))

ggsave("Fig_SharedUniaueLineages.SizeIsSamplesNumber.December.png",
plot= Fig_pangolin,path="Figures/", width = 30, height = 19, dpi=500, units = "cm")
