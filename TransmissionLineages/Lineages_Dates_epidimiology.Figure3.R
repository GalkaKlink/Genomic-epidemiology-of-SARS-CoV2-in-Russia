library(ggplot2)
library(Hmisc)#Contains many functions useful for data analysis
library(grDevices) #The R Graphics Devices and Support for Colours and Fonts
library(grid) # rewrite of the graphics layout capabilities, plus some support for interaction
library(data.table) #provides an enhanced version of data.frame; useful for huge datasets
library(ggpubr) #facilitates the creation of beautiful ggplot2-based graphs for researcher with non-advanced programming backgrounds
library(tidyverse)   #includes many packages, including ggplot2


##Read initial files
Stem<-read.csv("Stemclusters_5clades.csv", sep = "\t", header = F)    #Group Name Singletone/Stem
Singletons<-read.csv("Singletons_5clades.csv", sep = "\t", header = F) #Singletones Name
Lineages<-read.csv("Lineages_5clades.csv", sep = "\t", header = F) #InternalNode ExternalNode - output of Python script. All other inputs here are made from LOG of Python script
Lineages_stem_dates<-read.csv("cluster_roots.dates.csv", sep = "\t", header = F)#InternalName date/NoRussiaOnStem

pang_clades = read.table(file="full_cat_pangolin_clades.out",sep = "\t", header = F)
pang_clades = pang_clades[,c(1,5)]
colnames(pang_clades) = c("Sample","pangolin")
pang_clades$Sample = sub("\\|.*", "", as.character(pang_clades$Sample))
pang_clades$Sample = sub(".*?\\/", "", as.character(pang_clades$Sample))

##Add_columns
Lineages_stem_dates$Type<-ifelse(Lineages_stem_dates$V2 == "No Russia on stem",
                                 "LineageNoStem", "Lineage")
##
Lineages$ID_Type<-Lineages$V1    #ID_type == name of lineage founder
Lineages_wd<-merge(x =Lineages, y = Lineages_stem_dates, by.x = "V1", by.y = "V1") #merge()= merge two data frames by common columns or row names
##
Singletons$ID_Type<-rep("SingletonNoStem", length(Singletons$V1))
Singletons$Type<-rep("SingletonNoStem", length(Singletons$V1))
Stem$Type<-Stem$V3
###Set columns_names
colnames(Stem)<-c("ID","Presample","ID_Type","Type")
colnames(Singletons)<-c("ID","Presample","ID_Type","Type")
colnames(Lineages_wd)<-c("ID","Presample","ID_Type","StemDate","Type")

####
#Merge
lct <- Sys.getlocale("LC_TIME")
Sys.setlocale("LC_TIME", "C")   #as.Date doesn't want to work in some locales. A locale object tries to capture all the defaults that can vary between countries.
Lin_Stem_Singl<-rbind(Stem,Singletons,Lineages_wd[,c(1,2,3,5)])
Lin_Stem_Singl_date<-separate(data = Lin_Stem_Singl, col = Presample, into = c("Sample", "Date"), sep = "\\|")  #Separate column "Presample" to 2 columns
Lin_Stem_Singl_date$Date2<-as.Date(Lin_Stem_Singl_date$Date, format = "%b-%d") 

vect = vector()
for (i in 1:nrow(Lin_Stem_Singl_date)) {
  dat = Lin_Stem_Singl_date[i,6]
  d <- as.POSIXlt(as.Date(dat))
d$year <- d$year-1
d = as.character(as.Date(d))
vect = c(vect,d)
}
Lin_Stem_Singl_date$Date2 = vect
Lin_Stem_Singl_date$Date2 = as.Date(Lin_Stem_Singl_date$Date2)

Lin_Stem_Singl_date = merge(x=Lin_Stem_Singl_date, y=pang_clades,by="Sample")

#vect1 = c("samp","no","no","empty1","empty1","2021-01-01","no")
#vect2 = c("samp","no","no","empty2","empty2","2021-01-01","no")
#vect3 = c("samp","no","no","empty3","empty2","2021-01-01","no")

for (i in 1:30) {
vect1 = c("samp","no","no",i,i,"2021-01-01","no")
Lin_Stem_Singl_date=rbind(Lin_Stem_Singl_date,vect1)
}


#########
###################################
#Common theme
common = theme_classic()+
  theme(axis.text.x = element_text(size =9), 
        axis.text.y = element_text(size =9),
        axis.title.x = element_text(size =12), 
        axis.title.y = element_text(size =12,angle=90),
        legend.text = element_text(size =12),
        legend.title = element_text(size =12))
        #legend.position="none")

##################
###Plot dots in line
# define the summary function
f <- function(x) {
  r <- quantile(x, probs = c(0.00, 0.25, 0.5, 0.75, 1))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}
#########
#Order of y labels
CovMetadataSub<-Lin_Stem_Singl_date[Lin_Stem_Singl_date$ID %like% "NODE",][,c("ID","Date2")] #select columns ID and Date2 for rows with "NODE" in ID
Dates_additional<-CovMetadataSub %>% group_by(ID) %>% 
  summarise(Date2 = min(as.Date(Date2,format = "%Y-%m-%d"), na.rm = T))  #%>% is used to chain functions together (package magrittr, included in tidyverse)
#Sys.setlocale("LC_TIME", lct) #to return in default locale

##Add this dates back to stem date df
Lineages_stem_dates_fin<-merge(x =Lineages_stem_dates, y =Dates_additional,
                               by.x ="V1", by.y ="ID")

##
Dates_additional_sorted<-Dates_additional[order(Dates_additional$Date2),]
lineages_order<-Dates_additional_sorted$ID
lineages_order<-c("1","2","3","4","5","Stem","6","7","8","9","10","Singleton","11","12","13","14","15","SingletonNoStem","21","22","23","24","25",lineages_order)
############
##
##Plot intervals and lineages
cols = c("A.*" = "#c4c3d0", "B.*" = "#708090", "B.1" = "#37A1ED", "B.1.*" = "#63C76D","B.1.1" = "#F9D23C", "B.1.1.*" = "#ADE053", "B.1.1.163" = "#0A14A6","B.1.1.184" = "#7641F2","B.1.1.294" = "#EB1E2C","B.1.1.31" = "#FD8930","B.1.5.*" = "#CC5BE3")

lines<-ggplot() +
  #stat_summary(data=Lin_Stem_Singl_date, 
  #             aes(x = Date2,
  #                 y = factor(ID_Type,
  #                            level =lineages_order),
  #                 fill = ID_Type),
  #             alpha =0.4,fill = "#bdbdbd", width = 0.2,
  #             fun.data = f, geom="boxplot")+
  geom_line (data=Lin_Stem_Singl_date, 
               aes(x = Date2,
                   y = factor(ID_Type,
                              level =lineages_order)),alpha = 0.4)+
  geom_point(data=Lineages_stem_dates_fin,
             aes(y = factor(V1,
                            level =lineages_order), 
                 x=as.Date(V2, "%Y-%m-%d")),
             color = "#bdbdbd")+
  geom_count(data=Lin_Stem_Singl_date,
             aes(x = Date2,
                 y = factor(ID_Type,
                            level =lineages_order),
                 fill = pangolin),
             alpha = 0.8, shape =21)+
  scale_fill_manual(values = cols,name="PANGOLIN lineages") +
  geom_segment(data = Lineages_stem_dates_fin,
               linetype =2,
               aes(y = factor(V1,
                                  level =lineages_order), 
                   x = as.Date(V2, "%Y-%m-%d"),
                   yend = V1,
                   xend = as.Date(Date2, "%Y-%m-%d")),
               color = "#bdbdbd")+
  scale_size(name="Number of samples",range=c(2,6),
             breaks = c(1, 2, 4, 10))+#max_size = 10)+

   
scale_x_date(name="Sampling date", 
               limits=c(as.Date("2020-02-23", "%Y-%m-%d"),
                        as.Date("2020-11-30", "%Y-%m-%d")),
               date_breaks = "7 days",
               date_labels = "%b-%d") +
  common +
  theme(axis.text.y = element_blank(), 
        axis.title.y = element_text(size =12),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle = 90)) +
   ylab("Russian transmission lineages") +
   guides(color = FALSE,
         fill= guide_legend(override.aes = list(size=6)))
lines


#save
ggsave("E:/SARS_CoV2/RESULTS_DEFINITELY_FOR_PAPER/Timeline/LineagesDates_PangolinColor.png",plot= lines,
       width = 20, height = 35, dpi = 350, units = "cm")



