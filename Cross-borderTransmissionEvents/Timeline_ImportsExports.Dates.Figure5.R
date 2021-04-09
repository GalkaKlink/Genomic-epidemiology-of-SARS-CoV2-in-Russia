library(ggplot2)
#library(rstudioapi)#Safely Access the RStudio API
library(Hmisc)#Contains many functions useful for data analysis
library(grDevices) #The R Graphics Devices and Support for Colours and Fonts
library(grid) # rewrite of the graphics layout capabilities, plus some support for interaction
library(data.table) #provides an enhanced version of data.frame; useful for huge datasets
library(ggpubr) #facilitates the creation of beautiful ggplot2-based graphs for researcher with non-advanced programming backgrounds
library(tidyverse)   #includes many packages, including ggplot2
library(egg)

lct <- Sys.getlocale("LC_TIME")
Sys.setlocale("LC_TIME", "C")

#input files are results of nonrussian_lineages_from_russian_dates.strict.with_seq_number.with_countries.pl ("exports") and 
#russian_lineages_from_nonrussian_dates.strict.with_seq_number.with_countries.pl ("imports"), merged for 5 trees
exports_events = read.table(file="EXPORTS_JAN.PaperDefinition.txt")
imports_events = read.table(file="IMPORTS_JAN.PaperDefinition.txt")


dates_file1 = read.table(file="B_1_1_nodes_dates_old_labels.debugged.csv",sep=",",header=TRUE)
colnames(dates_file1) = c("int_id","temp","dates_file_date")
dates_file1 = dates_file1[,c(1,3)]
exp_b11 = subset(exports_events,exports_events$V9=="B.1.1")
imp_b11 = subset(imports_events,imports_events$V9=="B.1.1")
exports_events1 = merge(x=exp_b11,y=dates_file1,by.x="V1",by.y="int_id",all.x=TRUE,all.y=FALSE)
imports_events1 = merge(x=imp_b11,y=dates_file1,by.x="V1",by.y="int_id",all.x=TRUE,all.y=FALSE)
exports_events1$group = "exports"
imports_events1$group = "imports"


dates_file2 = read.table(file="B_1_x_nodes_dates_old_labels.debugged.csv",sep=",",header=TRUE)
colnames(dates_file2) = c("int_id","temp","dates_file_date")
dates_file2 = dates_file2[,c(1,3)]
exp_b1x = subset(exports_events,exports_events$V9=="B.1.x")
imp_b1x = subset(imports_events,imports_events$V9=="B.1.x")
exports_events2 = merge(x=exp_b1x,y=dates_file2,by.x="V1",by.y="int_id",all.x=TRUE,all.y=FALSE)
imports_events2 = merge(x=imp_b1x,y=dates_file2,by.x="V1",by.y="int_id",all.x=TRUE,all.y=FALSE)
exports_events2$group = "exports"
imports_events2$group = "imports"


dates_file3 = read.table(file="B_1_GH_nodes_dates_old_labels.debugged.csv",sep=",",header=TRUE)
colnames(dates_file3) = c("int_id","temp","dates_file_date")
dates_file3 = dates_file3[,c(1,3)]
exp_b1gh = subset(exports_events,exports_events$V9=="B.1.GH")
imp_b1gh = subset(imports_events,imports_events$V9=="B.1.GH")
exports_events3 = merge(x=exp_b1gh,y=dates_file3,by.x="V1",by.y="int_id",all.x=TRUE,all.y=FALSE)
imports_events3 = merge(x=imp_b1gh,y=dates_file3,by.x="V1",by.y="int_id",all.x=TRUE,all.y=FALSE)
exports_events3$group = "exports"
imports_events3$group = "imports"


dates_file4 = read.table(file="B_x_nodes_dates_old_labels.debugged.csv",sep=",",header=TRUE)
colnames(dates_file4) = c("int_id","temp","dates_file_date")
dates_file4 = dates_file4[,c(1,3)]
exp_bx = subset(exports_events,exports_events$V9=="B.x")
imp_bx = subset(imports_events,imports_events$V9=="B.x")
exports_events4 = merge(x=exp_bx,y=dates_file4,by.x="V1",by.y="int_id",all.x=TRUE,all.y=FALSE)
imports_events4 = merge(x=imp_bx,y=dates_file4,by.x="V1",by.y="int_id",all.x=TRUE,all.y=FALSE)
exports_events4$group = "exports"
imports_events4$group = "imports"

dates_file5 = read.table(file="A_nodes_dates_old_labels.debugged.csv",sep=",",header=TRUE)
colnames(dates_file5) = c("int_id","temp","dates_file_date")
dates_file5 = dates_file5[,c(1,3)]
exp_a = subset(exports_events,exports_events$V9=="A")
imp_a = subset(imports_events,imports_events$V9=="A")
exports_events5= merge(x=exp_a,y=dates_file5,by.x="V1",by.y="int_id",all.x=TRUE,all.y=FALSE)
imports_events5 = merge(x=imp_a,y=dates_file5,by.x="V1",by.y="int_id",all.x=TRUE,all.y=FALSE)
exports_events5$group = "exports"
imports_events5$group = "imports"


all = rbind(exports_events1,imports_events1,exports_events2,imports_events2,exports_events3,imports_events3,exports_events4,imports_events4,exports_events5,imports_events5)
all$group = ifelse(all$group=="imports","inbound transmissions","outbound transmissions")
all$dates_file_date = as.Date(all$dates_file_date,format="%Y-%m-%d")
all = all[complete.cases(all), ] 
all_g = all[,c(1,2,6,12)]
all_d = all[,c(1,2,11,12)]
colnames(all_g)= c("node","size","date","group")
all_g$date=as.Date(all_g$date)
all_g=all_g[order(all_g$date),]
order = all_g$node
colnames(all_d)= c("node","size","date","group")
all1=rbind(all_g,all_d)
all1$date=as.Date(all1$date)
all1$size=ifelse(all1$size>100,100,all1$size)
all1$group=as.factor(all1$group)
common = theme_classic()+
  theme(axis.text.x = element_text(size =16,angle= 90), 
        axis.text.y = element_blank(),
        #axis.text.y = element_text(size =16),
        axis.title.x = element_text(size =16), 
        axis.title.y = element_blank(),
        #axis.title.y = element_text(size=16),
        legend.text = element_text(size=16),
        legend.title = element_text(size=16))

common1 = theme_classic()+
  theme(axis.text.x = element_text(size =16,angle= 90), 
        axis.text.y = element_text(size =16),
        axis.title.x = element_text(size =16), 
        axis.title.y = element_text(size=16),
        legend.text = element_blank(),
        legend.title = element_blank())        


ggp = ggplot(data=all1,aes(x=date,y=factor(node,level = order)))+
geom_rect(aes(xmin = as.Date("2020-04-01"), xmax = as.Date("2020-08-01"), ymin = -Inf, ymax = Inf),
                   fill = "bisque", alpha = 0.2)+
geom_line(data=all1,aes(x=date,y=factor(node,level = order),size=size,color=group),alpha=1)+
scale_x_date(name="date", 
               limits=c(as.Date(min(all1$date), "%Y-%m-%d"),
                        as.Date(max(all1$date), "%Y-%m-%d")),
               date_breaks = "1 month",
               date_labels = "%Y-%m-%d") +
scale_size(labels=c("1","10",">100"),
  
             breaks = c(1, 10, 100),name="number of samples \n after an event")+ 
scale_color_manual(values = c("inbound transmissions" = "#d7191c","outbound transmissions" = "#fdae61"),name="")+
common
          
             
            ggp = ggp + theme(axis.title.x = element_blank(),axis.text.x = element_blank(),legend.title=element_text(size=16))
            
             
avia=read.table(file="Avia_and_Sequences.txt",header=TRUE)
avia$month=as.Date(avia$month)
av = ggplot(avia,aes(x=month,y=International)) + geom_ribbon(aes(ymin=0,ymax=International),fill="grey",color="black") +  common1 + scale_x_date(name="Earliest sample date", 
               limits=c(as.Date(min(all1$date), "%Y-%m-%d"),
                        as.Date(max(all1$date), "%Y-%m-%d")),
               date_breaks = "1 month",
               date_labels = "%Y-%m-%d")+
ylab("passenger traffic")+
theme(axis.title.x = element_blank(),axis.text.x = element_blank(),plot.margin = unit(c(3,3,3,3), "lines"))

seq = ggplot(avia,aes(x=month,y=seq)) + geom_line(color="black") + common1 + scale_x_date(name="date", 
               limits=c(as.Date(min(all1$date), "%Y-%m-%d"),
                        as.Date(max(all1$date), "%Y-%m-%d")),
               date_breaks = "1 month",
               date_labels = "%b-%d")+
ylab("number \n of samples")
#scale_x_datetime(labels = date_format("%b"))
grid.newpage()
figure = egg::ggarrange(ggp, av, seq, heights = c(0.9, 0.05,0.05))
ggsave("ImportsExports.Dates.WithExtraData.Final.png",plot= figure,
       width = 20, height = 50, dpi = 350, units = "cm") 

