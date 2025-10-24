#####This code is written by PV. It is meant to read in the output from PV_minian package, integrate this data with telemetry measurements and annotated events, and visualize the data all together.

library(tidyverse)
library(ggpubr)
library(ggprism)
library(ggridges)
library(patchwork)
library(rjson)
library(nlme)
library(chron)
library(gtools)
summarize<-dplyr::summarize
source("C:/Users/paulv/Documents/GitHub/minian_output_analysis/functions.R")

##### Things to set manually:
#Paths/directories:
exp_direc<-"C:/Users/paulv/Box/correalab/Member Folders/Paul Vander/Experiments"
output_dir<-"C:/Users/paulv/Box/correalab/Member Folders/Paul Vander/Data/Torpor project cross-experiment analyses/Miniscope"
direcs<-c("250417_circulating_E2_torpor_miniscope/pre-OVX_torpor")
dir_struct<-c("test","mouse","start_date","session","start_time","camera") #set this to the levels of metadata stored in the directory structure of the data. Set in order from the first/higher/root directory level, to the last/lowest/final branch directory level
separator<-"/" #set depending on OS

#Data to exclude:
mice_to_exclue<-c()
session_ids_to_exclude<-c()

#Global graph settings:
ms<-theme_prism()

##### Read in miniscope data:
#Initalize variables
A<-tibble()
C<-tibble()
S<-tibble()
YrA<-tibble()
motion<-tibble()
ts<-tibble()

#Read in csv files
setwd(exp_direc)
for (dir in direcs){
  # print(paste0(".",dir))
  for (mouse in list.dirs(dir, full.names=F,recursive=F)){
    # print(paste0("..",mouse))
    for (start_date in list.dirs(paste(dir,mouse,sep=separator),full.names=F,recursive=F)){
      # print(paste0("...",start_date))
      for (session in list.dirs(paste(dir,mouse,start_date,sep=separator),full.names=F,recursive=F)){
        # print(paste0("....",session))
        for (start_time in list.dirs(paste(dir,mouse,start_date,session,sep=separator),full.names=F,recursive=F)){
          if ("A.csv" %in% list.files(paste(dir,mouse,start_date,session, start_time, "minian",sep=separator))){
            path = paste(dir,mouse,start_date,session,start_time,sep=separator)
            print(paste0("Reading ", path))
            exp<-strsplit(strsplit(dir,separator)[[1]][1],"_")[[1]][1]
            bad_cells<-read_cell_label(paste(path,"minian","cell_label.csv", sep=separator))
            A<-rbind(A,read_csv_minian(paste(path,"minian","A.csv", sep=separator))%>%filter(A>0))
            C<-rbind(C, read_csv_minian(paste(path,"minian","C.csv", sep=separator)))
            S<-rbind(S, read_csv_minian(paste(path,"minian","S.csv", sep=separator)))
            YrA<-rbind(YrA, read_csv_minian(paste0(path,"minian","YrA.csv")))
            motion<-rbind(motion, read_motion(path))
            ts<-rbind(ts, read_timestamps(paste(path,"My_V4_Miniscope","timeStamps.csv", sep=separator)))
          }
        }
      }
    }
  }
}

#Create dataframe with all timeseries data
df<-merge(C,S)%>%
  merge(YrA)%>%
  merge(ts)%>%
  merge(motion)%>%
  mutate(start_chron = chron(dates.=start_date, times.=start_time, format=c(dates="y_m_d",times = "h_m_s"),out.format=c(date="m-d-y",time="h:m:s")),
                                                 chron = start_chron + (time_ms / (1000*60*60*24)),
                                                 init_unit_id = unit_id, #preserve original unit_id labels)
                                                 unit_id = factor(unit_id)%>%as.numeric()%>%factor() #renumbers unit_id starting from 1
                                                 )%>%
  scale_temporal()

#Clean up variables and memory
rm(C,S,YrA,ts)
gc()

##### Read in telemetry data
read<-c()
t_df<-tibble()
for (dir in direcs){
  telem_file<-paste(strsplit(direcs[1],separator)[[1]][1],"telem data","telem data combined.csv", sep=separator)
  if (telem_file %in% read){next}
  print(telem_file)
  t_df<-rbind(t_df, read_telemetry_data(telem_file))
  read <- c(read, telem_file)
}

##### Merge telemetry and miniscope data
#Create time bin column as basis for merging
df<-df%>%mutate(chron_bin = cut(chron, seq(min(t_df$chron)-1/(24*60*2),max(t_df$chron)+(1/(24*60*2)), (1/(24*60)))))
t_df<-t_df%>%mutate(chron_bin = cut(chron, seq(min(t_df$chron)-(1/(24*60*2)),max(t_df$chron)+(1/(24*60*2)), (1/(24*60)))))

#Preserve all miniscope data, and assign one temperature/act to all frames in the given range
df<-merge(df, t_df, by=c("mouse","chron_bin"), all=T)
df<-rbind(df%>%filter(is.na(chron.x))%>%rename(chron=chron.y)%>%select(!chron.x), #chron.x is time measured from miniscope data, chron.y is time measured from telemetry data (this line assigns "chron" to miniscope data time, when it is present, and telemetry data time when it isn't)
          df%>%filter(!is.na(chron.x))%>%rename(chron=chron.x)%>%select(!chron.y))

##### Graph
setwd(output_dir)

#Lines with motion
ls<-list(scale_x_continuous(expand=c(0,0)),
         theme(text=element_text(size=32)))

p1<-ggplot(df%>%filter(!is.na(YrA)),aes(x=frame,y=unit_id))+
  geom_hline(yintercept = 0:length(unique(df$unit_id)),color="grey")+
  geom_ridgeline(aes(height=scaled_YrA),scale=0.9,linewidth=0.35,fill=NA)+
  ms+ls
p2<-ggplot(df, aes(x=frame,y=motion_distance))+
  geom_line()+
  ms+ls
save_png_large("line plot and motion",plot=p2/p1+plot_layout(heights=c(1,10)),w=32,h=25)

p1<-ggplot(df%>%filter(!is.na(YrA),start_time=="01_55_00"),aes(x=frame,y=unit_id))+
  geom_hline(yintercept = 0:length(unique(df$unit_id)),color="grey")+
  geom_ridgeline(aes(height=scaled_YrA),scale=0.9,fill=NA)+
  ms+ls+scale_x_continuous(expand=c(0,0),breaks=seq(10000,20000,2500))
p2<-ggplot(df%>%filter(start_time=="01_55_00"), aes(x=frame,y=motion_distance))+
  geom_line()+
  ms+ls
save_png_large("line plot and motion subset",plot=p2/p1+plot_layout(heights=c(1,10)),w=18,h=25)

ggplot(df, aes(x=frame,y=temp))+
  geom_line()+
  ms
