#####This code is written by PV. It is meant to read in the output from PV_minian package, integrate this data with telemetry measurements and annotated events, and visualize the data all together.

library(tidyverse)
library(ggpubr)
library(ggprism)
library(ggridges)
library(patchwork)
library(rjson)
library(nlme)
library(gtools)
library(corrr)
summarize<-dplyr::summarize
source("C:/Users/paulv/Documents/GitHub/minian_output_analysis/functions.R")

##### Things to set manually:
#Paths/directories:
exp_direc<-"C:/Users/paulv/Box/correalab/Member Folders/Paul Vander/Experiments"
output_dir<-"C:/Users/paulv/Box/correalab/Member Folders/Paul Vander/Data/Torpor project cross-experiment analyses/Miniscope"
direcs<-c("250417_circulating_E2_torpor_miniscope/pre-OVX_torpor")
dir_struct<-c("test","mouse","start_date","session","start_time","camera") #set this to the levels of metadata stored in the directory structure of the data. Set in order from the first/higher/root directory level, to the last/lowest/final branch directory level
separator<-"/" #set depending on OS

#Global graph settings:
ms<-list(theme_prism())

##### Read in miniscope data:
### Initalize variables
A<-tibble()
C<-tibble()
S<-tibble()
YrA<-tibble()
motion<-tibble()
ts<-tibble()
bad_frames<-list()

### Read in csv files
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
            path <- paste(dir,mouse,start_date,session,start_time,sep=separator)
            msrun_dir <- grep("msRun",list.files(path),value = T)
            print(paste0("Reading ", path))
            exp<-strsplit(strsplit(dir,separator)[[1]][1],"_")[[1]][1]
            bad_cells<-read_cell_label(paste(path,"minian","cell_label.csv", sep=separator))
            A<-rbind(A,read_csv_minian(paste(path,"minian","A.csv", sep=separator))%>%filter(A>0))
            C<-rbind(C, read_csv_minian(paste(path,"minian","C.csv", sep=separator)))
            S<-rbind(S, read_csv_minian(paste(path,"minian","S.csv", sep=separator)))
            YrA<-rbind(YrA, read_csv_minian(paste(path,"minian","YrA.csv", sep=separator)))
            motion<-rbind(motion, read_motion(paste(path,msrun_dir,sep=separator)))
            ts<-rbind(ts, read_timestamps(paste(path,"My_V4_Miniscope","timeStamps.csv", sep=separator)))
            bad_frames[[path]]<-setdiff(unique(motion$frame), unique(C$frame)) ##Detects frame in motion, but not C. This should be equal to bad_frames set in minian
            motion<-motion%>%filter(frame %nin% bad_frames[[path]])
          }
        }
      }
    }
  }
}

### Create dataframe with all timeseries data
# Merge timeseries data and uniquely identify each unit/session with an ID
df<-merge(C,S)%>%
  merge(YrA)%>%
  merge(ts)%>%
  merge(motion)%>%
  mutate(start_ts = ymd_hms(paste(start_date, start_time)),
         miniscope_ts = (start_ts + milliseconds(time_ms)),
         init_unit_id = unit_id, #preserve original unit_id labels
         unit_id = factor(unit_id)%>%as.numeric()%>%factor() #renumbers unit_id starting from 1
         )%>%
  mutate(unit_id_id = paste0(mouse,"_",start_date,"_",session,"_",unit_id),
         session_id = paste0(mouse,"_",start_date,"_",session))%>%
  scale_temporal()

A<-A%>%mutate(unit_id_id = paste0(mouse,"_",start_date,"_",session,"_",unit_id),
              session_id = paste0(mouse,"_",start_date,"_",session))

#Create new column for total time within each session
df<-df%>%
  group_by(session_id)%>%
  mutate(session_start_ts = min(miniscope_ts))%>%
  ungroup()%>%
  mutate(session_time_minutes = (interval(session_start_ts, miniscope_ts)%>%as.period(unit = "seconds")%>%as.numeric())/60)

# Clean up variables and memory
rm(C,S,YrA,ts,motion)
gc()

##### Read in telemetry data
read<-c()
t_df<-tibble()
for (dir in direcs){
  telem_file<-paste(strsplit(direcs[1],separator)[[1]][1],"telem data","telem data combined.csv", sep=separator)
  if (telem_file %in% read){next}
  print(telem_file)
  t_df<-rbind(t_df, read_telemetry_data(telem_file)) ##Ignore warning about additional pieces
  read <- c(read, telem_file)
}

##### Merge telemetry and miniscope data
# Create time bin column as basis for merging
ts_seq<-seq(min(t_df$telem_ts)-seconds(30),max(t_df$telem_ts)+seconds(30), seconds(60)) #assumes a 1-minute (60 second) sampling interval
df<-df%>%mutate(ts_bin = cut(miniscope_ts, ts_seq))
t_df<-t_df%>%mutate(ts_bin = cut(telem_ts, ts_seq))

# Preserve all miniscope data, and assign one temperature/act to all frames in the given range
df<-merge(df, t_df, by=c("mouse","ts_bin"), all=T)
df<-df%>%mutate(ts = ifelse(is.na(miniscope_ts), telem_ts, miniscope_ts)) #set timestamp to miniscope timestamp (when available), otherwise use telem timestamp

# Create 1-minute bins in miniscope data
sumdf<-df%>%
  filter(!is.na(YrA))%>%
  group_by(ts_bin,unit_id_id)%>%
  arrange(ts_bin)%>%
  mutate(mean_C=mean(C), mean_S=mean(S), mean_YrA=mean(YrA), mean_motion_distance=mean(motion_distance))%>%
  ungroup()%>%
  distinct(ts_bin,unit_id_id, .keep_all = T)%>%
  scale_temporal_bin()

##### Graph
gc()
setwd(output_dir)

### Lines with motion (quality check)
ls<-list(scale_x_continuous(expand=c(0,0), breaks=NULL),
         theme(text=element_text(size=28),plot.title = element_text(size=28)),
         labs(x="Time (minutes)"),
         facet_wrap(vars(start_time), nrow=1,scales="free_x"),
         theme(strip.text.x = element_blank()))
ridge_set<-list(ridgeline_guide(),
                ridgeline(aes(height=scaled_YrA)),
                scale_fill_viridis_c(),
                labs(title="Normalized YrA (detrended + demixed GCaMP signal)", y="Cell ID",x=element_blank()))
motion_set<-list(geom_line(),
                 labs(y="pixels",x=element_blank(),title="Motion distance"))
temp_set<-list(geom_line(linewidth = 1),
               scale_x_continuous(expand=c(0,0),breaks=seq(0,5000,4)),
               labs(y="Deg. C",title="Core body tempeature"))

# All frames
p1<-ggplot(df%>%filter(!is.na(YrA)),aes(x=session_time_minutes,y=unit_id))+ms+ls+ridge_set
p2<-ggplot(df%>%filter(!is.na(motion_distance)), aes(x=session_time_minutes,y=motion_distance))+ms+ls+motion_set
p3<-ggplot(df%>%filter(!is.na(YrA)), aes(x=session_time_minutes, y=temp))+ms+ls+temp_set
save_png_large("line plot and motion and temp",plot=p2/p1/p3+plot_layout(heights=c(1,10,1)),w=32,h=25)

# Subset of frames
s1<-p1+filter(p1$data, start_time=="01_55_00"|start_time=="02_56_59")
s2<-p2+filter(p1$data, start_time=="01_55_00"|start_time=="02_56_59")
s3<-p3+filter(p1$data, start_time=="01_55_00"|start_time=="02_56_59")
s2/s1/s3+plot_layout(heights=c(1,10,1))
save_png_large("line plot and motion and temp subset",plot=s2/s1/s3+plot_layout(heights=c(1,10,1)),w=18,h=25)

###YrA - body temperature relationship
data<-sumdf%>%filter(!is.na(mean_YrA))
data<-merge(data, data%>%group_by(unit_id)%>%summarise(cor = cor(temp, scaled_mean_YrA, method = "pearson")))%>%mutate(unit_id_cor = paste0(unit_id," (",round(cor,digits = 2),")"))
data$unit_id<-factor(data$unit_id, levels = data%>%arrange(cor)%>%pull(unit_id)%>%unique())
p1<-ggplot(data, aes(x=temp, y=scaled_mean_YrA))+
  xy_point2(alpha=0.5)+
  regression_line()+
  ms+
  labs(y="Normalized + binned (mean) YrA", title="Cell ID", x="Core body temperature (Deg. C)")+
  theme(text = element_text(size=24))+
  facet_wrap(vars(unit_id), axes="all")#+theme(strip.text.x = element_blank())
p1
save_png_large("Normalized and binned YrA by body temperature", w=32,h=18)

