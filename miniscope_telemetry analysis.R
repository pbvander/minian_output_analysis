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
source("C:/Users/General Correa Lab/Documents/GitHub/minian_output_analysis/functions.R")

##### Things to set manually:
#Paths/directories:
exp_direc<-"C:/Users/General Correa Lab/Box/correalab/Member Folders/Paul Vander/Experiments"
output_dir<-"C:/Users/General Correa Lab/Box/correalab/Member Folders/Paul Vander/Data/Torpor project cross-experiment analyses/Miniscope"
direcs<-c("250417_circulating_E2_torpor_miniscope/pre-OVX_torpor",
          "250417_circulating_E2_torpor_miniscope/post-OVX_torpor",
          "251013_circulating_E2_torpor_miniscope/pre-ovx_torpor",
          "251013_circulating_E2_torpor_miniscope/post-ovx_torpor")
dir_struct<-c("test","mouse","start_date","session","start_time","camera") #set this to the levels of metadata stored in the directory structure of the data. Set in order from the first/higher/root directory level, to the last/lowest/final branch directory level
separator<-"/" #set depending on OS

#Exclusions:
session_id_type_to_exclude<-c("MT29_2025_06_11_session1_heat" #All frames from F0 baseline were excluded due to excess motion
                              )

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
         
         )%>%
  mutate(unit_id_id = paste0(mouse,"_",start_date,"_",session,"_",unit_id),
         session_id = paste0(mouse,"_",start_date,"_",session))%>%
  scale_temporal()

A<-A%>%mutate(unit_id_id = paste0(mouse,"_",start_date,"_",session,"_",unit_id),
              session_id = paste0(mouse,"_",start_date,"_",session))

# Add metadata on session type
read<-c()
session_mdf<-tibble()
for (dir in direcs){
  file<-paste(strsplit(dir[1],separator)[[1]][1],"session_type_metadata.csv", sep=separator)
  if (file %in% read){next}
  print(file)
  session_mdf<-rbind(session_mdf, read_csv(file, show_col_types=F)%>%mutate(across(everything(), as.character)))
  read<-c(read, file)
}
df<-df%>%merge(session_mdf, all.x=T)%>%mutate(session_id_type = paste0(session_id,"_",session_type))%>%filter(session_id_type %nin% session_id_type_to_exclude)

# Add event-based metadata
event_mdf<-tibble()
read<-c()
for (dir in direcs){
  file<-paste(strsplit(dir[1],separator)[[1]][1],"event_metadata.csv", sep=separator)
  if (file %in% read){next}
  print(file)
  event_mdf<-rbind(event_mdf, read_event_metadata(file))
  read<-c(read,file)
}
event_mdf<-event_mdf%>%select(event, session_id, frame)%>%pivot_wider(names_from = event,values_from = frame)
df<-df%>%
  merge(event_mdf, all.x=T)%>%
  mutate("male_interaction" = case_when(session_type != "male_interaction" ~ NA,
                                        session_type == "male_interaction" ~ ifelse(frame%>%between(male_added,male_removed), 1, 0)))

#Create new column for total time within each session
df<-df%>%
  group_by(session_id)%>%
  mutate(session_start_ts = min(miniscope_ts))%>%
  ungroup()%>%
  mutate(session_time_minutes = (interval(session_start_ts, miniscope_ts)%>%as.period(unit = "seconds")%>%as.numeric())/60, #time within the entire "session" 
         time_seconds = interval(start_ts, miniscope_ts)%>%as.period(unit = "seconds")%>%as.numeric()) #time in seconds within each "start_time" 

# Clean up variables and memory
rm(C,S,YrA,ts,motion)
gc()

##### Read in telemetry data
read<-c()
t_df<-tibble()
for (dir in direcs){
  telem_file<-paste(strsplit(dir,separator)[[1]][1],"telem data","telem data combined.csv", sep=separator)
  metadata_file<-paste(strsplit(dir,separator)[[1]][1],"telem data","metadata.csv", sep=separator)
  if (telem_file %in% read){next}
  print(telem_file)
  t_df<-rbind(t_df, read_telemetry_data(telem_file, metadata_file, format="starr-lifesci")) ##Ignore warning about additional pieces
  read <- c(read, telem_file)
}

##### Read in ambient temperature data
read<-c()
at_df<-tibble()
for (dir in direcs){
  at_file<-paste(strsplit(dir,separator)[[1]][1],"telem data","ambient temperature.csv", sep=separator)
  if (at_file %in% read){next}
  print(at_file)
  at_df<-rbind(at_df, read_ambient_data(at_file, upsample_interval_seconds = 60))
  read <- c(read, at_file)
}

##### Merge telemetry and miniscope data
# Create time bin column as basis for merging
ts_seq<-seq(min(t_df$telem_ts)-seconds(30),max(t_df$telem_ts)+seconds(30), seconds(60)) #assumes a 1-minute (60 second) sampling interval
df<-df%>%mutate(ts_bin = cut(miniscope_ts, ts_seq))
t_df<-t_df%>%mutate(ts_bin = cut(telem_ts, ts_seq))
at_df<-at_df%>%mutate(ts_bin = cut(ambient_ts, ts_seq))

# Preserve all miniscope data, and assign one temperature/act/ambient temperature to all frames in the given range
df<-merge(df, t_df, all=T)
df<-merge(df, at_df%>%select(!c(session_type,start_date)), all=T)
df<-df%>%mutate(ts = ifelse(is.na(miniscope_ts), telem_ts, miniscope_ts)) #set timestamp to miniscope timestamp (when available), otherwise use telem timestamp

##### Calculate dF/F0
df<-df%>%mutate(f0 = case_when(
    session_type == "cage_change" & time_seconds < 300 ~ T, #5-minute baseline before cage change
    session_type == "torpor" & aligned_time<30 & temp>35 ~ T, #day 1 torpor (euthermia timepoints)
    session_type == "torpor" & temp>35 ~ T, #day 2 torpor (euthermia timepoints)
    session_type == "cold" & time_seconds < 300 ~ T, #5-minute baseline before changing temperature
    session_type == "heat" & time_seconds < 300 ~ T, #5-minute baseline before changing temperature
    session_type == "male_interaction" & time_seconds < 300 ~ T, #5-minute baseline before adding male
    T ~ F
  )
)

f0_df<-df%>%filter(f0)%>%group_by(unit_id_id,session_id,session_type)%>%summarize(mean_f0 = mean(YrA), sd_f0 = sd(YrA))
check<-f0_df%>%group_by(session_id,session_type)%>%count()%>%nrow() - df%>%filter(!is.na(session_id))%>%group_by(session_id,session_type)%>%count()%>%nrow() 
print(check) #check that f0 is calculated for all sesssion_id values. This should equal 0
if (check !=0){stop("f0 not matching number of session_ids")}
df<-df%>%merge(f0_df, all.x=T)%>%mutate(df_f0 = (YrA - mean_f0) / mean_f0,
                                        z = (YrA - mean_f0) / sd_f0)

# Create 1-minute bins in miniscope data
sumdf<-df%>%
  filter(!is.na(YrA))%>%
  group_by(ts_bin,unit_id_id)%>%
  arrange(ts_bin)%>%
  mutate(mean_C=mean(C), mean_S=mean(S), mean_YrA=mean(YrA), mean_motion_distance=mean(motion_distance), mean_df_f0=mean(df_f0), mean_z=mean(z))%>%
  ungroup()%>%
  distinct(ts_bin,unit_id_id, .keep_all = T)%>%
  scale_temporal_bin()

##### Write output
setwd(output_dir)
write_output_rds(df)
write_output_rds(A)
write_output_rds(sumdf)

##### Checkpoint (resume here if above has run)
setwd(output_dir)
df<-read_rds("./output/df.rds")
A<-read_rds("./output/A.rds")
sumdf<-read_rds("./output/sumdf.rds")

##### Graph
gc()

for (id in df%>%filter(!is.na(session_id))%>%pull(session_id_type)%>%unique()){
  print(id)

  ### Lines with motion (quality check)
  ls<-list(scale_x_continuous(expand=c(0,0), breaks=NULL),
           theme(text=element_text(size=28),plot.title = element_text(size=28)),
           labs(x="Time (minutes)"),
           facet_wrap(vars(start_time), nrow=1,scales="free_x"),
           theme(strip.text.x = element_blank()))
  ridge_set<-list(ridgeline_guide(),
                  ridgeline(aes(height=scaled_YrA)),
                  scale_fill_viridis_c(),
                  labs(title="Detrended + demixed signal (F) / max F", y="Cell ID",x=element_blank()))
  motion_set<-list(geom_line(),
                   labs(y="pixels",x=element_blank(),title="Motion distance"))
  temp_set<-list(geom_line(linewidth = 1),
                 scale_x_continuous(expand=c(0,0),breaks=seq(0,5000,4)),
                 labs(y="Deg. C",title="Core body tempeature"))
  ambient_temp_set<-list(geom_line(linewidth = 1),
                         labs(y="Deg. C",title="Ambient tempeature",x=element_blank()))
  male_interaction_set<-list(geom_line(linewidth = 1),
                         labs(y="",title="Male social stimulus",x=element_blank()),
                         scale_y_continuous(breaks=c(0,1)))
  
  # All frames
  p1<-ggplot(df%>%filter(!is.na(YrA), session_id_type==id),aes(x=session_time_minutes,y=unit_id))+ms+ls+ridge_set
  p2<-ggplot(df%>%filter(!is.na(motion_distance), session_id_type==id), aes(x=session_time_minutes,y=motion_distance))+ms+ls+motion_set
  p3<-ggplot(df%>%filter(!is.na(YrA), session_id_type==id), aes(x=session_time_minutes, y=temp))+ms+ls+temp_set
  if (grepl("torpor",id)){
    save_png_large(paste("line plot and motion and temp",id),plot=p2/p1/p3+plot_layout(heights=c(1,10,1)),w=32,h=25)}
  if (grepl("heat",id) | grepl("cold",id)){
    p4<-ggplot(df%>%filter(!is.na(YrA), session_id_type==id), aes(x=session_time_minutes, y=ambient_temp_interpolated))+ms+ls+ambient_temp_set
    save_png_large(paste("line plot and motion and temp",id),plot=p2/p1/p4/p3+plot_layout(heights=c(1,10,1,1)),w=32,h=25)}
  if (grepl("male_interaction",id)){
    p4<-ggplot(df%>%filter(!is.na(YrA), session_id_type==id), aes(x=session_time_minutes, y=male_interaction))+ms+ls+male_interaction_set
    save_png_large(paste("line plot and motion and temp",id),plot=p2/p1/p4/p3+plot_layout(heights=c(1,10,1,1)),w=32,h=25)}
  
  # Subset of frames
  # s1<-p1+filter(p1$data, start_time=="01_55_00"|start_time=="02_56_59")
  # s2<-p2+filter(p1$data, start_time=="01_55_00"|start_time=="02_56_59")
  # s3<-p3+filter(p1$data, start_time=="01_55_00"|start_time=="02_56_59")
  # s2/s1/s3+plot_layout(heights=c(1,10,1))
  # save_png_large(paste("line plot and motion and temp subset",session_id_type),plot=s2/s1/s3+plot_layout(heights=c(1,10,1)),w=18,h=25)
}

### dF/F0 - body temperature relationship
## By temperature value
data<-sumdf%>%filter(!is.na(mean_df_f0), session_type=="torpor")
data<-merge(data, data%>%group_by(unit_id_id)%>%summarise(cor = cor(temp, mean_df_f0, method = "pearson")))%>%mutate(unit_id_id_cor = paste0(unit_id_id," (",round(cor,digits = 2),")"))
data$unit_id_id<-factor(data$unit_id_id, levels = data%>%ungroup()%>%arrange(cor)%>%distinct(unit_id_id,cor)%>%pull(unit_id_id))

p<-ggplot(data, aes(x=temp, y=mean_df_f0))+
  xy_point2(alpha=0.5)+
  regression_line()+
  ms+
  labs(y="dF/F0", title="Cell ID", x="Core body temperature (Deg. C)")+
  theme(text = element_text(size=24))+
  facet_wrap(vars(unit_id_id), axes="all",scales="free_y")#+theme(strip.text.x = element_blank())
p
save_png_large("df_f0 by body temperature", w=40,h=25)

## Torpor vs non-torpor bins
# Examine torpor status labeling
p<-ggplot(t_df%>%filter(aligned_time%>%between(-6,51)),aes(x=aligned_time,y=temp))+
  continuous_line(aes(group=mouse,color=torpor_status))+
  facet_wrap(vars(mouse,trial),axes="all")+
  scale_color_manual(values=colors)+
  ms
p
save_plot("torpor status labels",w=18,h=12)

#Non-torpor vs deep torpor
set<-list(geom_violin(aes(fill=torpor_status)),
              point_indiv(alpha=0.25,size=2, position=position_jitter(width=0.25,height=0,seed=123)),
              scale_fill_manual(values=colors),
              facet_wrap(vars(unit_id_id),scales="free_y",axes="all"),
              labs(x=element_blank(),y="dF/F0"),
              scale_x_discrete(breaks=c()))

p<-ggplot(sumdf%>%filter(!is.na(mean_df_f0), torpor_status %in% c("deep_torpor","non-torpor")), aes(x=torpor_status,y=mean_df_f0))+ms+set
p
save_plot("df_f0 non-torpor vs deep torpor",w=20,h=15)

#All torpor status bins
p<-ggplot(sumdf%>%filter(!is.na(mean_df_f0)), aes(x=torpor_status,y=mean_df_f0))+ms+set
p
save_plot("df_f0 by torpor status",w=20,h=15)

### dF/F0 - ambient temperature relationship
## By temeprature
data<-sumdf%>%filter(!is.na(mean_df_f0), !is.na(ambient_temp_interpolated))
data<-merge(data, data%>%group_by(unit_id_id)%>%summarise(cor = cor(ambient_temp_interpolated, mean_df_f0, method = "pearson")))%>%mutate(unit_id_id_cor = paste0(unit_id_id," (",round(cor,digits = 2),")"))
data$unit_id_id<-factor(data$unit_id_id, levels = data%>%ungroup()%>%arrange(cor)%>%distinct(unit_id_id,cor)%>%pull(unit_id_id))
p<-ggplot(data, aes(x=ambient_temp_interpolated, y=mean_df_f0))+
  xy_point2(alpha=0.5)+
  regression_line()+
  ms+
  labs(y="dF/F0", title="Cell ID", x="Ambient temperature (Deg. C)")+
  theme(text = element_text(size=24))+
  facet_wrap(vars(unit_id_id), axes="all",scales="free_y")#+theme(strip.text.x = element_blank())
p
save_png_large("df_f0 by ambient temperature", w=40,h=25)

## Temperature bins
set<-list(geom_violin(aes(fill=ambient_temp_bin)),
          point_indiv(alpha=0.25,size=2, position=position_jitter(width=0.25,height=0,seed=123)),
          scale_fill_manual(values=c(colors[3],colors[1],colors[2])),
          facet_wrap(vars(unit_id_id),scales="free_y",axes="all"),
          labs(x=element_blank(),y="dF/F0"),
          scale_x_discrete(breaks=c()))

sumdf<-sumdf%>%mutate(ambient_temp_bin=factor(ambient_temp_bin,levels=c("4C","22C","38C")))
p<-ggplot(sumdf%>%filter(!is.na(mean_df_f0),ambient_temp_bin!="NA"), aes(x=ambient_temp_bin, y=mean_df_f0))+ms+set
p
save_plot("df_f0 by ambient temperature bin",w=20,h=15)

### dF/F0 - male social stimulus relationship
p<-ggplot(sumdf%>%filter(!is.na(male_interaction))%>%mutate(male_interaction=factor(male_interaction)), aes(x=male_interaction,y=mean_df_f0))+
  geom_violin(aes(fill=male_interaction))+
  point_indiv(alpha=0.25,size=2, position=position_jitter(width=0.25,height=0,seed=123))+
  scale_fill_manual(values=c(colors[1],colors[5]))+
  facet_wrap(vars(unit_id_id),axes="all",scales="free_y")+
  ms
p
save_plot("df_f0 by male interaction",w=20,h=15)
