#####This code is written by PV. It is meant to read in the output from PV_minian package, integrate this data with telemetry measurements and annotated events, and visualize the data all together.
##Load packages ----
library(tidyverse) 
library(ggpubr)
library(ggprism)
library(ggridges)
library(patchwork)
library(rjson)
library(nlme)
library(gtools)
library(corrr)
library(pROC)
library(rstatix)
library(sperrorest)
library(ggnewscale)
library(GGally)
library(ggforce)
library(spatstat.geom)
# summarize<-dplyr::summarize
source("C:/Users/paulv/Documents/GitHub/minian_output_analysis/functions.R")
# filter<-dplyr::filter
options(dplyr.summarise.inform = FALSE)

##### Things to set manually:
#Paths/directories:
exp_direc<-"C:/Users/paulv/Box/correalab/Member Folders/Paul Vander/Experiments"
output_dir<-"C:/Users/paulv/Box/correalab/Member Folders/Paul Vander/Data/Torpor project cross-experiment analyses/Miniscope"
direcs<-c("250417_circulating_E2_torpor_miniscope/pre-OVX_torpor",
          "250417_circulating_E2_torpor_miniscope/post-OVX_torpor",
          "251013_circulating_E2_torpor_miniscope/pre-ovx_torpor",
          "251013_circulating_E2_torpor_miniscope/post-ovx_torpor",
          "260108_circulating_E2_torpor_miniscope/pre-ovx_torpor",
          "260108_circulating_E2_torpor_miniscope/post-ovx_torpor")
dir_struct<-c("test","mouse","start_date","session","start_time","camera") #set this to the levels of metadata stored in the directory structure of the data. Set in order from the first/higher/root directory level, to the last/lowest/final branch directory level
separator<-"/" #set depending on OS
single_session<-"MT29_2025_05_22_session1_torpor" #best session for creating figures

#Exclusions:
session_id_type_to_exclude<-c("MT29_2025_05_23_session1_heat",  "MT34_2025_12_20_session1_cold", "MT35_2026_03_07_session1_cold", # These animals do not have data from the full ambient temperature challenge (some issue arose during heat and/or cold exposure that led to exclusion of one "arm" of the ambient temperature challenge). They have been excluded so that all mice represented in the analysis have the same stimuli exposure
                              "MT30_2025_05_23_session1_torpor", "MT34_2025_11_25_session1_torpor" #sessions where torpor was very shallow and/or infrequent (min temp >33C AND <10 timepoints where temp <34). Can re-include these if cross-registration is implemented, since there was torpor on other torpor recording day for these
)

#Parameters
shuffle_iterations<-200

#Misc
lon<-7 #clock time in hours at ZT0 when lights come on
ftime<-10 #clock time in hours when fasting was started

#Settings for tuning analyses
target_cols<-c("temp_cor_torpor", 
               "ambient_temp_interpolated_cor_ambient","male_interaction_auc")
target_cols_binary<-c("temp_cor_sig_torpor",
                      "ambient_temp_interpolated_cor_sig_ambient","male_interaction_auc_sig")
labs<-c("T-Core",
        "T-Amb", "Social")

#Global graph settings:
ms<-list(theme_prism(),
         theme(text=element_text(size=12),
               strip.text = element_text(size=12,face="bold")))
theme_pie <- theme(axis.line=element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),axis.ticks = element_blank(),axis.text.x = element_blank(), legend.title = element_blank())
lw<-0.5
ls<-list(scale_x_continuous(expand=c(0,0), breaks=NULL),
         theme(text=element_text(size=12),plot.title = element_text(size=12,margin=margin(t=2,b=2,l=0,r=0,unit="pt")),panel.spacing=unit(0.06,"inches")),
         labs(x="Time (minutes)"),
         facet_wrap(vars(start_time), nrow=1, scales="free_x", space="free_x"),
         theme(strip.text.x = element_blank()))
ridge_set<-list(ridgeline_guide(),
                ridgeline(aes(height=scaled_YrA)),
                scale_fill_viridis_c(),
                labs(title="Cells (F / max F)", y=element_blank(),x=element_blank()),
                scale_y_discrete(breaks=seq(0,1000,5)))
motion_set<-list(geom_line(),
                 labs(y="pixels",x=element_blank(),title="Motion distance"))
temp_set<-list(geom_line(linewidth = lw),
               scale_x_continuous(expand=c(0,0),breaks=seq(0,5000,5)),
               labs(y=element_blank(),title="T-Core (Deg. C)"))
ambient_temp_set<-list(geom_line(linewidth = lw),
                       labs(y=element_blank(),title="T-Amb (Deg. C)",x=element_blank()))
male_interaction_set<-list(geom_line(linewidth = lw),
                           labs(y="",title="Male social stimulus",x=element_blank()),
                           scale_y_continuous(breaks=c(0,1)))
fed_set<-list(geom_line(linewidth=lw),
              labs(y="",x=element_blank(),title="Fed/fasted status"))
injection_set<-list(geom_line(linewidth=lw),
                    labs(y="",title="Injection"))
group_gonad_scale<-c("black","#56B4E9","black","#D55E00")
pellet_scale<-c("black","#56B4E9","#D55E00")
post_ovx_scale<-c("#56B4E9","#D55E00")
post_ovx_scale2<-c("#56B4E9","#D55E00")
cell_type_scale<-c("grey50","#F0E442","#0072B2")
cell_type_scale2<-c("black","#F0E442","#0072B2")
ambient_temp_bin_scale<-c("#0072B2","black","#E69F00")
male_interaction_scale<-c("black","#F0E442")

##### Read in miniscope data:
### Initalize variables
bad_frames<-list()
sumdf<-tibble()
A_all<-tibble()
event_df<-tibble()

### Read and prepare metadata
setwd(exp_direc)

#Session type metadata
read<-c()
session_mdf<-tibble()
for (dir in direcs){
  file<-paste(strsplit(dir[1],separator)[[1]][1],"session_type_metadata.csv", sep=separator)
  if (file %in% read){next}
  print(file)
  session_mdf<-rbind(session_mdf, read_csv(file, show_col_types=F)%>%mutate(across(everything(), as.character)))
  read<-c(read, file)
}

#Event-based metadata
read<-c()
event_mdf<-tibble()
for (dir in direcs){
  file<-paste(strsplit(dir[1],separator)[[1]][1],"event_metadata.csv", sep=separator)
  if (file %in% read){next}
  print(file)
  event_mdf<-rbind(event_mdf, read_event_metadata(file))
  read<-c(read,file)
}
event_metadata<-event_mdf%>%select(event, session_id,start_time, event_ts)
event_mdf<-event_metadata%>%pivot_wider(names_from = event,values_from = event_ts)

#Telemetry data
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
ts_seq<-seq(min(t_df$telem_ts)-seconds(30),max(t_df$telem_ts)+seconds(30), seconds(60)) #assumes a 1-minute (60 second) sampling interval
t_df<-t_df%>%mutate(ts_bin = cut(telem_ts, ts_seq))

#Ambient temperature data
read<-c()
at_df<-tibble()
for (dir in direcs){
  at_file<-paste(strsplit(dir,separator)[[1]][1],"telem data","ambient temperature.csv", sep=separator)
  if (at_file %in% read){next}
  print(at_file)
  at_df<-rbind(at_df, read_ambient_data(at_file, upsample_interval_seconds = 60))
  read <- c(read, at_file)
}
at_df<-at_df%>%mutate(ts_bin = cut(ambient_ts, ts_seq))

### Read in Miniscope data from each session and process individually ----
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
            # #Read data
            path <- paste(dir,mouse,start_date,session,start_time,sep=separator)
            msrun_dir <- grep("msRun",list.files(path),value = T)
            print(paste0("Reading ", path))
            exp<-strsplit(strsplit(dir,separator)[[1]][1],"_")[[1]][1]
            bad_cells<-read_cell_label(paste(path,"minian","cell_label.csv", sep=separator))
            A<-read_csv_minian(paste(path,"minian","A.csv", sep=separator))%>%filter(A>0)
            C<-read_csv_minian(paste(path,"minian","C.csv", sep=separator))
            S<-read_csv_minian(paste(path,"minian","S.csv", sep=separator))
            YrA<-read_csv_minian(paste(path,"minian","YrA.csv", sep=separator))
            motion<-read_motion(paste(path,msrun_dir,sep=separator))
            ts<-read_timestamps(paste(path,"My_V4_Miniscope","timeStamps.csv", sep=separator))
            bad_frames[[path]]<-setdiff(unique(motion$frame), unique(C$frame)) ##Detects frame in motion, but not C. This should be equal to bad_frames set in minian
            motion<-motion%>%filter(frame %nin% bad_frames[[path]])
            
            #Convert to one dataframe
            df<-merge(C,S)%>%
              merge(YrA)%>%
              merge(ts)%>%
              merge(motion)%>%
              mutate(start_ts = ymd_hms(paste(start_date, start_time)),
                     miniscope_ts = (start_ts + milliseconds(time_ms)))%>%
              mutate(unit_id_id = paste0(mouse,"_",start_date,"_",session,"_",unit_id),
                     session_id = paste0(mouse,"_",start_date,"_",session))%>%
              scale_temporal()
            
            A<-A%>%mutate(unit_id_id = paste0(mouse,"_",start_date,"_",session,"_",unit_id),
                          session_id = paste0(mouse,"_",start_date,"_",session))
            
            # Clean up variables and memory
            rm(C,S,YrA,ts,motion)
            gc()
            
            #Write intermediate output (Checkpoint 1)
            setwd(output_dir)
            if ("./output/intermediate" %nin% list.dirs()){dir.create("./output/intermediate")} #create directory if it doesn't exist
            write_output_rds(df, direc="./output/intermediate/", name=paste0("1_",gsub(separator,"-",path)))
            setwd(exp_direc)
            
            # Add metadata on session type
            df<-df%>%merge(session_mdf, all.x=T)%>%mutate(session_id_type = paste0(session_id,"_",session_type))%>%filter(session_id_type %nin% session_id_type_to_exclude)
            if (nrow(df)>0){  #skip the rest of this code if the entire session_id is excluded
              
              # Add event-based metadata
              df<-df%>%merge(event_mdf, all.x=T)
              df<-df%>%mutate("male_interaction" = case_when(session_type != "male_interaction" ~ NA,
                                                             session_type == "male_interaction" ~ ifelse(miniscope_ts > male_added & miniscope_ts < male_removed, 1, 0)),
                              # "first_interaction_aligned_time_seconds" = case_when (session_type != "male_interaction" ~ NA,
                              # session_type == "male_interaction" ~ as.duration(miniscope_ts - first_interaction)%>%as.numeric()),
                              # "first_interact" = case_when(session_type != "male_interaction" ~ NA,
                              # session_type == "male_interaction" ~ ifelse(first_interaction_aligned_time_seconds < 1, T, F)),
                              "injection" = case_when(session_type != "E2_injection" ~ NA,
                                                      session_type == "E2_injection" & miniscope_ts >= E2_injection ~ "E2",
                                                      session_type == "E2_injection" & miniscope_ts >= veh_injection ~ "Veh",
                                                      session_type == "E2_injection" & miniscope_ts >= sham_injection ~ "Sham",
                                                      session_type == "E2_injection" & miniscope_ts < sham_injection ~ "Baseline")%>%factor(levels=c("Baseline","Sham","Veh","E2")),
                              "injection_aligned_time_seconds" = case_when(session_type != "E2_injection" ~ NA,
                                                                           session_type == "E2_injection" & injection=="E2" ~ as.duration(miniscope_ts - E2_injection)%>%as.numeric(),
                                                                           session_type == "E2_injection" & injection=="Veh" ~ as.duration(miniscope_ts - veh_injection)%>%as.numeric(),
                                                                           session_type == "E2_injection" & injection=="Sham" ~ as.duration(miniscope_ts - sham_injection)%>%as.numeric(),
                                                                           session_type == "E2_injection" & injection=="Baseline" ~ as.duration(miniscope_ts - start_ts)%>%as.numeric()),
                              "fed_status" = case_when(session_type != "torpor" ~ NA,
                                                       session_type == "torpor" & miniscope_ts < fed ~ "Fasting",
                                                       session_type == "torpor" & miniscope_ts > first_bite ~ "Eating",
                                                       session_type == "torpor" & miniscope_ts >= fed & miniscope_ts <= first_bite ~ "Fed",
                                                       T ~ "Fasting")%>%factor(levels=c("Fasting","Fed","Eating"),),
                              "fed_aligned_time_seconds" = case_when(session_type != "torpor" ~ NA,
                                                                     session_type == "torpor" ~ as.duration(miniscope_ts - fed)%>%as.numeric()),
                              "cage_change_status" = case_when(session_type != "cage_change" ~ NA,
                                                               session_type == "cage_change" ~ ifelse(miniscope_ts < cage_change, "Home","New")))
              
              #Create new column for total time within each session
              df<-df%>%
                group_by(session_id)%>%
                mutate(session_start_ts = min(miniscope_ts))%>%
                ungroup()%>%
                mutate(session_time_minutes = (interval(session_start_ts, miniscope_ts)%>%as.period(unit = "seconds")%>%as.numeric())/60, #time within the entire "session"
                       time_seconds = interval(start_ts, miniscope_ts)%>%as.period(unit = "seconds")%>%as.numeric()) #time in seconds within each "start_time"
              
              ##### Merge telemetry and miniscope data
              # Create time bin column as basis for merging
              df<-df%>%mutate(ts_bin = cut(miniscope_ts, ts_seq))
              
              # Preserve all miniscope data, and assign one temperature/act/ambient temperature to all frames in the given range
              #Filter data to match mouse and date (avoid duplicates for later)
              mous<-mouse
              start_d<-start_date
              t_d<-t_df%>%filter(mouse==mous, mdy(date)==ymd(start_d))
              at_d<-at_df%>%filter(mouse==mous, start_date==start_d)
              
              #Merge
              df<-merge(df, t_d, all=T)
              df<-merge(df, at_d%>%select(!c(session_type,start_date)), all=T)
              if (df%>%filter(session_type %in% c("cold","heat"), is.na(ambient_temp_interpolated))%>%nrow() > 0){warning("NAs present in ambient_temp_interpolated column")}
              df<-df%>%mutate(ts = ifelse(is.na(miniscope_ts), telem_ts, miniscope_ts)) #set timestamp to miniscope timestamp (when available), otherwise use telem timestamp
              
              ##### Calculate dF/F0
              df<-df%>%mutate(f0 = case_when(
                session_type == "cage_change" & time_seconds < 300 ~ T, #5-minute baseline before cage change
                session_type == "torpor" & aligned_time<30 & temp>35 ~ T, #day 1 torpor (euthermia timepoints)
                session_type == "torpor" & temp>35 ~ T, #day 2 torpor (euthermia timepoints)
                session_type == "cold" & time_seconds < 300 ~ T, #5-minute baseline before changing temperature
                session_type == "heat" & time_seconds < 300 ~ T, #5-minute baseline before changing temperature
                session_type == "male_interaction" & time_seconds < 300 ~ T, #5-minute baseline before adding male
                session_type == "E2_injection" & time_seconds < 300 ~ T, #5-minute baseline before starting injections
                T ~ F)
              )
              
              f0_df<-df%>%filter(f0)%>%group_by(unit_id_id,session_id,session_type)%>%summarize(mean_f0 = mean(YrA), sd_f0 = sd(YrA))
              check<-f0_df%>%group_by(session_id,session_type)%>%count()%>%nrow() - df%>%filter(!is.na(session_id))%>%group_by(session_id,session_type)%>%count()%>%nrow()
              print(check) #check that f0 is calculated for all sesssion_id values. This should equal 0
              if (check !=0){stop("f0 not matching number of session_ids")}
              df<-df%>%merge(f0_df, all.x=T)%>%mutate(df_f0 = (YrA - mean_f0) / mean_f0,
                                                      z = (YrA - mean_f0) / sd_f0)
              
              ###Compile event_df
              setwd(output_dir)
              events<-event_metadata%>%filter(session_id == paste0(mouse,"_",start_date,"_",session), !is.na(event_ts))
              for (.event in events$event){
                print(.event)
                .ts<-events%>%filter(event==.event)%>%pull(event_ts)
                d<-df%>%
                  filter(miniscope_ts %>% between(.ts-minutes(5), .ts+minutes(5)))%>%
                  mutate(event=.event, event_ts = .ts, event_aligned_time_seconds = as.duration(miniscope_ts - event_ts)%>%as.numeric())
                if (nrow(d)>0){
                  suppressMessages(#gets rid of "scale for x is already present message
                    p<-ggplot(d, aes(x=event_aligned_time_seconds, y=unit_id))+ms+ls+ridge_set+
                      scale_x_continuous(expand=c(0,0),breaks=seq(-240,240,120))+
                      labs(title=.event,x="Event aligned time (seconds)")+
                      geom_vline(xintercept=0,linewidth=0.5)
                  )
                  save_plot(paste0(mouse,"_",start_date,"_",session,"_",.event),plot=p,w=4,h=6)
                  
                  event_df<-rbind(event_df,d)
                }else{print(paste("Skipping",.event))}
              }
              setwd(exp_direc)
              
              ##### Checkpoint 2
              setwd(output_dir)
              write_output_rds(df, direc="./output/intermediate/", name=paste0("2_",gsub(separator,"-",path)))
              setwd(exp_direc)
              
              # Create 1-minute bins in miniscope data
              sumd<-df%>%
                filter(!is.na(YrA))%>%
                group_by(ts_bin,unit_id_id)%>%
                arrange(ts_bin)%>%
                mutate(C_bin=mean(C), S_bin=mean(S), YrA_bin=mean(YrA), mean_motion_distance=mean(motion_distance), df_f0_bin=mean(df_f0), z_bin=mean(z), sd_f0_bin=mean(sd_f0))%>%
                ungroup()%>%
                distinct(ts_bin,unit_id_id, .keep_all = T)%>%
                scale_temporal_bin()
              sumdf<-rbind(sumdf,sumd)
              
              # Compile A
              A_all<-rbind(A_all,A)
              
              ## Graph full data for all sessions
              setwd(output_dir)
              for (id in df%>%filter(!is.na(session_id))%>%pull(session_id_type)%>%unique()){
                # if (paste0("line plot and motion and temp ",id,".png") %in% list.files("./output")){next}
                print(id)
                data<-df%>%filter(!is.na(YrA), session_id_type==id)
                
                ### Lines with motion (quality check)
                torpor_y<-get_torpor_y(data)
                # male_interaction_set<-c(male_interaction_set,geom_vline(xintercept = (data%>%filter(first_interact)%>%pull(session_time_minutes))[1], linewidth=lw,color="red"))
                
                # All frames
                p1<-ggplot(data, aes(x=session_time_minutes, y=unit_id))+ms+ls+ridge_set
                p2<-ggplot(data, aes(x=session_time_minutes, y=motion_distance))+ms+ls+motion_set
                suppressMessages(p3<-ggplot(data, aes(x=session_time_minutes, y=temp))+ms+ls+temp_set)
                
                if (grepl("torpor",id)){
                  p3<-p3+geom_hline(linewidth=lw, yintercept = 31, linetype="dashed")
                  aligned_tim<-interval(data$fstart[1], ymd_hms(paste0(start_date,"_00_00_01")))%>%time_length("hours")
                  if (aligned_tim < 30){ #day 1 torpor (no re-feed)
                    p<-p2/p1/p3+plot_layout(heights=c(1,10,1))}
                  if (aligned_tim > 30){  #day 2 torpor with re-feed
                    p4<-ggplot(data, aes(x=session_time_minutes, y=fed_status,group=mouse))+ms+ls+fed_set
                    p<-p2/p1/p4/p3+plot_layout(heights=c(1,10,1,1))}
                  save_png_large(paste("line plot and motion and temp",id),plot=p,w=32,h=25)}
                if (grepl("heat",id) | grepl("cold",id)){
                  p4<-ggplot(data, aes(x=session_time_minutes, y=ambient_temp_interpolated))+ms+ls+ambient_temp_set
                  save_png_large(paste("line plot and motion and temp",id),plot=p2/p1/p4/p3+plot_layout(heights=c(1,10,1,1)),w=32,h=25)}
                if (grepl("male_interaction",id)){
                  p4<-ggplot(data, aes(x=session_time_minutes, y=male_interaction))+ms+ls+male_interaction_set
                  save_png_large(paste("line plot and motion and temp",id),plot=p2/p1/p4/p3+plot_layout(heights=c(1,10,1,1)),w=32,h=25)}
                if (grepl("E2_injection", id)){
                  p4<-ggplot(data, aes(x=session_time_minutes, y=injection, group=mouse))+ms+ls+injection_set
                  save_png_large(paste("line plot and motion and temp",id),plot=p2/p1/p4/p3+plot_layout(heights=c(1,10,1,1)),w=32,h=25)}
                
                # Subset of frames
                # s1<-p1+filter(p1$data, start_time=="01_55_00"|start_time=="02_56_59")
                # s2<-p2+filter(p1$data, start_time=="01_55_00"|start_time=="02_56_59")
                # s3<-p3+filter(p1$data, start_time=="01_55_00"|start_time=="02_56_59")
                # s2/s1/s3+plot_layout(heights=c(1,10,1))
                # save_png_large(paste("line plot and motion and temp subset",session_id_type),plot=s2/s1/s3+plot_layout(heights=c(1,10,1)),w=18,h=25)
                gc()
              }
              setwd(exp_direc)
            }else{print(paste("Skipping",mouse,start_date,session))}
          }
        }
      }
    }
  }
}
###Checkpoint 3 ----
#Write
setwd(output_dir)
write_output(sumdf)
write_output(t_df)
write_output(A_all)
write_output(event_df)

#Read
setwd(output_dir)
sumdf<-read_rds("./output/sumdf.rds")
t_df<-read_rds("./output/t_df.rds")
A_all<-read_rds("./output/A_all.rds")
event_df<-read_rds("./output/event_df.rds")

##### Single-cell analysis
unit_df<-unit_analysis(sumdf%>%filter(!is.na(df_f0_bin)), roc_session_type = c("torpor","heat","cold","male_interaction"), shuf_iters=shuffle_iterations)

# Downsample to equalize temperature sampling (taking average (quantitative measures) or mode (cell type classification/"sig" columns) of downsampling iterations).
#torpor and post-ovx data only
unit_df_torpor_ovx_ds<-tibble()
pb <- txtProgressBar(min = 0, max = shuffle_iterations, style = 3)
for (i in 1:shuffle_iterations){
  d<-sumdf%>%
    filter(!is.na(df_f0_bin),session_type=="torpor",gonad=="ovx")%>%
    equalize_data_temporal(verbose=F)%>%
    unit_analysis(roc_session_type = c("torpor"), shuf_iters = shuffle_iterations, verbose=F)%>%
    mutate(iteration=i)
  unit_df_torpor_ovx_ds<-rbind(unit_df_torpor_ovx_ds,d)
  setTxtProgressBar(pb, i)
}
close(pb)
unit_df_torpor_ovx_ds_sum<-unit_df_torpor_ovx_ds%>%
  group_by(unit_id_id)%>%
  mutate(across(c(torpor_auc, torpor_fc, temp_cor_torpor, temp_slope_torpor, temp_change1_cor_torpor, temp_change1_slope_torpor), ~mean(.x)),
         across(c(torpor_auc_sig,temp_cor_sig_torpor,temp_change1_cor_sig_torpor), ~get_mode(.x)))%>%
  distinct(unit_id_id,.keep_all = T)

##### Population-level analysis
### Linear model
torpor_lm_ls<-lm_analysis(sumdf%>%filter(!is.na(df_f0_bin)), id_col="telem_ts", .session_type = "torpor", response="temp", predictor="df_f0_bin", cv_folds=5, shuf_iters=shuffle_iterations)
torpor_w_tempchange1_lm_ls<-lm_analysis(sumdf%>%filter(!is.na(df_f0_bin)), id_col="telem_ts", .session_type = "torpor", response="temp", predictor="df_f0_bin", additional_x_var = "temp_change1", cv_folds=5, shuf_iters=shuffle_iterations)
ambient_lm_ls<-lm_analysis(sumdf%>%filter(!is.na(df_f0_bin)), id_col="ambient_ts", .session_type = c("heat","cold"), response="ambient_temp_interpolated", predictor="df_f0_bin", cv_folds=5, shuf_iters=shuffle_iterations)

#Cross-training entries/arousals
sessions<-sumdf%>%group_by(session_id,torpor_status)%>%summarize(n=n_distinct(telem_ts))%>%filter(torpor_status %in% c("entry","arousal"),n>=4)%>%group_by(session_id)%>%count()%>%filter(n==2)%>%pull(session_id) #session_ids with both arousal and entry timepoints
torpor_arousal_entry_lm_ls<-lm_analysis(sumdf%>%filter((!is.na(df_f0_bin)), session_id %in% sessions, torpor_status %in% c("entry","arousal")), id_col="telem_ts", .session_type = "torpor", response="temp", predictor="df_f0_bin", cv_folds=c("entry","arousal"), partition_type="entry_arousal", partition_col="torpor_status", shuf_iters=shuffle_iterations)

#combine data
lm_df<-merge(torpor_lm_ls$lm_df, ambient_lm_ls$lm_df,all=T)%>%merge(torpor_w_tempchange1_lm_ls$lm_df,all=T)%>%merge(torpor_w_tempchange1_lm_ls$lm_add_x_var_coef_df,all=T)%>%merge(torpor_arousal_entry_lm_ls$lm_df, all=T)%>% #combine data
  merge(sumdf%>%ungroup()%>%distinct(session_id,.keep_all = T),all.x=T) #add metadata
unit_df<-merge(torpor_lm_ls$lm_coef_df, ambient_lm_ls$lm_coef_df,all=T)%>%merge(torpor_w_tempchange1_lm_ls$lm_coef_df,all=T)%>%merge(torpor_arousal_entry_lm_ls$lm_coef_df, all=T)%>% #combine data
  merge(unit_df,all=T) #add coefficients from population model to unit_df

# Gonad-intact different cell types
#torpor
data<-sumdf%>%filter(!is.na(df_f0_bin),gonad=="intact")%>%merge(unit_df%>%select(unit_id_id, session_id,mouse,temp_cor_torpor, temp_cor_sig_torpor), all.x=T)

cell_type_lm_df<-tibble()
cell_type_predict_df<-tibble()
for (cell_type in unique(data$temp_cor_sig_torpor)){
  print(cell_type)
  lm<-lm_analysis(data%>%filter(temp_cor_sig_torpor==cell_type), id_col="telem_ts", .session_type = "torpor", response="temp", predictor="df_f0_bin", cv_folds=5, shuf_iters=shuffle_iterations,verbose=F)
  cell_type_predict_df<-rbind(cell_type_predict_df,(lm$predict_df)%>%mutate(temp_cor_sig_torpor=cell_type))
  cell_type_lm_df<-rbind(cell_type_lm_df, (lm$lm_df)%>%mutate(temp_cor_sig_torpor=cell_type))
}

#ambient
data<-sumdf%>%filter(!is.na(df_f0_bin),gonad=="intact")%>%merge(unit_df%>%select(unit_id_id, session_id,mouse,ambient_temp_interpolated_cor_ambient, ambient_temp_interpolated_cor_sig_ambient), all.x=T)

ambient_cell_type_lm_df<-tibble()
ambient_cell_type_predict_df<-tibble()
for (cell_type in data%>%filter(!is.na(ambient_temp_interpolated_cor_sig_ambient))%>%pull(ambient_temp_interpolated_cor_sig_ambient)%>%unique()){
  print(cell_type)
  lm<-lm_analysis(data%>%filter(ambient_temp_interpolated_cor_sig_ambient==cell_type), id_col="telem_ts", .session_type = c("cold","heat"), response="ambient_temp_interpolated", predictor="df_f0_bin", cv_folds=5, shuf_iters=shuffle_iterations,verbose=F)
  ambient_cell_type_predict_df<-rbind(ambient_cell_type_predict_df,(lm$predict_df)%>%mutate(ambient_temp_interpolated_cor_sig_ambient=cell_type))
  ambient_cell_type_lm_df<-rbind(ambient_cell_type_lm_df, (lm$lm_df)%>%mutate(ambient_temp_interpolated_cor_sig_ambient=cell_type))
}

# Downsample to equalize temperature sampling (taking average (quantitative measures) or mode (cell type classification/"sig" columns) of downsampling iterations)
#torpor and post-ovx data only
lm_df_torpor_ovx_ds<-tibble()
lm_predict_df_torpor_ovx_ds<-tibble()
lm_coef_df_torpor_ovx_ds<-tibble()
pb <- txtProgressBar(min = 0, max = shuffle_iterations, style = 3)
for (i in 1:shuffle_iterations){
  #Run analysis
  lm_ls<-sumdf%>%filter(!is.na(df_f0_bin),gonad=="ovx")%>%
    equalize_data_temporal(verbose=F)%>%
    lm_analysis(id_col="telem_ts", .session_type = "torpor", response="temp", predictor="df_f0_bin", cv_folds=5, shuf_iters=shuffle_iterations,verbose=F)
  
  #Add to dataframes
  lm_df_torpor_ovx_ds<-rbind(lm_df_torpor_ovx_ds, (lm_ls$lm_df)%>%mutate(iteration=i))
  lm_predict_df_torpor_ovx_ds<-rbind(lm_predict_df_torpor_ovx_ds, (lm_ls$predict_df)%>%mutate(iteration=i))
  lm_coef_df_torpor_ovx_ds<-rbind(lm_coef_df_torpor_ovx_ds,(lm_ls$lm_coef_df)%>%mutate(iteration=i))
  
  #Update progress bar
  setTxtProgressBar(pb, i)
}
close(pb)
lm_df_torpor_ovx_ds_sum<-lm_df_torpor_ovx_ds%>%
  group_by(session_id)%>%
  mutate(temp_mean_cor_torpor_ = mean(temp_mean_cor_torpor_), 
         temp_cor_sig_torpor_ = get_mode(temp_cor_sig_torpor_))%>%
  distinct(session_id,.keep_all = T)%>% 
  merge(sumdf%>%ungroup()%>%distinct(session_id,.keep_all = T),all.x=T) #add metadata
lm_predict_df_torpor_ovx_ds_sum<-lm_predict_df_torpor_ovx_ds%>%
  group_by(session_id,id)%>%
  mutate(predicted = mean(predicted))%>%
  distinct(session_id,id,.keep_all = T)%>% 
  merge(sumdf%>%ungroup()%>%distinct(session_id,.keep_all = T),all.x=T) #add metadata
lm_coef_df_torpor_ovx_ds_sum<-lm_coef_df_torpor_ovx_ds%>%
  group_by(unit_id_id)%>%
  summarize(temp_MeanLmCoef_torpor_ovx_ds = mean(temp_MeanLmCoef_torpor))
unit_df<-merge(unit_df,lm_coef_df_torpor_ovx_ds_sum, all=T)

### PCA
pca_ls<-pca(sumdf%>%filter(!is.na(df_f0_bin)), predictor = "df_f0_bin", dims=c("unit_id_id","telem_ts"))
pca_time<-merge(pca_ls$telem_ts%>%select(-session_id_type)%>%mutate(telem_ts=ymd_hms(telem_ts)), sumdf%>%ungroup()%>%distinct(telem_ts,session_id,.keep_all = T), all.x=T)
pca_cell<-merge(pca_ls$unit_id_id, unit_df%>%select(-session_id,-session_id_type), all.x=T)

##### Checkpoint 3
#Write new
setwd(output_dir)
write_rds(torpor_lm_ls, "./output/torpor_lm_ls.rds")
write_rds(torpor_w_tempchange1_lm_ls, "./output/torpor_w_tempchange1_lm_ls.rds")
write_rds(ambient_lm_ls, "./output/ambient_lm_ls.rds")
write_output(lm_df)
write_output(lm_df_torpor_ovx_ds_sum)
write_output(lm_predict_df_torpor_ovx_ds_sum)
write_output(cell_type_lm_df)
write_output(cell_type_predict_df)
write_output(unit_df)
write_output(unit_df_torpor_ovx_ds_sum)
write_rds(pca_ls,"./output/pca_ls.rds")
write_output(pca_time)
write_output(pca_cell)
write_output(t_df)
write_output(A_all)
write_output(ambient_cell_type_lm_df)
write_output(ambient_cell_type_predict_df)

#Read all
setwd(output_dir)
sumdf<-read_rds("./output/sumdf.rds")
torpor_lm_ls<-read_rds("./output/torpor_lm_ls.rds")
torpor_w_tempchange1_lm_ls<-read_rds("./output/torpor_w_tempchange1_lm_ls.rds")
ambient_lm_ls<-read_rds("./output/ambient_lm_ls.rds")
lm_df<-read_rds("./output/lm_df.rds")
lm_df_torpor_ovx_ds_sum<-read_rds("./output/lm_df_torpor_ovx_ds_sum.rds")
lm_predict_df_torpor_ovx_ds_sum<-read_rds("./output/lm_predict_df_torpor_ovx_ds_sum.rds")
cell_type_lm_df<-read_rds("./output/cell_type_lm_df.rds")
cell_type_predict_df<-read_rds("./output/cell_type_predict_df.rds")
unit_df<-read_rds("./output/unit_df.rds")
unit_df_torpor_ovx_ds_sum<-read_rds("./output/unit_df_torpor_ovx_ds_sum.rds")
pca_ls<-read_rds("./output/pca_ls.rds")
pca_time<-read_rds("./output/pca_time.rds")
pca_cell<-read_rds("./output/pca_cell.rds")
t_df<-read_rds("./output/t_df.rds")
A_all<-read_rds("./output/A_all.rds")
event_df<-read_rds("./output/event_df.rds")
ambient_cell_type_lm_df<-read_rds("./output/ambient_cell_type_lm_df.rds")
ambient_cell_type_predict_df<-read_rds("./output/ambient_cell_type_predict_df.rds")

########### Graph ###########
### Session timing schematic
record_df<-tibble("start"=c(seq(15,26,1), seq(39,44,1)))%>%mutate("stop"=start+0.167)%>%rbind(tibble("start"=c(45.1,49,50.5, 53),"stop"=c(45.8,50,51.5,53.5)))
fddf<-data.frame("xmax"=c((lon-ftime)+(-5:30*24)),"xmin"=c((lon-ftime)+(-5:30*24)-12))###use this df below to create automatic shading below 

p<-ggplot(tibble("x"=seq(-2,72), "y"=1),aes(x))+ms+  
  scale_x_continuous(breaks=seq(0,60,1),#breaks=c(0, seq(15,26,1), seq(39,45,1),seq(49,52,1)),
                     labels = function(x) ifelse(x %% 3 == 0, x, "") ,
                     expand = c(0.01,0.01),
                     name="Hours")+
  scale_y_continuous(expand=c(0,0),breaks=NULL)+
  geom_rect(inherit.aes=F, data=fddf, aes(xmin=xmin,xmax=xmax),ymin=0,ymax=2,fill="grey")+
  geom_rect(inherit.aes=F, data=record_df, aes(xmin=start,xmax=stop),ymin=0,ymax=2,fill="red",alpha=0.5)+
  geom_vline(xintercept = c(14.75, 26.25, 38.75), linewidth=0.8,linetype="22")+
  theme(axis.line.y = element_blank(),
        axis.title.x = element_text(margin=margin(t=2,b=0,l=0,r=0,unit="pt")))
p+coord_cartesian(xlim=c(0,54))
save_plot("torpor session timing schematic one plot",w=12,h=1.1)
p+coord_cartesian(xlim=c(-1,0.9)) + 
  p+coord_cartesian(xlim=c(14,27)) + 
  p+coord_cartesian(xlim=c(38,47)) + 
  p+coord_cartesian(xlim=c(48,52.5)) + 
  p+coord_cartesian(xlim=c(53,53.5)) +
  plot_layout(widths=c(2,13,9,4.5,0.5), axis_titles = "collect_x")
save_plot("torpor session timing schematic split plots", w=12,h=1.1)

### Number of cells per group/session
counts<-sumdf%>%filter(!is.na(df_f0_bin))%>%ungroup()%>%distinct(unit_id_id,.keep_all = T)%>%group_by(session_id,pellet)%>%count()
t_test(counts%>%ungroup()%>%mutate(pellet=as.character(pellet))%>%filter(pellet!="pre-OVX"), n~pellet)

p<-ggplot(counts, aes(x=pellet,y=n))+
  point_errorbar()+
  point_summary(aes(color=pellet,shape=pellet))+
  point_indiv()+
  scale_color_manual(values=pellet_scale)+
  labs(x=element_blank(),y="Number of cells per session")+
  ms
p

### dF/F0 - body temperature relationship (single unit analysis) ----
### Number of observations
##Set up data
#All data
data<-sumdf%>%filter(session_type=="torpor",gonad=="ovx")%>%mutate(temp_bin1=cut(temp,seq(min(sumdf$temp)%>%floor(),max(sumdf$temp)%>%ceiling(),1)))
timepoints_per_pellet<-data%>%ungroup()%>%distinct(telem_ts, mouse,.keep_all = T)%>%mutate(pellet=droplevels(pellet))%>%group_by(pellet,temp_bin1,.drop = F)%>%summarize(timepoints_per_pellet=n())
mice_per_pellet<-data%>%ungroup()%>%distinct(mouse,pellet,.keep_all = T)%>%group_by(pellet)%>%summarize(mice_per_pellet=n())
timepoints_per_pellet_per_mouse<-merge(timepoints_per_pellet,mice_per_pellet,all=T)%>%mutate(timepoints_per_mouse=timepoints_per_pellet/mice_per_pellet)

#downsampled
data_downsampled<-sumdf%>%filter(session_type == "torpor",gonad=="ovx")%>%equalize_data_temporal()
timepoints_per_pellet_ds<-data_downsampled%>%ungroup()%>%distinct(telem_ts, mouse,.keep_all = T)%>%mutate(pellet=droplevels(pellet))%>%group_by(pellet,temp_bin1)%>%summarize(timepoints_per_pellet=n())
mice_per_pellet_ds<-data_downsampled%>%ungroup()%>%distinct(mouse,pellet,.keep_all = T)%>%group_by(pellet)%>%summarize(mice_per_pellet=n())
timepoints_per_pellet_per_mouse_ds<-merge(timepoints_per_pellet_ds,mice_per_pellet_ds,all=T)%>%mutate(timepoints_per_mouse=timepoints_per_pellet/mice_per_pellet)

#Numbers of rows in data (# of cells x # of timepoints)
p<-ggplot(data, aes(x=temp))+
  geom_histogram(aes(fill=pellet),position = position_dodge(),breaks=seq(min(sumdf$temp)%>%floor(), max(sumdf$temp)%>%ceiling(), 1))+
  scale_x_continuous(expand=c(0,0),breaks=seq(20,50,1))+
  scale_y_continuous(expand=c(0,0),limits=c(0,max(data%>%group_by(pellet,temp_bin1)%>%count()%>%pull(n))*1.1))+
  labs(x="Core temperature",y="Observations")+
  scale_fill_manual(values=post_ovx_scale)+
  ms
p
save_plot("torpor observations by temp and pellet", w=7,h=5)
p+data_downsampled
save_plot("torpor observations by temp and pellet downsampled", w=7,h=5)

#Number of timepoints
p<-ggplot(timepoints_per_pellet, aes(x=temp_bin1, y=timepoints_per_pellet))+
  geom_col(aes(fill=pellet),position = position_dodge())+
  scale_y_continuous(expand=c(0,0),limits=c(0,max(timepoints_per_pellet$timepoints_per_pellet)*1.1))+
  scale_x_discrete(limits=c(levels(data$temp_bin1[1]),levels(data$temp_bin1)[length(levels(data$temp_bin1))]), labels = function(x) ifelse(as.numeric(substr(x, 2, 3)) %% 2 == 0, substr(x,2,3), ""))+
  labs(x="T-Core (Deg. C)",y="Timepoints")+
  scale_fill_manual(values=post_ovx_scale)+
  ms+
  theme(legend.position = "none",plot.title = element_text(size=12,margin=margin(b=3,unit="pt")))
p
p+labs(title="Observed")
save_plot("torpor timepoint by temp and pellet", w=2.4,h=2)
p+timepoints_per_pellet_per_mouse_ds+labs(title="Downsampled")
save_plot("torpor timepoints by temp and pellet downsampled", w=2.4,h=2)
as_ggplot(get_legend(p+(p$data)%>%mutate(pellet=factor(pellet,levels=c("OVX+Veh","OVX+E2"),labels=c("OVX+Vehicle","OVX+E2")))+
                       theme(legend.position = "top",legend.direction = "vertical", legend.text = element_text(size=12,face="bold"))))
save_plot("pellet legend",w=4,h=0.5)

#Number of timepoints per mouse
p<-ggplot(timepoints_per_pellet_per_mouse, aes(x=temp_bin1, y=timepoints_per_mouse))+
  geom_col(aes(fill=pellet),position = position_dodge())+
  scale_y_continuous(expand=c(0,0),limits=c(0,max(timepoints_per_pellet_per_mouse$timepoints_per_mouse)*1.1))+
  scale_x_discrete(limits=c(levels(data$temp_bin1[1]),levels(data$temp_bin1)[length(levels(data$temp_bin1))]), labels = function(x) ifelse(as.numeric(substr(x, 2, 3)) %% 2 == 0, substr(x,2,3), ""))+
  labs(x="T-Core (Deg. C)",y="Timepoints per mouse")+
  scale_fill_manual(values=post_ovx_scale2,name="Treatment")+
  ms+
  theme(legend.position = "none",plot.title = element_text(size=12,margin=margin(b=3,unit="pt")))
p+labs(title="Observed")
save_plot("torpor timepoints per mouse by temp and pellet", w=2.4,h=2)
p+timepoints_per_pellet_per_mouse_ds+labs(title="Downsampled")
save_plot("torpor timepoints per mouse by temp and pellet downsampled", w=3.2,h=2.4)

## dF/F0 By temperature value
# Graph all sessions
for (id in sumdf%>%filter(session_type=="torpor")%>%pull(session_id)%>%unique()){
  print(id)
  data<-sumdf%>%filter(!is.na(df_f0_bin), session_type=="torpor", session_id==id)%>%
    merge(unit_df%>%select(unit_id_id,session_id,temp_cor_torpor,temp_cor_sig_torpor), all.x=T)%>%
    mutate(unit_id = factor(unit_id, levels = unit_df%>%filter(session_id==id)%>%arrange(temp_cor_torpor)%>%pull(unit_id)))
  
  p<-ggplot(data, aes(x=temp, y=df_f0_bin))+
    regression_line(aes(color=temp_cor_sig_torpor))+
    xy_point2(alpha=0.2,size=1)+
    scale_color_manual(values=c(cell_type_scale[2],cell_type_scale[1],cell_type_scale[3]))+
    scale_y_continuous(expand=c(0.2,0.2))+
    ms+
    labs(y="dF/F0", title="Cell ID", x="Core temperature (Deg. C)")+
    theme(legend.position = "none")+
    facet_wrap(vars(unit_id),scales="free_y",axes="all")#+theme(strip.text.x = element_blank())
  save_plot(paste("df_f0 by body temperature",id),plot=p, w=15,h=12)
  save_plot(paste("df_f0 by body temperature no color",id),plot=p+scale_color_manual(values=c("black","black","black")),w=15,h=12)
}
p<-ggplot(sumdf%>%filter(!is.na(df_f0_bin), session_type=="torpor"), aes(x=temp, y=df_f0_bin))+
  xy_point2(alpha=0.2)+
  regression_line()+
  ms+
  labs(y="dF/F0", title=paste0("All cells, all sessions (Slope = ",(lm(temp~df_f0_bin, sumdf%>%filter(!is.na(df_f0_bin), session_type=="torpor"))$coefficients[[2]])%>%round(digits=4),")"), x="Core temperature (Deg. C)")+
  theme(text = element_text(size=24))
save_plot("df_f0 by body temperature all cells all sessions",plot=p,w=10,h=10)

# Graph all cells on one graph summarized per temp_bin1
data<-sumdf%>%
  filter(!is.na(df_f0_bin), session_type=="torpor")%>%
  merge(unit_df%>%select(unit_id_id, temp_cor_torpor,temp_cor_sig_torpor,session_id))%>%
  mutate(temp_bin1=cut(temp,breaks=seq(0,50,2), labels = seq(0,48,2)),
         unit_id_id=factor(unit_id_id, levels=unit_df%>%arrange(desc(pellet),desc(temp_cor_torpor))%>%pull(unit_id_id)%>%unique()))

set<-list(theme(text=element_text(size=12),
                plot.title = element_text(size=12,margin=margin(t=3,b=3,l=0,r=0,unit="pt")),
                plot.margin = margin(t=0, b=0, l=0, r=0, "inches")))
rect_label_set<-list(geom_tile(aes(fill=temp_cor_sig_torpor)),
                     scale_x_discrete(expand=c(0,0)),
                     scale_fill_manual(values=c(cell_type_scale[2],cell_type_scale[1],cell_type_scale[3])),
                     theme(axis.line=element_blank(),axis.text = element_blank(),legend.position = "none",axis.ticks = element_blank(),axis.title = element_blank()))
rect_label_pellet<-list(geom_tile(aes(fill=pellet)),
                        scale_fill_manual(values=post_ovx_scale2),
                        theme(axis.line=element_blank(),axis.text = element_blank(),legend.position = "none",axis.ticks = element_blank(),axis.title = element_blank()))

labels_intact<-ggplot(data%>%filter(gonad=="intact"),aes(x="",y=unit_id_id))+ms+set+rect_label_set
labels_ovx<-ggplot(data%>%filter(gonad=="ovx"),aes(x="",y=unit_id_id))+ms+set+rect_label_set
pellet_label_ovx<-ggplot(data%>%filter(gonad=="ovx"),aes(x="",y=unit_id_id))+ms+set+rect_label_pellet

p<-ggplot(data%>%filter(gonad=="intact"), aes(x=temp_bin1, y=unit_id_id))+
  labs(x="T-Core (Deg. C)", y=element_blank())+
  geom_tile(aes(fill=scaled_YrA_bin))+
  scale_y_discrete(breaks=c(),expand=c(0,0))+
  scale_x_discrete(expand=c(0,0),breaks=seq(20,40,1),labels= ~ ifelse(as.numeric(as.character(.x)) %% 2 == 0, .x, ""))+
  scale_fill_continuous(type = "viridis", breaks = c(0, 1), labels = c("Min", "Max"),name="")+
  ms+
  theme(panel.background = element_rect(fill="black"),
        # legend.title = element_text(),
        # legend.text = element_text(size=12,face="bold"),
        axis.title.y=element_text(margin=margin(r=3,unit="pt")),
        # axis.title.x = element_text(margin=margin(t=3,l=-3,unit="pt")),
        axis.title.x = element_text(size=12),
        legend.text = element_blank(),
        legend.position = "top",
        legend.key.width = unit(0.2,"inches"),
        legend.justification = 0,
        legend.box.spacing = unit(6,"pt"),
        axis.line.y = element_blank())
p+labels_intact+plot_layout(widths = c(15,1))
save_plot("df_f0 by temperature all intact cells",w=2.3,h=5)
p+data%>%filter(gonad=="ovx")+labels_ovx+pellet_label_ovx+plot_layout(widths = c(20,1,1))
save_plot("df_f0 by temperature all ovx cells",w=2.3,h=4)

# Graph all data for one session
set<-list(theme(text=element_text(size=12),
                plot.title = element_text(size=12,margin=margin(t=3,b=3,l=0,r=0,unit="pt")),
                plot.margin = margin(t=0, b=0, l=3, r=0, "pt")))

#all data
d<-read_rds("./output/intermediate/2_250417_circulating_E2_torpor_miniscope-pre-OVX_torpor-MT29-2025_05_22-session1-concatenated.rds")%>%
  filter(!is.na(scaled_YrA))%>%
  merge(unit_df%>%select(temp_cor_sig_torpor,temp_cor_torpor, unit_id_id,session_id,unit_id), all.x=T)
data<-d%>%mutate(unit_id = factor(unit_id, levels=d%>%arrange(desc(temp_cor_torpor))%>%pull(unit_id)%>%unique()))

p1<-ggplot(data, aes(x=session_time_minutes, y=unit_id))+ms+ls+ridge_set+set+scale_y_discrete(breaks=c())
p2<-ggplot(data, aes(x=session_time_minutes, y=temp))+ms+ls+temp_set+set+scale_y_continuous(breaks=seq(24,38,7),expand = c(0.2,0.2))+geom_hline(linewidth=lw, yintercept = 31, linetype="dashed")
p3<-ggplot(data%>%ungroup()%>%distinct(unit_id,.keep_all = T),aes(x="",y=unit_id))+ms+set+rect_label_set
save_plot("example session all data",plot=p1+p3+p2+plot_layout(heights=c(10,1),widths=c(60,1),ncol=2),w=13,h=8)

#example cells from a single session
data<-d%>%filter(start_time %in% c("01_21_59","04_58_09"))
data<-data%>%mutate(unit_id = factor(unit_id, levels=data%>%arrange(desc(temp_cor_torpor))%>%pull(unit_id)%>%unique()))

# rect_data<-data%>%group_by(temp_cor_sig_torpor,start_time)%>%summarise(ymin = min(as.numeric(unit_id)) - 0.1, ymax = max(as.numeric(unit_id)) + 0.9,
#                                                             xmin = min(session_time_minutes), xmax = max(session_time_minutes),
#                                                             start_time=first(start_time))
# rect_set<-list(geom_rect(data=rect_data, inherit.aes=F, fill=NA, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, color = temp_cor_sig_torpor)),
#                scale_color_manual(values=cell_type_scale[2:3]))

p1<-ggplot(data, aes(x=session_time_minutes, y=unit_id))+ms+ls+ridge_set+set+scale_y_discrete(breaks=c())+theme(axis.line=element_blank(),legend.position = "none")
p2<-ggplot(data, aes(x=session_time_minutes, y=temp))+ms+ls+temp_set+set+scale_y_continuous(breaks=c(24,38),expand = c(0.25,0.25))
p3<-ggplot(data%>%ungroup()%>%distinct(unit_id,.keep_all = T), aes(x="",y=unit_id))+ms+set+rect_label_set
save_plot("example session selected cells and timepoints",plot=p1+p3+p2+plot_layout(heights=c(15,1),widths=c(20,1),ncol=2),w=3.1,h=4.8)
# save_plot("example session selected cells and timepoints with rectangle",plot=(p1+rect_set)/p2+plot_layout(heights=c(10,1)),w=4,h=5)

# Temp-temp_change1 correlation during torpor
p<-ggplot(sumdf%>%filter(session_type=="torpor")%>%ungroup()%>%distinct(mouse,telem_ts,.keep_all = T), aes(x=temp_change1, y=temp))+
  xy_point2(alpha=0.2)+
  regression_line()+
  ms
p
save_plot("temp-temp_change1 relationship during torpor",w=5,h=5)
p+facet_wrap(vars(session_id))
save_plot("temp-temp_change1 relationship during torpor by session_id",w=7,h=6)

temp_tempchange_lm<-lm_analysis(sumdf%>%filter(session_type=="torpor")%>%ungroup()%>%distinct(mouse,telem_ts,.keep_all = T), id_col="telem_ts", .session_type = "torpor", response="temp", predictor="temp_change1", cv_folds=5, shuf_iters=shuffle_iterations)
temp_tempchange_lm$lm_df

# Correlation coefficient by gonad/E2 state or pellet
p<-ggplot(unit_df%>%filter(!is.na(temp_cor_sig_torpor)), aes(x=group_gonad, y=temp_cor_torpor))+
  geom_violin(aes(fill=group_gonad))+
  point_indiv()+
  scale_fill_manual(values=group_gonad_scale)+
  facet_wrap(vars(temp_cor_sig_torpor),axes="all")+
  ms
p
save_plot("torpor temperature correlation coefficient by cell type and group_gonad", w=12,h=8)

for (cell_type in unit_df%>%filter(!is.na(temp_cor_sig_torpor))%>%pull(temp_cor_sig_torpor)%>%unique()){
  if (cell_type=="neutral"){next}
  print(cell_type)
  anov_ds<-anova(lme(data=unit_df_torpor_ovx_ds_sum%>%filter(gonad=="ovx",temp_cor_sig_torpor==cell_type), fixed=temp_cor_torpor ~ pellet, random=~1|mouse))
  print(anov_ds)
  anov<-anova(lme(data=unit_df%>%filter(gonad=="ovx",temp_cor_sig_torpor==cell_type), fixed=temp_cor_torpor ~ pellet, random=~1|mouse))
  print(anov)
  anov_torpor_only<-anova(lme(data=unit_df%>%filter(gonad=="ovx",temp_cor_sig_torpor_TempBelow34==cell_type), fixed=temp_cor_torpor_TempBelow34 ~ pellet, random=~1|mouse))
  print(anov_torpor_only)
}

p1<-ggplot(unit_df%>%filter(temp_cor_sig_torpor!="neutral"), aes(x=pellet, y=abs(temp_cor_torpor), fill=pellet))+
  geom_violin(aes())+
  labs(y="|r|",x=element_blank())+
  coord_cartesian(ylim=c(0,1))+
  point_cell()+
  point_mouse()+
  scale_fill_manual(values=pellet_scale)+
  ms+
  theme(legend.position = "none",
        axis.title.y = element_text(margin=margin(r=3,unit="pt")),
        plot.title = element_text(size=12,face="bold",margin=margin(t=0,b=3,l=0,r=0,unit="pt"),hjust=0))
p<-p1+facet_wrap(vars(temp_cor_sig_torpor), labeller = as_labeller(tools::toTitleCase))
p
save_plot("torpor temperature correlation coefficient by cell type and pellet", w=12,h=8)
p1+unit_df%>%filter(temp_cor_sig_torpor!="neutral",gonad=="intact")+labs(x=element_blank())+scale_x_discrete(labels=c("Negative","Positive"))+scale_fill_manual(values=cell_type_scale[2:3])+aes(x=temp_cor_sig_torpor,fill=temp_cor_sig_torpor)+scale_x_discrete(labels=c("Neg.","Pos."))
save_plot("torpor temperature correlation coefficient by cell type intact",w=1.8,h=1.6)
p1+unit_df%>%filter(temp_cor_sig_torpor_TempBelow34!="neutral")+aes(y=abs(temp_cor_torpor_TempBelow34))+facet_wrap(vars(temp_cor_sig_torpor_TempBelow34))
save_plot("torpor temperature correlation coefficient by cell type and pellet torpor timepoints only",w=12,h=8)
p1+unit_df_torpor_ovx_ds_sum%>%filter(temp_cor_sig_torpor_TempBelow34!="neutral")+aes(y=abs(temp_cor_torpor_TempBelow34))+facet_wrap(vars(temp_cor_sig_torpor_TempBelow34))+scale_fill_manual(values=post_ovx_scale)
save_plot("torpor temperature correlation coefficient by cell type and pellet ovx downsample torpor timepoints only",w=12,h=8)
p+(unit_df%>%filter(temp_cor_sig_torpor!="neutral",gonad=="ovx"))+scale_fill_manual(values=post_ovx_scale)
save_plot("torpor temperature correlation coefficient by cell type and pellet ovx",w=4,h=3)
p+(unit_df_torpor_ovx_ds_sum%>%filter(temp_cor_sig_torpor!="neutral"))+scale_fill_manual(values=post_ovx_scale)+labs(title="Treatment ns (both cell types)")+scale_x_discrete(labels=c("OVX+Vehicle","OVX+E2"),guide=guide_axis(n.dodge=2))
save_plot("torpor temperature correlation coefficient by cell type and pellet downsampled", w=3.2,h=2)

# Correlation type frequencies
#Grouped by pellet
data<-transform_data_piegraph(unit_df, animal_var = "pellet", cell_var = "temp_cor_sig_torpor")%>%
  mutate(temp_cor_sig_torpor=factor(temp_cor_sig_torpor,levels=c("neutral","negative","positive")))
data_torpor_only<-transform_data_piegraph(unit_df, animal_var = "pellet", cell_var = "temp_cor_sig_torpor_TempBelow34")%>%
  mutate(temp_cor_sig_torpor_TempBelow34=factor(temp_cor_sig_torpor_TempBelow34,levels=c("neutral","negative","positive")))%>%rename(temp_cor_sig_torpor=temp_cor_sig_torpor_TempBelow34)
data_torpor_only_ds<-transform_data_piegraph(unit_df_torpor_ovx_ds_sum, animal_var = "pellet", cell_var = "temp_cor_sig_torpor_TempBelow34")%>%
  mutate(temp_cor_sig_torpor_TempBelow34=factor(temp_cor_sig_torpor_TempBelow34,levels=c("neutral","negative","positive")))%>%rename(temp_cor_sig_torpor=temp_cor_sig_torpor_TempBelow34)
data_downsampled<-transform_data_piegraph(unit_df_torpor_ovx_ds_sum, animal_var = "pellet", cell_var = "temp_cor_sig_torpor")%>%
  mutate(temp_cor_sig_torpor=factor(temp_cor_sig_torpor,levels=c("neutral","negative","positive")))
data_group_gonad<-transform_data_piegraph(unit_df, animal_var = "group_gonad", cell_var = "temp_cor_sig_torpor")%>%
  mutate(temp_cor_sig_torpor=factor(temp_cor_sig_torpor,levels=c("neutral","negative","positive")))

chisq<-data%>%
  filter(pellet!="pre-OVX")%>%
  select(-percent)%>%
  pivot_wider(names_from = temp_cor_sig_torpor, values_from = n)%>%
  ungroup()%>%
  select(-pellet)%>%
  mutate(across(everything(), ~replace_na(.x,0)))%>%
  as.matrix()%>%
  chisq.test()

chisq_torpor_only<-data_torpor_only%>%
  filter(pellet!="pre-OVX")%>%
  select(-percent)%>%
  pivot_wider(names_from = temp_cor_sig_torpor, values_from = n)%>%
  ungroup()%>%
  select(-pellet)%>%
  mutate(across(everything(), ~replace_na(.x,0)))%>%
  as.matrix()%>%
  chisq.test()

chisq_ds<-data_downsampled%>%
  filter(pellet!="pre-OVX")%>%
  select(-percent)%>%
  pivot_wider(names_from = temp_cor_sig_torpor, values_from = n)%>%
  ungroup()%>%
  select(-pellet)%>%
  mutate(across(everything(), ~replace_na(.x,0)))%>%
  as.matrix()%>%
  chisq.test()

chisq_ds_to<-data_torpor_only_ds%>%
  filter(pellet!="pre-OVX")%>%
  select(-percent)%>%
  pivot_wider(names_from = temp_cor_sig_torpor, values_from = n)%>%
  ungroup()%>%
  select(-pellet)%>%
  mutate(across(everything(), ~replace_na(.x,0)))%>%
  as.matrix()%>%
  chisq.test()

pie<-ggplot(data, aes(x="", y=percent, fill=temp_cor_sig_torpor))+ms+theme_pie+
  theme(legend.position = "none",
        legend.title = element_text(hjust=0.5,margin=margin(t=0,b=3,l=0,r=0)),
        legend.title.position = "top",
        legend.text = element_text(size=12,face="bold"),
        plot.title = element_text(size=12,margin=margin(b=3,unit="pt"),hjust=0))+
  geom_bar(stat="identity", width=1,color="white",position = position_stack(reverse=T)) +
  scale_fill_manual(values=cell_type_scale, labels=tools::toTitleCase,name="Core body temperature (T-Core) correlation")+
  coord_polar("y", start=0)+
  geom_text(aes(label = paste0(round(percent,digits=0),"%"),color=temp_cor_sig_torpor,x=1.1),position = position_stack(vjust=0.5,reverse = T), size=5, fontface="bold")+
  scale_color_manual(values=c("black","black","black"))+
  guides(color="none")+
  facet_wrap(vars(pellet))
pie
save_plot("torpor temperature correlation types by pellet",w=8,h=5)
pie+data_group_gonad+facet_wrap(vars(group_gonad))
save_plot("torpor temperature correlation types by group_gonad",w=3,h=3)
pie+data%>%filter(pellet=="pre-OVX")+theme(strip.text = element_blank())
save_plot("torpor temperature correlation types pre-OVX",w=1.6,h=1.6)
pie+data_torpor_only
save_plot("torpor temperature correlation types by pellet torpor only",w=8,h=5)
pie+data_torpor_only_ds
save_plot("torpor temperature correlation types ovx by pellet downsample torpor only",w=8,h=5)
pie+data_downsampled%>%mutate(pellet=factor(pellet,levels=c("OVX+Veh","OVX+E2"),labels=c("OVX+Vehicle","OVX+E2")))+labs(title="Treatment ****")
save_plot("torpor temperature correlation types ovx by pellet downsampled",w=3,h=1.8)
get_legend(pie+theme(legend.position="top"))%>%as_ggplot()
save_plot("torpor cell type legend",w=4,h=2)

#Grouped by mouse and pellet
mouse_data<-transform_data_piegraph(unit_df, animal_var = c("mouse","pellet"), cell_var = "temp_cor_sig_torpor")%>%
  mutate(temp_cor_sig_torpor=factor(temp_cor_sig_torpor,levels=c("neutral","negative","positive")))

pie<-ggplot(mouse_data, aes(x="", y=percent, fill=temp_cor_sig_torpor)) +
  theme_prism()+
  theme_pie+
  geom_bar(stat="identity", width=1,color="white",position = position_stack(reverse=T)) +
  scale_fill_manual(values=cell_type_scale)+
  coord_polar("y", start=0)+
  geom_text(aes(label = paste0(round(percent,digits=0),"%"),color=temp_cor_sig_torpor,x=1.1),position = position_stack(vjust=0.5,reverse = T), size=5, fontface="bold")+
  scale_color_manual(values=c("grey90","black","black"))+
  guides(color="none")+
  facet_grid(vars(mouse),vars(pellet))
pie
save_plot("torpor temperature correlation types by mouse and pellet",w=8,h=5)

#Slope by gonad/E2 state
t_test(unit_df%>%filter(temp_cor_sig_torpor!="neutral", gonad=="ovx")%>%group_by(mouse,pellet,temp_cor_sig_torpor)%>%summarize(mean_slope=mean(temp_slope_torpor))%>%group_by(temp_cor_sig_torpor), mean_slope ~ pellet)
for (cell_type in unit_df%>%filter(!is.na(temp_cor_sig_torpor))%>%pull(temp_cor_sig_torpor)%>%unique()){
  if (cell_type=="neutral"){next}
  print(cell_type)
  anov<-anova(lme(data=unit_df%>%filter(gonad=="ovx",temp_cor_sig_torpor==cell_type), fixed=temp_slope_torpor ~ pellet, random=~1|mouse))
  print(anov)
}

t_test(unit_df_torpor_ovx_ds_sum%>%filter(temp_cor_sig_torpor!="neutral", gonad=="ovx")%>%group_by(mouse,pellet,temp_cor_sig_torpor)%>%summarize(mean_slope=mean(temp_slope_torpor))%>%group_by(temp_cor_sig_torpor), mean_slope ~ pellet)
for (cell_type in unique(unit_df_torpor_ovx_ds_sum$temp_cor_sig_torpor)){
  if (cell_type=="neutral"){next}
  print(cell_type)
  anov_ds<-anova(lme(data=unit_df_torpor_ovx_ds_sum%>%filter(gonad=="ovx",temp_cor_sig_torpor==cell_type), fixed=temp_slope_torpor ~ pellet, random=~1|mouse))
  print(anov_ds)
  anov<-anova(lme(data=unit_df%>%filter(gonad=="ovx",temp_cor_sig_torpor==cell_type), fixed=temp_slope_torpor ~ pellet, random=~1|mouse))
  print(anov)
  anov_torpor_only<-anova(lme(data=unit_df%>%filter(gonad=="ovx",temp_cor_sig_torpor_TempBelow34==cell_type), fixed=temp_slope_torpor_TempBelow34 ~ pellet, random=~1|mouse))
  print(anov_torpor_only)
  anov_torpor_only_ds<-anova(lme(data=unit_df_torpor_ovx_ds_sum%>%filter(gonad=="ovx",temp_cor_sig_torpor_TempBelow34==cell_type), fixed=temp_slope_torpor_TempBelow34 ~ pellet, random=~1|mouse))
  print(anov_torpor_only_ds)
}
(unit_df_torpor_ovx_ds_sum%>%filter(temp_cor_sig_torpor!="neutral")%>%group_by(mouse,pellet,temp_cor_sig_torpor)%>%summarize(mean_slope=mean(temp_slope_torpor)))%>%group_by(temp_cor_sig_torpor)%>%wilcox_test(mean_slope ~ pellet)%>%adjust_pvalue()


p1<-ggplot(unit_df%>%filter(temp_cor_sig_torpor!="neutral"), aes(x=pellet, y=abs(temp_slope_torpor),fill=pellet))+
  geom_violin(aes())+
  point_cell()+
  point_mouse()+
  coord_cartesian(ylim=c(0,NA))+
  labs(x=element_blank(),y="|Slope|")+
  scale_fill_manual(values=pellet_scale)+
  ms+
  theme(legend.position = "none",
        axis.title.y = element_text(margin=margin(r=5,unit="pt")),
        plot.title = element_text(size=12,hjust=0,margin=margin(b=3,unit="pt")))
p<-p1+facet_wrap(vars(temp_cor_sig_torpor),axes="all",scales="free",labeller = as_labeller(tools::toTitleCase))
p
save_plot("torpor temperature slope by cell type and pellet", w=12,h=8)
p1+aes(y=temp_slope_torpor_TempBelow34)
save_plot("torpor temperature slope by cell type and pellet torpor only",w=12,h=8)
p1+unit_df_torpor_ovx_ds_sum%>%filter(temp_cor_sig_torpor_TempBelow34!="neutral")+aes(y=abs(temp_slope_torpor_TempBelow34))+scale_fill_manual(values=post_ovx_scale)+facet_wrap(vars(temp_cor_sig_torpor_TempBelow34))
save_plot("torpor temperature slope by cell type and pellet torpor only downsample",w=12,h=8)
p1+unit_df%>%filter(temp_cor_sig_torpor!="neutral",gonad=="intact")+labs(x=element_blank())+scale_x_discrete(labels=c("Negative","Positive"))+scale_fill_manual(values=cell_type_scale[2:3])+aes(x=temp_cor_sig_torpor,fill=temp_cor_sig_torpor)+scale_x_discrete(labels=c("Neg.","Pos."))
save_plot("torpor temperature slope by cell type intact",w=1.8,h=1.6)
p+unit_df_torpor_ovx_ds_sum%>%filter(temp_cor_sig_torpor!="neutral")+scale_fill_manual(values=post_ovx_scale)+labs(title="Treatment ns (both cell types)")+scale_x_discrete(labels=c("OVX+Vehicle","OVX+E2"),guide=guide_axis(n.dodge=2))
save_plot("torpor temperature slope by cell type and pellet downsampled",w=3.2,h=2)

## Torpor vs non-torpor bins
# Examine torpor status labeling
p<-ggplot(t_df%>%filter(aligned_time%>%between(-6,51)),aes(x=aligned_time,y=temp))+
  continuous_line(aes(group=mouse,color=torpor_status))+
  facet_wrap(vars(mouse,trial),axes="all")+
  scale_color_manual(values=colors)+
  ms
p
save_plot("torpor status labels",w=18,h=12)

# Non-torpor vs deep torpor
set<-list(geom_violin(aes(fill=torpor_status)),
          point_indiv(alpha=0.25,size=2, position=position_jitter(width=0.25,height=0,seed=123)),
          scale_fill_manual(values=colors),
          facet_wrap(vars(unit_id_id),scales="free_y",axes="all"),
          labs(x=element_blank(),y="dF/F0"),
          scale_x_discrete(breaks=c()))

p<-ggplot(sumdf%>%filter(!is.na(df_f0_bin), torpor_status %in% c("deep_torpor","non-torpor")), aes(x=torpor_status,y=df_f0_bin))+ms+set
p
save_plot("df_f0 non-torpor vs deep torpor",w=40,h=30)

# All torpor status bins
p<-ggplot(sumdf%>%filter(!is.na(df_f0_bin)), aes(x=torpor_status,y=df_f0_bin))+ms+set
p
save_plot("df_f0 by torpor status",w=40,h=30)

# ROC analysis
#Cell type frequencies
#By pellet
data<-transform_data_piegraph(unit_df, animal_var="pellet", cell_var = "torpor_auc_sig")%>%
  mutate(torpor_auc_sig=factor(torpor_auc_sig,levels=c("neutral","activated","suppressed")))

data_ds<-transform_data_piegraph(unit_df_torpor_ovx_ds_sum, animal_var="pellet", cell_var = "torpor_auc_sig")%>%
  mutate(torpor_auc_sig=factor(torpor_auc_sig,levels=c("neutral","activated","suppressed")))

chisq<-data%>%
  filter(pellet!="pre-OVX")%>%
  select(-percent)%>%
  pivot_wider(names_from = torpor_auc_sig, values_from = n)%>%
  ungroup()%>%
  select(-pellet)%>%
  mutate(across(everything(), ~replace_na(.x,0)))%>%
  as.matrix()%>%
  chisq.test()

pie<-ggplot(data, aes(x="", y=percent, fill=torpor_auc_sig)) +
  theme_prism()+
  theme_pie+
  geom_bar(stat="identity", width=1,color="white",position = position_stack(reverse=T)) +
  scale_fill_manual(values=cell_type_scale)+
  coord_polar("y", start=0)+
  geom_text(aes(label = paste0(round(percent,digits=0),"%"),color=torpor_auc_sig,x=1.1),position = position_stack(vjust=0.5,reverse = T), size=5, fontface="bold")+
  scale_color_manual(values=c("grey90","black","black"))+
  guides(color="none")+
  facet_wrap(vars(pellet))
pie
save_plot("torpor roc types by pellet",w=8,h=5)
pie+data_ds
save_plot("torpor roc types by pellet downsample",w=5,h=3)

#Grouped by mouse and pellet
mouse_data<-transform_data_piegraph(unit_df, animal_var = c("mouse","pellet"), cell_var = "torpor_auc_sig")%>%
  mutate(torpor_auc_sig=factor(torpor_auc_sig,levels=c("neutral","activated","suppressed")))

pie<-ggplot(mouse_data, aes(x="", y=percent, fill=torpor_auc_sig)) +
  theme_prism()+
  theme_pie+
  geom_bar(stat="identity", width=1,color="white",position = position_stack(reverse=T)) +
  scale_fill_manual(values=cell_type_scale)+
  coord_polar("y", start=0)+
  geom_text(aes(label = paste0(round(percent,digits=0),"%"),color=torpor_auc_sig,x=1.1),position = position_stack(vjust=0.5,reverse = T), size=5, fontface="bold")+
  scale_color_manual(values=c("grey90","black","black"))+
  guides(color="none")+
  facet_grid(vars(mouse),vars(pellet))
pie
save_plot("torpor roc types by pellet and mouse",w=8,h=5)

#Fold-change between status bins
p<-ggplot(unit_df%>%filter(!is.na(torpor_auc_sig))%>%mutate(torpor_auc_sig=factor(torpor_auc_sig,levels=c("neutral","activated","suppressed"))),aes(x=torpor_auc,y=abs(log2(torpor_fc))))+
  xy_point3(aes(color=torpor_auc_sig))+
  scale_color_manual(values=cell_type_scale2)+
  ms+
  theme(legend.position = "none")
p
save_plot("torpor status bins fc",w=10,h=7)
p+unit_df_torpor_ovx_ds_sum%>%filter(!is.na(torpor_auc_sig))%>%mutate(torpor_auc_sig=factor(torpor_auc_sig,levels=c("neutral","activated","suppressed")))
save_plot("torpor status bins fc downsample",w=8,h=5)

p<-ggplot(unit_df%>%filter(!is.na(torpor_auc_sig),torpor_auc_sig!="neutral"),aes(x=pellet,y=(log2(torpor_fc))))+
  geom_violin(aes(fill=pellet))+
  point_summary(aes(color=mouse),position=position_jitter(width=0.05,height=0,seed=123))+
  point_indiv()+
  labs(x=element_blank(),y="Fold Change")+
  scale_fill_manual(values=pellet_scale)+
  facet_wrap(vars(torpor_auc_sig),axes="all",scales="free")+
  ms+
  theme(legend.position = "none")
p
save_plot("torpor roc fold change",w=5,h=4)
p+unit_df_torpor_ovx_ds_sum%>%filter(!is.na(torpor_auc_sig),torpor_auc_sig!="neutral")+scale_fill_manual(values=post_ovx_scale)
save_plot("torpor roc fold change downsample",w=6,h=4)

### dF-F0 - body temperature (population level analysis) ----
##LM significance by pellet
data<-transform_data_piegraph(lm_df,"pellet","temp_cor_sig_torpor_")
data_ovx_ds<-transform_data_piegraph(lm_df_torpor_ovx_ds_sum,"pellet","temp_cor_sig_torpor_")

pie<-ggplot(data, aes(x="", y=percent, fill=temp_cor_sig_torpor_)) +
  theme_prism()+
  theme_pie+
  geom_bar(stat="identity", width=1,color="white",position = position_stack(reverse=T)) +
  scale_fill_manual(values=cell_type_scale)+
  coord_polar("y", start=0)+
  geom_text(aes(label = paste0(round(percent,digits=0),"%"),color=temp_cor_sig_torpor_,x=1.1),position = position_stack(vjust=0.5,reverse = T), size=5, fontface="bold")+
  scale_color_manual(values=c("grey90","black","black"))+
  guides(color="none")+
  facet_grid(vars(pellet))
pie
save_plot("lm correlation significance by pellet",w=4,h=5)
pie+data_ovx_ds
save_plot("lm correlation significance by pellet downsample",w=4,h=5)

##LM accuracy
lm_anova<-anova(lme(data=lm_df%>%filter(temp_cor_sig_torpor_=="significant"),
                    fixed=temp_mean_cor_torpor_ ~ pellet,
                    random=~1|mouse))
lm_anova_ds<-anova(lme(data=lm_df_torpor_ovx_ds_sum%>%filter(temp_cor_sig_torpor_=="significant"),
                       fixed=temp_mean_cor_torpor_ ~ pellet,
                       random=~1|mouse))

p<-ggplot(lm_df%>%filter(temp_cor_sig_torpor_=="significant"), aes(x=pellet,y=temp_mean_cor_torpor_))+
  geom_violin(aes(fill=pellet))+
  point_mouse()+
  point_cell()+
  labs(x=element_blank(),y="r")+
  scale_fill_manual(values=pellet_scale)+
  ms+
  coord_cartesian(ylim=c(0,1))+
  theme(legend.position = "none",
        plot.title=element_text(size=12,hjust=0,margin=margin(b=3,unit="pt")),
        plot.subtitle = element_text(size=12,hjust=0.5,margin=margin(t=3,b=3,unit="pt")))
p
save_plot("lm correlation coefficient by pellet",w=6,h=6)
p+lm_df%>%filter(temp_cor_sig_torpor_=="significant",gonad=="intact")
save_plot("lm correlation coefficient intact",w=2,h=2)
p+lm_df_torpor_ovx_ds_sum%>%filter(temp_cor_sig_torpor_=="significant")+scale_fill_manual(values=post_ovx_scale)+labs(title="Treatment ns",subtitle="All cells")+scale_x_discrete(labels=c("OVX+Vehicle","OVX+E2"))
save_plot("lm correlation coefficient by pellet downsample", w=2.7,h=2)

data<-cell_type_lm_df%>%
  merge(lm_df%>%select(session_id,mouse),all.x=T)%>%
  rbind(lm_df%>%filter(gonad=="intact")%>%select(session_id,temp_mean_cor_torpor_,temp_cor_sig_torpor_,mouse)%>%mutate(temp_cor_sig_torpor="All"))%>%
  mutate(temp_cor_sig_torpor=factor(temp_cor_sig_torpor,levels=c("All","neutral","negative","positive"),labels=c("All","Neutral","Negative","Positive")))
p<-ggplot(data, aes(x=temp_cor_sig_torpor,y=temp_mean_cor_torpor_))+
  geom_violin(aes(fill=temp_cor_sig_torpor),scale="width",width=0.95)+
  point_summary(aes(color=mouse),position=position_jitter(width=0.15,height=0,seed=123),size=2.8,shape=21,stroke=1)+
  point_indiv(size=1.5,position=position_jitter(width=0.25,height=0,seed=321))+
  labs(x=element_blank(),y="r")+
  scale_fill_manual(values=c("black",cell_type_scale))+
  scale_y_continuous(breaks=seq(0,1,0.5),labels=seq(0,1,0.5))+
  ms+
  coord_cartesian(ylim=c(0,1))+
  theme(legend.position = "none",
        axis.title.y = element_text(margin = margin(r=3, unit="pt")))+
  scale_x_discrete(labels=c("All","Neutral","Neg.","Pos."), guide=guide_axis(n.dodge = 2))
p
save_plot("lm correlation by cell type",w=2,h=1.84)

##LM predictions
for (sess in sumdf%>%filter(session_type=="torpor")%>%pull(session_id)%>%unique()){
  data<-(torpor_lm_ls$predict_df)%>%filter(session_id==sess)
  data_ds<-lm_predict_df_torpor_ovx_ds_sum%>%filter(session_id==sess)
  data_cell_type<-cell_type_predict_df%>%
    filter(session_id==sess)%>%
    mutate(fold=as.numeric(fold))%>%
    bind_rows(data%>%mutate(temp_cor_sig_torpor="All"))%>%
    mutate(temp_cor_sig_torpor=factor(temp_cor_sig_torpor,levels=c("All","neutral","negative","positive"),labels=c("All","Neutral","Negative","Positive")))
  
  p<-ggplot(data, aes(x=predicted,y=true))+
    labs(x="Predicted", y="Observed",title="T-Core (all cells)")+
    geom_abline(slope=1,linetype="dashed",linewidth=1)+
    regression_line(linewidth=1,alpha=0.5,color="blue")+
    xy_point2(alpha=0.2,size=1)+
    ms+theme(plot.title=element_text(size=12,margin=margin(t=0,r=0,b=3,l=-15,"pt")),
             panel.spacing = unit(3,"pt"),
             axis.title.x = element_text(margin = margin(t = 3,unit="pt")),
             axis.title.y = element_text(margin = margin(r=3, unit="pt")))
  p
  save_plot(paste("predicted temp",sess,"torpor"),w=1.9,h=1.9)
  p+data_ds
  save_plot(paste("predicted temp",sess,"torpor downsampled"),w=6,h=6)
  p+data_cell_type+facet_wrap(vars(temp_cor_sig_torpor),axes="all",nrow=1)+labs(title=element_blank())
  save_plot(paste("predicted temp",sess,"torpor cell types"),w=4.7,h=1.8)
}

##LM weight
p<-ggplot(unit_df%>%filter(!is.na(temp_cor_sig_torpor)), aes(x=temp_cor_sig_torpor,y=temp_MeanLmCoef_torpor,fill=temp_cor_sig_torpor))+
  geom_violin()+
  scale_fill_manual(values=cell_type_scale)+
  ms
p+facet_wrap(vars(pellet))
save_plot("lm coefficient by cell type and pellet",w=5,h=3)

### dF/F0 - ambient temperature relationship ----
## Plot ambient temperature challenge schematic
ggplot(sumdf%>%filter(session_id=="MT35_2026_02_10_session1", session_type %in% c("cold","heat"))%>%mutate(session_type=factor(session_type,levels=c("heat","cold"))),aes(x=session_time_minutes,y=ambient_temp_interpolated))+
  labs(x="Time (minutes)",y="Ambient temperature (Deg. C)")+
  continuous_line()+
  ms+theme(panel.spacing=unit(1,"in"))+
  facet_wrap(vars(session_type),scales="free_x")
save_plot("ambient temperature schematic",w=7,h=5)

##Observations by pellet group
#Set up data
data<-sumdf%>%filter(session_type %in% c("cold","heat"),gonad=="ovx")%>%mutate(ambient_temp_bin1=cut(ambient_temp_interpolated,seq(3,39,1)))
timepoints_per_pellet<-data%>%ungroup()%>%distinct(telem_ts, mouse,.keep_all = T)%>%mutate(pellet=droplevels(pellet))%>%group_by(pellet,ambient_temp_bin1,.drop = F)%>%summarize(timepoints_per_pellet=n())
mice_per_pellet<-data%>%ungroup()%>%distinct(mouse,pellet,.keep_all = T)%>%group_by(pellet)%>%summarize(mice_per_pellet=n())
timepoints_per_pellet_per_mouse<-merge(timepoints_per_pellet,mice_per_pellet,all=T)%>%mutate(timepoints_per_mouse=timepoints_per_pellet/mice_per_pellet)

#Numbers of rows in data (# of cells x # of timepoints)
p<-ggplot(data, aes(x=ambient_temp_interpolated))+
  geom_histogram(aes(fill=pellet),position = position_dodge(),breaks=seq(3, 39, 1))+
  scale_x_continuous(expand=c(0,0),breaks=seq(3,39,1))+
  scale_y_continuous(expand=c(0,0))+
  labs(x="Ambient temperature",y="Observations")+
  scale_fill_manual(values=post_ovx_scale)+
  ms
p
save_plot("ambient temp observations by ambient_temp and pellet", w=7,h=5)

#Number of timepoints
p<-ggplot(timepoints_per_pellet, aes(x=ambient_temp_bin1, y=timepoints_per_pellet))+
  geom_col(aes(fill=pellet),position = position_dodge())+
  scale_y_continuous(expand=c(0,0),limits=c(0,max(timepoints_per_pellet$timepoints_per_pellet)*1.1))+
  scale_x_discrete(limits=c(levels(data$ambient_temp_bin1[1]),levels(data$ambient_temp_bin1)[length(levels(data$ambient_temp_bin1))]))+
  labs(x="Ambient temperature",y="Timepoints")+
  scale_fill_manual(values=post_ovx_scale)+
  ms
p
save_plot("ambient temp timepoints by temp and pellet", w=7,h=5)

#Number of timepoints per mouse
p<-ggplot(timepoints_per_pellet_per_mouse, aes(x=ambient_temp_bin1, y=timepoints_per_mouse))+
  geom_col(aes(fill=pellet),position = position_dodge())+
  scale_y_continuous(expand=c(0,0),limits=c(0,max(timepoints_per_pellet_per_mouse$timepoints_per_mouse)*1.1))+
  scale_x_discrete(limits=c(levels(data$ambient_temp_bin1[1]),levels(data$ambient_temp_bin1)[length(levels(data$ambient_temp_bin1))]))+
  labs(x="Ambient temperature",y="Timepoints per mouse")+
  scale_fill_manual(values=post_ovx_scale)+
  ms
p
save_plot("ambient temp timepoints per mouse by temp and pellet", w=7,h=5)

# Graph all cells on one graph summarized per ambient_temp_interpolated_bin1
data<-sumdf%>%
  filter(!is.na(df_f0_bin), session_type %in% c("heat","cold"))%>%
  group_by(unit_id_id)%>%
  mutate(scaled_df_f0_bin = scales::rescale(df_f0_bin))%>%
  ungroup()%>%
  merge(unit_df%>%select(unit_id_id, ambient_temp_interpolated_cor_ambient,ambient_temp_interpolated_cor_sig_ambient,session_id))%>%
  mutate(ambient_temp_interpolated_bin1=cut(ambient_temp_interpolated,breaks=c(c(3.9,6.1),seq(8,34,2),c(36.1,38.1)), labels = seq(5,37,2)),
         unit_id_id=factor(unit_id_id, levels=unit_df%>%arrange(desc(pellet),desc(ambient_temp_interpolated_cor_ambient))%>%pull(unit_id_id)%>%unique()))

set<-list(theme(text=element_text(size=12),
                plot.title = element_text(size=12,margin=margin(t=3,b=3,l=0,r=0,unit="pt")),
                plot.margin = margin(t=0, b=0, l=0, r=0, "inches")))
rect_label_set<-list(geom_tile(aes(fill=ambient_temp_interpolated_cor_sig_ambient)),
                     scale_x_discrete(expand=c(0,0)),
                     scale_fill_manual(values=c(cell_type_scale[2],cell_type_scale[1],cell_type_scale[3])),
                     theme(axis.line=element_blank(),axis.text = element_blank(),legend.position = "none",axis.ticks = element_blank(),axis.title = element_blank()))
rect_label_pellet<-list(geom_tile(aes(fill=pellet)),
                        scale_fill_manual(values=post_ovx_scale2),
                        theme(axis.line=element_blank(),axis.text = element_blank(),legend.position = "none",axis.ticks = element_blank(),axis.title = element_blank()))

labels_intact<-ggplot(data%>%filter(gonad=="intact"),aes(x="",y=unit_id_id))+ms+set+rect_label_set
labels_ovx<-ggplot(data%>%filter(gonad=="ovx"),aes(x="",y=unit_id_id))+ms+set+rect_label_set
pellet_label_ovx<-ggplot(data%>%filter(gonad=="ovx"),aes(x="",y=unit_id_id))+ms+set+rect_label_pellet

p<-ggplot(data%>%filter(gonad=="intact"), aes(x=ambient_temp_interpolated_bin1, y=unit_id_id))+
  labs(x="T-Amb (Deg. C)", y=element_blank())+
  geom_tile(aes(fill=scaled_df_f0_bin))+
  scale_y_discrete(breaks=c(),expand=c(0,0))+
  scale_x_discrete(expand=c(0,0),breaks= ~ ., labels= ~ ifelse((as.numeric(as.character(.x))-1) %% 4 == 0, .x, ""))+
  scale_fill_continuous(type = "viridis", breaks = c(0, 1), labels = c("Min", "Max"),name="")+
  ms+
  theme(panel.background = element_rect(fill="black"),
        # legend.title = element_text(),
        # legend.text = element_text(size=12,face="bold"),
        axis.title.y=element_text(margin=margin(r=3,unit="pt")),
        axis.title.x = element_text(size=12,margin=margin(t=4,l=-25,unit="pt")),
        legend.text = element_blank(),
        legend.position = "top",
        legend.key.width = unit(0.2,"inches"),
        legend.justification = 0.2,
        legend.box.spacing = unit(6,"pt"),
        axis.line.y = element_blank())
p+labels_intact+plot_layout(widths = c(15,1))
save_plot("df_f0 by ambient temperature all intact cells",w=2,h=3.8)
p+data%>%filter(gonad=="ovx")+labels_ovx+pellet_label_ovx+plot_layout(widths = c(20,1,1))
save_plot("df_f0 by ambient temperature all ovx cells",w=2.3,h=4.4)

## Plot an example session
set<-list(theme(text=element_text(size=12),
                plot.title = element_text(size=12,margin=margin(t=3,b=3,l=0,r=0,unit="pt")),
                plot.margin = margin(t=0, b=-5, l=3, r=0, "pt"),
                axis.title.y = element_text(margin(t=3,b=-5,l=0,r=0,unit="pt")),
                panel.spacing = unit(0.05,"inches")))

#all data
d<-read_rds("./output/intermediate/2_251013_circulating_E2_torpor_miniscope-pre-ovx_torpor-MT31-2025_11_23-session1-concatenated.rds")%>%
  filter(!is.na(scaled_YrA), session_type %in% c("cold","heat"))%>%
  group_by(unit_id_id)%>%
  mutate(scaled_df_f0 = scales::rescale(df_f0))%>%
  ungroup()%>%
  mutate(session_time_minutes = session_time_minutes-min(session_time_minutes))%>%
  merge(unit_df%>%select(ambient_temp_interpolated_cor_sig_ambient,ambient_temp_interpolated_cor_ambient, unit_id_id,session_id,unit_id), all.x=T)
data<-d%>%mutate(unit_id = factor(unit_id, levels=d%>%arrange(desc(ambient_temp_interpolated_cor_ambient))%>%pull(unit_id)%>%unique()))

p1<-ggplot(data, aes(x=session_time_minutes, y=unit_id))+ms+ls+ridge_set+set+scale_y_discrete(breaks=c())+aes(height=scaled_df_f0)
p2<-ggplot(data, aes(x=session_time_minutes, y=ambient_temp_interpolated))+ms+ls+ambient_temp_set+set+scale_y_continuous(expand = c(0.2,0.2))
p3<-ggplot(data, aes(x=session_time_minutes, y=temp))+ms+ls+temp_set+set+scale_y_continuous(expand = c(0.2,0.2))+labs(title="Core body temperature")
p4<-ggplot(data%>%ungroup()%>%distinct(unit_id,.keep_all = T),aes(x="",y=unit_id))+ms+set+rect_label_set
save_plot("example session all data ambient",plot=p1+p4+p2+plot_spacer()+p3+plot_layout(heights=c(12,1,1),widths=c(60,1),ncol=2),w=13,h=8)

#subset of all data
data<-d%>%
  group_by(session_type)%>%
  filter(session_time_minutes < min(session_time_minutes + 5) | session_time_minutes > max(session_time_minutes-5))%>%
  mutate(start_time=cut(session_time_minutes, breaks=c(0,6,80,130,250)))%>%
  filter(!is.na(start_time))
data<-data%>%mutate(unit_id = factor(unit_id, levels=data%>%arrange(desc(ambient_temp_interpolated_cor_ambient))%>%pull(unit_id)%>%unique()))

p1<-ggplot(data, aes(x=session_time_minutes, y=unit_id))+ms+ls+ridge_set+set+scale_y_discrete(breaks=c())+aes(height=scaled_df_f0)+theme(axis.line=element_blank(),legend.position = "none")
p2<-ggplot(data, aes(x=session_time_minutes, y=ambient_temp_interpolated))+ms+ls+ambient_temp_set+set+scale_y_continuous(expand = c(0.25,0.25),breaks=c(5,37))
p3<-ggplot(data, aes(x=session_time_minutes, y=temp))+ms+ls+temp_set+set+scale_y_continuous(expand = c(0.25,0.25),breaks=seq(34,40,4))+scale_x_continuous(expand=c(0.01,0.01),breaks=c(0,3,46,49,102,105,181,184))
p4<-ggplot(data%>%ungroup()%>%distinct(unit_id,.keep_all = T),aes(x="",y=unit_id))+ms+set+rect_label_set
save_plot("example session subset data ambient",plot=p1+p4+p2+plot_spacer()+p3+plot_layout(heights=c(7,1,1),widths=c(20,1),ncol=2),w=2.5,h=3.3)

(p4+(p4$data)%>%mutate(ambient_temp_interpolated_cor_sig_ambient=factor(ambient_temp_interpolated_cor_sig_ambient,levels=c("neutral","negative","positive")))+
    scale_fill_manual(name="Ambient temperature (T-Amb) correlation",values=cell_type_scale,labels=tools::toTitleCase)+
    theme(legend.position = "top",
          legend.text = element_text(size=12,face="bold"),
          legend.title=element_text(hjust=0.5),
          legend.title.position = "top"))%>%get_legend()%>%as_ggplot()
save_plot("ambient temperature cell type legend",w=4,h=0.5)

## By temeprature
for (id in sumdf%>%filter(session_type %in% c("cold","heat"))%>%pull(session_id)%>%unique()){
  print(id)
  data<-sumdf%>%filter(!is.na(df_f0_bin), !is.na(ambient_temp_interpolated), session_id==id)
  data<-data%>%mutate(unit_id_id = factor(unit_id_id, levels = unit_df%>%arrange(ambient_temp_interpolated_cor_ambient)%>%pull(unit_id_id)))
  p<-ggplot(data, aes(x=ambient_temp_interpolated, y=df_f0_bin))+
    xy_point2(alpha=0.5)+
    regression_line()+
    ms+
    labs(y="dF/F0", title="Cell ID", x="Ambient temperature (Deg. C)")+
    theme(text = element_text(size=24))+
    facet_wrap(vars(unit_id_id), axes="all",scales="free_y")#+theme(strip.text.x = element_blank())
  p
  save_plot(paste("df_f0 by ambient temperature",id), w=40,h=25)
  
  ## Temperature bins
  set<-list(geom_violin(aes(fill=ambient_temp_bin)),
            point_indiv(alpha=0.25,size=2, position=position_jitter(width=0.25,height=0,seed=123)),
            scale_fill_manual(values=c(ambient_temp_bin_scale)),
            facet_wrap(vars(unit_id_id),scales="free_y",axes="all"),
            labs(x=element_blank(),y="dF/F0"),
            scale_x_discrete(breaks=c()))
  
  p<-ggplot(data%>%filter(ambient_temp_bin!="NA"), aes(x=ambient_temp_bin, y=df_f0_bin))+ms+set
  p
  save_plot(paste("df_f0 by ambient temperature bin",id),w=20,h=15)
}

##Cell type frequencies
data<-transform_data_piegraph(unit_df, animal_var="pellet",cell_var="ambient_temp_interpolated_cor_sig_ambient")%>%
  mutate(ambient_temp_interpolated_cor_sig_ambient = factor(ambient_temp_interpolated_cor_sig_ambient, levels=c("neutral","negative","positive"),labels=c("Neutral","Negative","Positive")))

chisq<-data%>%
  filter(pellet!="pre-OVX")%>%
  select(-percent)%>%
  pivot_wider(names_from = ambient_temp_interpolated_cor_sig_ambient, values_from = n)%>%
  ungroup()%>%
  select(-pellet)%>%
  mutate(across(everything(), ~replace_na(.x,0)))%>%
  as.matrix()%>%
  chisq.test()

pie<-ggplot(data, aes(x="", y=percent, fill=ambient_temp_interpolated_cor_sig_ambient))+ms+theme_pie+
  theme(legend.position = "none",
        legend.title = element_text(hjust=0.5,margin=margin(t=0,b=3,l=0,r=0)),
        legend.title.position = "top",
        legend.text = element_text(size=12,face="bold"),
        plot.title=element_text(size=12,hjust=0,margin=margin(b=3,unit="pt")))+
  geom_bar(stat="identity", width=1,color="white",position = position_stack(reverse=T)) +
  scale_fill_manual(values=cell_type_scale, labels=tools::toTitleCase,name="Core body temperature (T-Core) correlation")+
  coord_polar("y", start=0)+
  geom_text(aes(label = paste0(round(percent,digits=0),"%"),color=ambient_temp_interpolated_cor_sig_ambient,x=1.1),position = position_stack(vjust=0.5,reverse = T), size=5, fontface="bold")+
  scale_color_manual(values=c("black","black","black"))+
  guides(color="none")+
  facet_wrap(vars(pellet))
pie
save_plot("ambient temperature correlation types by pellet",w=4,h=2)
pie+data%>%filter(pellet=="pre-OVX")+theme(strip.text = element_blank())
save_plot("ambient tempertaure correlation types intact",w=1.5,h=1.5)
pie+data%>%filter(pellet!="pre-OVX")%>%mutate(pellet=factor(pellet,levels=c("OVX+Veh","OVX+E2"),labels=c("OVX+Vehicle","OVX+E2")))+labs(title="Treatment ns")
save_plot("ambient temperature correlation types by pellet ovx",w=3,h=1.8)

##correlation coefficient
t_test(unit_df%>%filter(ambient_temp_interpolated_cor_sig_ambient!="neutral", gonad=="ovx")%>%group_by(mouse,pellet,ambient_temp_interpolated_cor_sig_ambient)%>%summarize(mean_coef=mean(ambient_temp_interpolated_cor_ambient))%>%group_by(ambient_temp_interpolated_cor_sig_ambient), mean_coef ~ pellet)
for (cell_type in unit_df%>%filter(!is.na(ambient_temp_interpolated_cor_sig_ambient))%>%pull(ambient_temp_interpolated_cor_sig_ambient)%>%unique()){
  if (cell_type=="neutral"){next}
  print(cell_type)
  anov<-anova(lme(data=unit_df%>%filter(gonad=="ovx",ambient_temp_interpolated_cor_sig_ambient==cell_type), fixed=ambient_temp_interpolated_cor_ambient ~ pellet, random=~1|mouse))
  print(anov)
}

p<-ggplot(unit_df%>%filter(ambient_temp_interpolated_cor_sig_ambient!="neutral"), aes(x=ambient_temp_interpolated_cor_sig_ambient, y=abs(ambient_temp_interpolated_cor_ambient),fill=ambient_temp_interpolated_cor_sig_ambient))+
  geom_violin()+
  point_cell()+
  point_mouse()+
  coord_cartesian(ylim=c(0,1))+
  labs(x=element_blank(),y="|r|")+
  scale_fill_manual(values=cell_type_scale[2:3])+
  scale_x_discrete(breaks= ~., labels=tools::toTitleCase)+
  scale_y_continuous(breaks=c(0,0.5,1))+
  ms+
  theme(legend.position = "none",
        plot.title=element_text(size=12,hjust=0,margin=margin(b=3,unit="pt")))
p+facet_wrap(vars(pellet),axes="all")
save_plot("ambient temperature coefficient by cell type and pellet", w=12,h=8)
p+(p$data)%>%filter(pellet=="pre-OVX")+scale_x_discrete(labels=c("Neg.","Pos."))
save_plot("ambient temeprature coefficient by cell type intact",w=1.6,h=1.2)
p+(p$data)%>%filter(pellet!="pre-OVX")%>%mutate(pellet=factor(pellet, levels=c("OVX+Veh","OVX+E2"),labels=c("OVX+Vehicle","OVX+E2")))+aes(x=pellet,fill=pellet)+facet_wrap(vars(ambient_temp_interpolated_cor_sig_ambient),labeller = as_labeller(tools::toTitleCase))+scale_fill_manual(values=post_ovx_scale)+labs(title="Treatment ns (both cell types)")+scale_x_discrete(guide=guide_axis(n.dodge=2))
save_plot("ambient temperature coefficient by cell type and pellet ovx",w=3.2,h=2)

##slope
t_test(unit_df%>%filter(ambient_temp_interpolated_cor_sig_ambient!="neutral", gonad=="ovx")%>%group_by(mouse,pellet,ambient_temp_interpolated_cor_sig_ambient)%>%summarize(mean_slope=mean(ambient_temp_interpolated_slope_ambient))%>%group_by(ambient_temp_interpolated_cor_sig_ambient), mean_slope ~ pellet)
for (cell_type in unit_df%>%filter(!is.na(ambient_temp_interpolated_cor_sig_ambient))%>%pull(ambient_temp_interpolated_cor_sig_ambient)%>%unique()){
  if (cell_type=="neutral"){next}
  print(cell_type)
  anov<-anova(lme(data=unit_df%>%filter(gonad=="ovx",ambient_temp_interpolated_cor_sig_ambient==cell_type), fixed=ambient_temp_interpolated_slope_ambient ~ pellet, random=~1|mouse))
  print(anov)
}

p<-ggplot(unit_df%>%filter(ambient_temp_interpolated_cor_sig_ambient!="neutral"), aes(x=ambient_temp_interpolated_cor_sig_ambient, y=abs(ambient_temp_interpolated_slope_ambient),fill=ambient_temp_interpolated_cor_sig_ambient))+
  geom_violin()+
  point_cell()+
  point_mouse()+
  labs(x=element_blank(),y="|Slope|")+
  scale_fill_manual(values=cell_type_scale[2:3])+
  scale_x_discrete(breaks= ~., labels=tools::toTitleCase)+
  ms+
  theme(legend.position = "none")
p+facet_wrap(vars(pellet),axes="all")
save_plot("ambient temperature slope by cell type and pellet", w=12,h=8)
p+(p$data)%>%filter(pellet=="pre-OVX")+scale_x_discrete(labels=c("Neg.","Pos."))
save_plot("ambient temeprature slope by cell type intact",w=1.6,h=1.2)
p+(p$data)%>%filter(pellet!="pre-OVX")%>%mutate(pellet=factor(pellet, levels=c("OVX+Veh","OVX+E2"),labels=c("OVX+Vehicle","OVX+E2")))+aes(x=pellet,fill=pellet)+facet_wrap(vars(ambient_temp_interpolated_cor_sig_ambient),axes="all",labeller = as_labeller(tools::toTitleCase))+scale_fill_manual(values=post_ovx_scale)+scale_x_discrete(guide=guide_axis(n.dodge=2))
save_plot("ambient temperature slope by cell type and pellet ovx",w=3.2,h=2)

##ambient LM results ----
#LM significance by pellet
data<-transform_data_piegraph(lm_df,"pellet","ambient_temp_interpolated_cor_sig_ambient_")

pie<-ggplot(data, aes(x="", y=percent, fill=ambient_temp_interpolated_cor_sig_ambient_)) +
  theme_prism()+
  theme_pie+
  geom_bar(stat="identity", width=1,color="white",position = position_stack(reverse=T)) +
  scale_fill_manual(values=cell_type_scale)+
  coord_polar("y", start=0)+
  geom_text(aes(label = paste0(round(percent,digits=0),"%"),color=ambient_temp_interpolated_cor_sig_ambient_,x=1.1),position = position_stack(vjust=0.5,reverse = T), size=5, fontface="bold")+
  scale_color_manual(values=c("grey90","black","black"))+
  guides(color="none")+
  facet_grid(vars(pellet))
pie
save_plot("ambient lm correlation significance by pellet",w=4,h=5)

##LM accuracy
ambient_lm_anova<-anova(lme(data=lm_df%>%filter(ambient_temp_interpolated_cor_sig_ambient_=="significant", pellet!="pre-OVX"),
                            fixed=ambient_temp_interpolated_mean_cor_ambient_ ~ pellet,
                            random=~1|mouse))

p<-ggplot(lm_df%>%filter(ambient_temp_interpolated_cor_sig_ambient_=="significant"), aes(x=pellet,y=ambient_temp_interpolated_mean_cor_ambient_))+
  geom_violin(aes(fill=pellet))+
  point_mouse()+
  labs(x=element_blank(),y="r")+
  scale_fill_manual(values=pellet_scale)+
  ms+
  coord_cartesian(ylim=c(0,1))+
  theme(legend.position = "none",
        plot.title=element_text(size=12,margin=margin(b=3,unit="pt"),hjust=0))
p
save_plot("ambient lm correlation coefficient by pellet",w=6,h=6)
p+(p$data)%>%filter(pellet=="pre-OVX")+scale_x_discrete(breaks=c())
save_plot("ambient lm correlation coefficient intact",w=1.4,h=2)
p+(p$data)%>%filter(pellet!="pre-OVX")+scale_x_discrete(labels=c("OVX+Vehicle", "OVX+E2"),guide=guide_axis(n.dodge=2))
save_plot("ambient lm accuracy ovx",w=2.2,h=2)
p+(p$data)%>%filter(pellet!="pre-OVX")+coord_cartesian(ylim=c(0.95,1))+labs(title="Treatment ns")+scale_fill_manual(values=post_ovx_scale)+scale_x_discrete(labels=c("OVX+Vehicle", "OVX+E2"),guide=guide_axis(n.dodge = 2))
save_plot("ambient lm accuracy ovx zoom y",w=2.2,h=2)

#intact cell types
data<-ambient_cell_type_lm_df%>%
  merge(lm_df%>%select(session_id,mouse),all.x=T)%>%
  bind_rows(lm_df%>%filter(gonad=="intact",!is.na(ambient_temp_interpolated_mean_cor_ambient_))%>%
              select(session_id,ambient_temp_interpolated_mean_cor_ambient_,ambient_temp_interpolated_cor_sig_ambient_,mouse)%>%
              mutate(ambient_temp_interpolated_cor_sig_ambient="All"))%>%
  mutate(ambient_temp_interpolated_cor_sig_ambient=factor(ambient_temp_interpolated_cor_sig_ambient,levels=c("All","neutral","negative","positive"),labels=c("All","Neutral","Negative","Positive")))

p<-ggplot(data, aes(x=ambient_temp_interpolated_cor_sig_ambient,y=ambient_temp_interpolated_mean_cor_ambient_))+
  geom_violin(aes(fill=ambient_temp_interpolated_cor_sig_ambient),scale="width",width=0.95)+
  point_summary(aes(color=mouse),position=position_jitter(width=0.15,height=0,seed=123),size=2.8,shape=21,stroke=1)+
  # point_indiv(size=1.5,position=position_jitter(width=0.25,height=0,seed=321))+
  labs(x=element_blank(),y="r")+
  scale_x_discrete(guide = guide_axis(n.dodge=2),labels=c("All","Neutral","Neg.","Pos."))+
  scale_fill_manual(values=c("black",cell_type_scale))+
  scale_y_continuous(breaks=seq(0,1,0.5),labels=seq(0,1,0.5))+
  ms+
  coord_cartesian(ylim=c(0,1))+
  theme(legend.position = "none",
        axis.title.y = element_text(margin = margin(r=3, unit="pt")))
p
save_plot("ambient lm correlation by cell type",w=1.7,h=1.6)

##LM predictions (ambient)
for (sess in sumdf%>%filter(session_type %in% c("cold","heat"))%>%pull(session_id)%>%unique()){
  data<-(ambient_lm_ls$predict_df)%>%filter(session_id==sess)
  data_cell_type<-ambient_cell_type_predict_df%>%
    filter(session_id==sess)%>%
    mutate(fold=as.numeric(fold))%>%
    bind_rows(data%>%mutate(ambient_temp_interpolated_cor_sig_ambient="All"))%>%
    mutate(ambient_temp_interpolated_cor_sig_ambient=factor(ambient_temp_interpolated_cor_sig_ambient,levels=c("All","neutral","negative","positive"),labels=c("All","Neutral","Negative","Positive")))
  
  p<-ggplot(data, aes(x=predicted,y=true))+
    labs(x="Predicted", y="Observed",title="T-Amb (all cells)")+
    geom_abline(slope=1,linetype="dashed",linewidth=1)+
    regression_line(linewidth=1,alpha=0.5,color="blue")+
    xy_point2(alpha=0.2,size=1)+
    ms+theme(plot.title=element_text(size=12,margin=margin(t=0,r=0,b=3,l=-15,"pt")),
             panel.spacing = unit(3,"pt"),
             axis.title.x = element_text(margin = margin(t = 3,unit="pt")),
             axis.title.y = element_text(margin = margin(r=3, unit="pt")))
  p
  save_plot(paste("predicted temp",sess,"ambient"),w=1.7,h=1.8)
  p+data_cell_type+facet_wrap(vars(ambient_temp_interpolated_cor_sig_ambient),axes="all",nrow=1)+labs(title=element_blank())
  save_plot(paste("predicted temp",sess,"ambient cell types"),w=4.7,h=1.8)
}

### dF/F0 - male social stimulus relationship ----
##plot all cells on one graph summarized by timepoint bin
data<-sumdf%>%
  filter(!is.na(df_f0_bin), session_type=="male_interaction")%>%
  merge(unit_df%>%select(male_interaction_auc, male_interaction_auc_sig, male_interaction_fc,unit_id_id, unit_id),all.x=T)%>%
  group_by(unit_id_id)%>%
  mutate(session_time_minutes = as.duration(telem_ts-ymd_hms(paste(start_date,start_time)))%>%as.numeric()/60)
male_entry_minutes<-data%>%filter(male_interaction==1)%>%group_by(session_id)%>%summarize(male_entry_minutes=min(session_time_minutes))
data<-data%>%merge(male_entry_minutes,all.x=T)%>%
  mutate(time_bin = case_when(male_interaction==0 & session_time_minutes<10 ~ "Pre-male",
                              male_interaction==0 & session_time_minutes>10 ~ "Post-male",
                              male_interaction==1 & session_time_minutes < male_entry_minutes+1.1 ~ "0-1 min.",
                              male_interaction==1 & session_time_minutes < male_entry_minutes+6 ~ "1-5 min.",
                              male_interaction==1 & session_time_minutes ~ "5-10 min.")%>%
           factor(levels=c("Pre-male","0-1 min.","1-5 min.","5-10 min.","Post-male")),
         unit_id_id = factor(unit_id_id, levels=data%>%arrange(desc(pellet),male_interaction_auc)%>%pull(unit_id_id)%>%unique()))

set<-list(theme(text=element_text(size=12),
                plot.title = element_text(size=12,margin=margin(t=3,b=3,l=0,r=0,unit="pt")),
                plot.margin = margin(t=0, b=0, l=0, r=0, "inches")))
rect_label_set<-list(geom_tile(aes(fill=male_interaction_auc_sig)),
                     scale_x_discrete(expand=c(0,0)),
                     scale_fill_manual(values=c(cell_type_scale[2],cell_type_scale[1],cell_type_scale[3])),
                     theme(axis.line=element_blank(),axis.text = element_blank(),legend.position = "none",axis.ticks = element_blank(),axis.title = element_blank()))
rect_label_pellet<-list(geom_tile(aes(fill=pellet)),
                        scale_fill_manual(values=post_ovx_scale2),
                        theme(axis.line=element_blank(),axis.text = element_blank(),legend.position = "none",axis.ticks = element_blank(),axis.title = element_blank()))

labels_intact<-ggplot(data%>%filter(gonad=="intact"),aes(x="",y=unit_id_id))+ms+set+rect_label_set
labels_ovx<-ggplot(data%>%filter(gonad=="ovx"),aes(x="",y=unit_id_id))+ms+set+rect_label_set
pellet_label_ovx<-ggplot(data%>%filter(gonad=="ovx"),aes(x="",y=unit_id_id))+ms+set+rect_label_pellet

p<-ggplot(data%>%filter(gonad=="intact"), aes(x=time_bin, y=unit_id_id))+
  labs(x=element_blank(), y=element_blank())+
  geom_tile(aes(fill=scaled_YrA_bin))+
  scale_y_discrete(breaks=c(),expand=c(0,0))+
  scale_x_discrete(expand=c(0,0))+
  scale_fill_continuous(type = "viridis", breaks = c(0, 1), labels = c("Min", "Max"),name="")+
  ms+
  theme(panel.background = element_rect(fill="black"),
        # legend.title = element_text(),
        # legend.text = element_text(size=12,face="bold"),
        axis.title.y=element_text(margin=margin(r=3,unit="pt")),
        axis.title.x = element_text(size=12,margin=margin(t=4,l=-25,unit="pt")),
        axis.text.x = element_text(angle = 270, hjust = 0.5, vjust = 0),
        legend.text = element_blank(),
        legend.position = "top",
        legend.key.width = unit(0.2,"inches"),
        legend.justification = 0,
        legend.box.spacing = unit(6,"pt"),
        axis.line.y = element_blank())
p+labels_intact+plot_layout(widths = c(15,1))
save_plot("df_f0 by male interactin time bin all intact cells",w=1.8,h=3.4)
p+data%>%filter(gonad=="ovx")+labels_ovx+pellet_label_ovx+plot_layout(widths = c(20,1,1))
save_plot("df_f0 by male interaction time bin all ovx cells",w=1.8,h=4.5)


##plot an example session
set<-list(theme(text=element_text(size=12),
                plot.title = element_text(size=12,margin=margin(t=3,b=3,l=0,r=0,unit="pt")),
                plot.margin = margin(t=0, b=0, l=-3, r=0, "pt"),
                panel.spacing = unit(0.0375,"inches"),
                axis.title.y=element_text(margin=margin(r=-35,unit="pt"))))

#all data
d<-read_rds("./output/intermediate/2_250417_circulating_E2_torpor_miniscope-pre-OVX_torpor-MT29-2025_05_23-session1-concatenated.rds")%>%
  filter(session_type=="male_interaction", !is.na(scaled_YrA))%>%
  ungroup()%>%
  mutate(session_time_minutes = session_time_minutes-min(session_time_minutes))%>%
  merge(unit_df%>%select(male_interaction_auc,male_interaction_auc_sig, male_interaction_fc, unit_id_id,session_id,unit_id), all.x=T)
data<-d%>%mutate(unit_id = factor(unit_id, levels=d%>%arrange(male_interaction_auc)%>%pull(unit_id)%>%unique()))

rect_label_set<-list(geom_tile(aes(fill=male_interaction_auc_sig)),
                     scale_x_discrete(expand=c(0,0)),
                     scale_fill_manual(values=c(cell_type_scale[2],cell_type_scale[1],cell_type_scale[3])),
                     theme(axis.line=element_blank(),axis.text = element_blank(),legend.position = "none",axis.ticks = element_blank(),axis.title = element_blank()))

p1<-ggplot(data, aes(x=session_time_minutes, y=unit_id))+ms+ls+ridge_set+set+scale_y_discrete(breaks=c())+geom_vline(xintercept=min(data%>%filter(male_interaction==1)%>%pull(session_time_minutes)),alpha=0.5)
p2<-ggplot(data, aes(x=session_time_minutes, y=male_interaction))+ms+ls+male_interaction_set+set+scale_y_continuous(expand = c(0.2,0.2),breaks=c(0,1),labels=c("-","+"))
p3<-ggplot(data, aes(x=session_time_minutes, y=temp))+ms+ls+temp_set+set+scale_y_continuous(limits=c(35.9,38.1),breaks=c(36,38))+labs(title="T-Core (Deg. C)",y=element_blank())+scale_x_continuous(breaks=c(0,4,8),expand=c(0,0))
p4<-ggplot(data%>%ungroup()%>%distinct(unit_id,.keep_all = T),aes(x="",y=unit_id))+ms+set+rect_label_set+labs(x=element_blank(),y=element_blank())
save_plot("example session all data male_interaction",plot=p1+p4+p2+plot_spacer()+p3+plot_layout(heights=c(12,1,1),widths=c(60,1),ncol=2),w=13,h=8)

#subset
data<-data%>%filter(session_time_minutes<10)

p1<-p1+data+theme(axis.line=element_blank(),legend.position = "none")+labs(title="Cells (F / max F)")
p2<-p2+data
p3<-p3+data+scale_y_continuous(limits=c(34.9,38.1),breaks=c(35,38))
p4<-p4+data%>%ungroup()%>%distinct(unit_id,.keep_all = T)
save_plot("example session subset data male_interaction",plot=p1+p4+p2+plot_spacer()+p3+plot_layout(heights=c(10,1,1),widths=c(20,1),ncol=2),w=2.6,h=3.3)

##plot values by male_interaction for each session
for (id in sumdf%>%filter(session_type=="male_interaction")%>%pull(session_id)%>%unique()){
  print(id)
  data<-sumdf%>%filter(!is.na(male_interaction), session_id==id)%>%mutate(male_interaction=factor(male_interaction))
  p<-ggplot(data, aes(x=male_interaction,y=df_f0_bin))+
    geom_violin(aes(fill=male_interaction))+
    point_indiv(alpha=0.25,size=2, position=position_jitter(width=0.25,height=0,seed=123))+
    scale_fill_manual(values=male_interaction_scale)+
    facet_wrap(vars(unit_id_id),axes="all",scales="free_y")+
    ms
  p
  save_plot(paste("df_f0 by male interaction",id),w=20,h=15)
}


p<-ggplot(unit_df%>%mutate(male_interaction_auc_sig=factor(male_interaction_auc_sig,levels=c("neutral","activated","suppressed"))), aes(x=male_interaction_auc,y=abs(log2(male_interaction_fc))))+
  xy_point3(aes(color=male_interaction_auc_sig))+
  scale_color_manual(values=cell_type_scale2)+
  ms
p
save_plot("male interaction volcano plot",w=5,h=3)

##Cell type frequencies
data<-transform_data_piegraph(unit_df, animal_var="pellet",cell_var="male_interaction_auc_sig")%>%
  mutate(male_interaction_auc_sig = factor(male_interaction_auc_sig, levels=c("neutral","activated","suppressed"),labels=c("Neutral","Activated","Suppressed")))

chisq<-data%>%
  filter(pellet!="pre-OVX")%>%
  select(-percent)%>%
  pivot_wider(names_from = male_interaction_auc_sig, values_from = n)%>%
  ungroup()%>%
  select(-pellet)%>%
  mutate(across(everything(), ~replace_na(.x,0)))%>%
  as.matrix()%>%
  chisq.test()

pie<-ggplot(data, aes(x="", y=percent, fill=male_interaction_auc_sig))+ms+theme_pie+
  theme(legend.position = "none",
        legend.title = element_text(hjust=0.5,margin=margin(t=0,b=3,l=0,r=0)),
        legend.title.position = "top",
        legend.text = element_text(size=12,face="bold"),
        plot.title=element_text(size=12,hjust=0,margin=margin(b=3,unit="pt")))+
  geom_bar(stat="identity", width=1,color="white",position = position_stack(reverse=T)) +
  scale_fill_manual(values=cell_type_scale, labels=tools::toTitleCase,name="Male social stimulus response")+
  coord_polar("y", start=0)+
  geom_text(aes(label = paste0(round(percent,digits=0),"%"),color=male_interaction_auc_sig,x=1.1),position = position_stack(vjust=0.5,reverse = T), size=5, fontface="bold")+
  scale_color_manual(values=c("black","black","black"))+
  guides(color="none")
pie+facet_wrap(vars(pellet))
save_plot("male response cell types by pellet",w=4,h=2)
pie+(pie$data)%>%filter(pellet=="pre-OVX")
save_plot("male response cell types intact",w=1.5,h=1.5)
pie+(pie$data)%>%filter(pellet!="pre-OVX")+facet_wrap(vars(pellet))+labs(title="Treatment p=0.09")
save_plot("male response cell types ovx by pellet",w=3,h=1.8)
(pie+theme(legend.position = "top"))%>%get_legend()%>%as_ggplot()
save_plot("male interaction cell type legend",w=4,h=0.5)

##AUROC
data<-unit_df%>%filter(male_interaction_auc_sig!="neutral")%>%
  mutate(male_interaction_auc_sig=factor(male_interaction_auc_sig,levels=c("activated","suppressed"),labels=c("Activated","Suppressed")))

for (cell_type in data%>%pull(male_interaction_auc_sig)%>%unique()){
  print(cell_type)
  anova<-anova(lme(data=data%>%filter(male_interaction_auc_sig==cell_type, !is.na(male_interaction_auc), pellet !="pre-OVX"), fixed = male_interaction_auc ~ pellet, random=~1|mouse))
  print(anova)
}

p<-ggplot(data, aes(x=pellet, y=2*abs(male_interaction_auc-0.5),fill=pellet))+
  geom_violin()+
  point_mouse()+
  point_cell()+
  coord_cartesian(ylim=c(0,1))+
  labs(x=element_blank(),y="ROC score (0-1)")+ #y="2 * |AUROC - 0.5|"
  scale_y_continuous(breaks=c(0,0.5,1))+
  scale_fill_manual(values=pellet_scale)+
  ms+
  theme(legend.position = "none",
        plot.title=element_text(size=12,hjust=0,margin=margin(b=3,unit="pt")))
p+facet_wrap(vars(male_interaction_auc_sig),axes="all",scales="free")
save_plot("male auc by cell type and pellet", w=12,h=8)
p+(p$data)%>%filter(pellet=="pre-OVX")+
  aes(x=male_interaction_auc_sig,fill=male_interaction_auc_sig)+
  scale_fill_manual(values=cell_type_scale[2:3])+
  scale_x_discrete(labels=c("Act.","Supp."))
save_plot("male auc intact by cell type",w=1.6,h=1.6)
p+(p$data)%>%filter(pellet!="pre-OVX")+facet_wrap(vars(male_interaction_auc_sig))+scale_fill_manual(values=post_ovx_scale)+labs(title="Treatment ns (both cell types)")+scale_x_discrete(labels=c("OVX+Vehicle","OVX+E2"),guide=guide_axis(n.dodge = 2))
save_plot("male auc by cell type and pellet ovx",w=3.2,h=2)

##fold-change
for (cell_type in data%>%pull(male_interaction_auc_sig)%>%unique()){
  print(cell_type)
  anova<-anova(lme(data=data%>%filter(male_interaction_auc_sig==cell_type, !is.na(male_interaction_auc), pellet !="pre-OVX"), fixed = male_interaction_fc ~ pellet, random=~1|mouse))
  print(anova)
}

p<-ggplot(data, aes(x=pellet, y=abs(log2(male_interaction_fc)), fill=pellet))+
  geom_violin()+
  point_mouse()+
  point_cell()+
  labs(x=element_blank(),y=expression(bold("|"*log[2]*"FC|")))+
  scale_fill_manual(values=pellet_scale)+
  ms+
  theme(legend.position = "none",
        plot.title = element_text(size=12,hjust=0,margin=margin(b=3,unit="pt")))
p+facet_wrap(vars(male_interaction_auc_sig),axes="all",scales="free")
save_plot("male log2fc auc by cell type and pellet", w=12,h=8)
p+(p$data)%>%filter(pellet=="pre-OVX")+
  aes(x=male_interaction_auc_sig, fill=male_interaction_auc_sig)+
  scale_fill_manual(values=cell_type_scale[2:3])+
  scale_x_discrete(labels=c("Act.","Supp."))
save_plot("male log2fc by cell type intact",w=1.6,h=1.6)
p+(p$data)%>%filter(pellet!="pre-OVX")+facet_wrap(vars(male_interaction_auc_sig),axes="all")+labs(title="Treatment ns (both cell types)")+scale_fill_manual(values=post_ovx_scale)+scale_x_discrete(labels=c("OVX+Vehicle","OVX+E2"), guide=guide_axis(n.dodge = 2))
save_plot("male log2fc by cell type and pellet ovx",w=3.2,h=2)

##Logistic regression
data<-event_df%>%filter(!is.na(df_f0),event=="male_added",!is.na(male_interaction))

auc_df<-tibble()
for (cell_type in c(unit_df%>%filter(!is.na(male_interaction_auc_sig))%>%pull(male_interaction_auc_sig)%>%unique(),"all")){
  cells<-unit_df%>%filter(male_interaction_auc_sig==cell_type)%>%pull(unit_id_id)%>%unique()
  if(cell_type=="all"){cells<-unique(unit_df$unit_id_id)}
  print(paste(cell_type,length(cells)))
  for (sid in unique(data$session_id)){
    print(sid)
    d<-data%>%filter(session_id==sid,unit_id_id %in% cells)%>%
      select(unit_id_id, frame, male_interaction, df_f0)%>%
      pivot_wider(values_from = df_f0, names_from = unit_id_id)%>%
      select(-frame)
    if(nrow(d)==0){
      print("No data, skipping")
      next}
    
    pcs<-partition_data(d, id_col = "frame", response="male_interaction", partition_type="standard", cv_folds=5, verbose=F)
    
    auc_values<-c()
    roc_df<-tibble()
    for (fold in 1:5){
      idx= pcs[["1"]][[fold]]
      train = d[idx$train,]
      test = d[idx$test,]
      
      model<-glm(male_interaction ~ ., data = train, family="binomial")
      predict<-predict(model, newdata=test, type="response")
      auc<-suppressMessages(roc(test%>%pull(male_interaction),predict)%>%auc())
      auc_values<-c(auc_values,auc)
      if (fold==1 & sid=="MT29_2025_05_23_session1"){
        p<-ggroc(roc(test%>%pull(male_interaction),predict), linewidth=1)+ms+
          scale_x_continuous(breaks=c(1,0.5,0),transform = "reverse")+
          scale_y_continuous(breaks=c(0,0.5,1))+
          labs(x="Specificity",y="Sensitivity",title="All cells")+
          geom_abline(slope=1,intercept=1,color="grey50",alpha=0.8,linewidth=0.6,linetype="dashed")+
          theme(plot.title = element_text(size=12,margin=margin(b=3,unit="pt")))
        save_plot(paste(sid,cell_type,"male interaction roc"),plot=p,w=1.6,h=1.4)
      }
    }
    auc_df<-rbind(auc_df, tibble("session_id"=sid, "mouse"=strsplit(sid,"_")[[1]][1], "male_interaction_auc_sig"=cell_type, mean_auc=mean(auc_values)))
  }
}
auc_df<-auc_df%>%
  mutate(male_interaction_auc_sig=factor(male_interaction_auc_sig,levels=c("all","neutral","activated","suppressed"),labels=c("All","Neutral","Activated","Suppressed")))%>%
  merge(lm_df%>%select(pellet,session_id,mouse),all.x=T)

t_test(auc_df%>%filter(pellet!="pre-OVX",male_interaction_auc_sig=="All")%>%mutate(pellet=as.character(pellet)), mean_auc~pellet)

p<-ggplot(auc_df, aes(x=male_interaction_auc_sig,y=mean_auc,fill=male_interaction_auc_sig))+ms+
  geom_violin(scale = "width",width=0.9)+
  point_mouse()+
  coord_cartesian(ylim=c(0,1))+
  scale_fill_manual(values=c("black",cell_type_scale))+
  labs(x=element_blank(),y="AUROC")+
  scale_x_discrete(labels=c("All","Neutral","Act.","Supp."),guide=guide_axis(n.dodge=2))+
  theme(legend.position = "none",
        plot.title=element_text(size=12,hjust=0,margin=margin(b=3,unit="pt")))
p+auc_df%>%filter(pellet=="pre-OVX")
save_plot("male interaction auc by cell type intact",w=2.4,h=1.4)
p+auc_df%>%filter(pellet!="pre-OVX")+facet_wrap(vars(pellet),axes="all")
save_plot("male interaction auc ovx by cell type and pellet",w=4.5,h=2.5)
p+auc_df%>%filter(pellet!="pre-OVX",male_interaction_auc_sig =="All")+aes(x=pellet,fill=pellet)+scale_fill_manual(values=post_ovx_scale)+scale_x_discrete(labels=c("OVX+Vehicle","OVX+E2"),guide = guide_axis(n.dodge = 2))
save_plot("male interaction auc ovx all cells by pellet",w=2.2,h=2)
save_plot("male interaction auc ovx all cells by pellet zoom y",plot=last_plot()+coord_cartesian(ylim=c(0.98,1))+labs(title="Treatment ns"),w=2.2,h=2)

####Comparison of tuning across stimuli ----
# Heat map
#Values
t=1
for (target in target_cols){
  heatmap_df_all<-unit_df%>%
    mutate(male_interaction_auc = 2 * (male_interaction_auc-0.5))%>%
    filter(!if_any(target_cols,is.na))%>%
    pivot_longer(cols=target_cols,names_to = "var",values_to = "val")%>%
    mutate(unit_id_id=factor(unit_id_id, levels=unit_df%>%arrange(!!sym(target_cols[t]))%>%pull(unit_id_id)%>%unique()),
           var=factor(var,levels=target_cols,labels=labs))
  
  heatmap_df_intact<-unit_df%>%
    filter(pellet=="pre-OVX")%>%
    mutate(male_interaction_auc = 2 * (male_interaction_auc-0.5))%>%
    filter(!if_any(target_cols,is.na))%>%
    pivot_longer(cols=target_cols,names_to = "var",values_to = "val")%>%
    mutate(unit_id_id=factor(unit_id_id, levels=unit_df%>%arrange(!!sym(target_cols[t]))%>%pull(unit_id_id)%>%unique()),
           var=factor(var,levels=target_cols,labels=labs))
  
  heatmap_df_ovx<-unit_df%>%
    filter(pellet!="pre-OVX")%>%
    mutate(male_interaction_auc = 2 * (male_interaction_auc-0.5))%>%
    filter(!if_any(target_cols,is.na))%>%
    pivot_longer(cols=target_cols,names_to = "var",values_to = "val")%>%
    mutate(unit_id_id=factor(unit_id_id, levels=unit_df%>%arrange(!!sym(target_cols[t]))%>%pull(unit_id_id)%>%unique()),
           var=factor(var,levels=target_cols,labels=labs))
  
  p<-ggplot(heatmap_df_intact, aes(x=var,y=unit_id_id,fill=val))+
    geom_tile()+
    labs(x=element_blank(),y="Cell ID")+
    scale_y_discrete(breaks=NULL)+
    scale_x_discrete(expand=c(0,0))+
    scale_fill_gradient2(high="red",low="blue",mid = "white",midpoint=0,name="Correlation / response")+
    scale_color_manual(values=pellet_scale)+
    ms+
    theme(axis.line.y = element_blank(),
          legend.text = element_text(size=12,face="bold"),
          legend.position = "top",
          legend.title.position = "top",
          legend.margin = margin(t=0,b=0,l=-5,r=0,unit="pt"),
          legend.title = element_text(hjust=0.5))
  p
  save_plot(paste("Tuning heatmap intact",target),w=2.8,h=3.6)
  p+heatmap_df_ovx%>%mutate(pellet=factor(pellet, levels=c("OVX+Veh","OVX+E2"),labels=c("OVX+Vehicle","OVX+E2")))+facet_wrap(vars(pellet),scales="free")+scale_x_discrete(expand=c(0,0),guide=guide_axis(n.dodge = 2))+theme(legend.position = "right",legend.title=element_blank())
  save_plot(paste("Tuning heatmap ovx",target),w=4.8,h=3.6)
  p+heatmap_df_all
  save_plot(paste("Tuning heatmap all",target),w=2.8,h=3.6)
  
  t=t+1
}

#scatter plots
combos<-combinations(length(target_cols),2,target_cols)
combos_binary<-combinations(length(target_cols_binary),2,target_cols_binary)
for (i in 1:length(combos[,1])){
  #Set up data
  data<-unit_df%>%
    mutate(male_interaction_auc = 2 * (male_interaction_auc-0.5))%>%
    filter(!is.na(!!sym(combos[i,1])))%>%
    filter(!is.na(!!sym(combos[i,2])))%>%
    filter(!!sym(combos_binary[i,1]) != "neutral" & !!sym(combos_binary[i,2]) != "neutral")%>%
    mutate(r = .data[[combos[i,2]]], pred =  .data[[combos[i,1]]])
  x_lab=labs[which(target_cols==combos[i,2])]
  y_lab=labs[which(target_cols==combos[i,1])]
  print(paste(combos[i,2],combos[i,1]))
  print(paste(x_lab,y_lab))
  
  #Calculate correlation and test for significance
  d<-data%>%filter(pellet=="pre-OVX")
  cor<-cor(d$r, d$pred, method="pearson")
  shuf_cor<-c()
  for (shuf in 1:shuffle_iterations){
    d$r_shuf <- sample(d$r, length(d$r),replace = F)
    shuf_cor<-c(shuf_cor,
                cor(d$r_shuf, d$pred, method="pearson"))
  }
  rank<-(c(cor,shuf_cor)%>%rank())[1]
  cor_sig<-case_when(rank<shuffle_iterations*0.00005 ~ "****",
                     rank>shuffle_iterations*0.99995 ~ "****",
                     rank<shuffle_iterations*0.0005 ~ "***",
                     rank>shuffle_iterations*0.9995 ~ "***",
                     rank<shuffle_iterations*0.005 ~ "**",
                     rank>shuffle_iterations*0.995 ~ "**",
                     rank<shuffle_iterations*0.025 ~ "*",
                     rank>shuffle_iterations*0.975 ~ "*",
                     T ~ " ns")
  print(paste(cor, cor_sig))
  
  #Plot
  p<-ggplot(data,aes(x=!!sym(combos[i,2]),y=!!sym(combos[i,1])))+
    regression_line()+
    xy_point2(size=2,alpha=0.3)+
    labs(x=x_lab, y=y_lab,title=paste0("Correlation",cor_sig))+
    scale_y_continuous(expand = c(0.1,0.1),breaks=c(-1,0,1))+
    scale_x_continuous(expand = c(0.1,0.1),breaks=c(-1,0,1))+
    ms+
    theme(plot.title=element_text(size=12,face="bold",hjust=0,margin=margin(t=0,b=3,l=-40,r=0)),
          axis.title = element_text(margin=margin(t=-5,r=-5,b=0,l=0,unit="pt")))
  p+(p$data)%>%filter(pellet=="pre-OVX")
  save_plot(paste("intact",combos[i,1],"-",combos[i,2]),w=1.6,h=1.6)
  p+(p$data)%>%filter(pellet!="pre-OVX")%>%mutate(pellet=factor(pellet,levels=c("OVX+Veh","OVX+E2"),labels=c("OVX+Vehicle","OVX+E2")))+facet_wrap(vars(pellet),axes="all")+labs(title=element_blank())
  save_plot(paste("ovx", combos[i,1], "-",combos[i,2]), w=4.2,h=1.8)
}

#Binary significance
heatmap_df<-unit_df%>%
  filter(!is.na(ambient_temp_interpolated_cor_ambient) & !is.na(temp_cor_torpor))%>%
  pivot_longer(cols=target_cols_binary,names_to = "var",values_to = "val")%>%
  mutate(unit_id_id=factor(unit_id_id, levels=unit_df%>%arrange(temp_cor_torpor)%>%pull(unit_id_id)%>%unique()),
         var=factor(var,levels=target_cols_binary,labels=labs),
         val=case_when(val=="positive"|val=="negative"|val=="neutral" ~ val,
                       val=="activated" ~ "positive",
                       val=="suppressed" ~ "negative"),
         val=factor(val,levels=c("neutral","positive","negative")))

p<-ggplot(heatmap_df, aes(x=var,y=unit_id_id,fill=val))+
  geom_tile()+
  scale_y_discrete(breaks=NULL)+
  # scale_x_discrete(labels = c("TCore during torpor","Ambient temp during heat/cold","Male interaction"))+
  scale_fill_manual(values=c("white","red","blue"))+
  scale_color_manual(values=pellet_scale)+
  ms
p
save_plot("Tuning heatmap binary",w=8,h=5)

#Frequencies (binary)
tune_col<-paste(labs,collapse = "_")
pie_data<-transform_data_piegraph(unit_df%>%drop_na(all_of(target_cols_binary))%>%unite(col=!!tune_col,target_cols_binary), "pellet", tune_col)
p<-ggplot(pie_data,aes(x=pellet,y=percent))+
  geom_col(aes(fill=!!sym(tune_col)))+
  ms
p

i=1
for (target in target_cols_binary){
  data<-transform_data_piegraph(unit_df%>%drop_na(all_of(target)), "pellet",target)%>%mutate(target_ = .data[[target]], target_=case_when(target_=="positive"|target_=="negative"|target_=="neutral" ~ target_,
                                                                                                                                          target_=="activated" ~ "positive",
                                                                                                                                          target_=="suppressed" ~ "negative"),
                                                                                             target_=factor(target_,levels=c("neutral","positive","negative")))
  
  pie<-ggplot(data, aes(x="", y=percent, fill=target_)) +
    theme_prism()+
    theme_pie+
    geom_bar(stat="identity", width=1,color="white",position = position_stack(reverse=T)) +
    scale_fill_manual(values=cell_type_scale)+
    coord_polar("y", start=0)+
    geom_text(aes(label = paste0(round(percent,digits=0),"%"),color=target_,x=1.1),position = position_stack(vjust=0.5,reverse = T), size=5, fontface="bold")+
    scale_color_manual(values=c("grey90","black","black"))+
    guides(color="none")+
    labs(title=labs[i])+
    facet_wrap(vars(pellet))
  pie
  save_plot(paste("correlation significance by pellet",target),w=7,h=5)
  
  i=i+1
}

## Compare lm results by variable
#Tcore only vs Tcore + change
data<-lm_df%>%select(temp_mean_cor_torpor_,temp_mean_cor_torpor_temp_change1,session_id)%>%
  pivot_longer(cols=c(temp_mean_cor_torpor_,temp_mean_cor_torpor_temp_change1), names_to = "var",values_to = "cor")%>%
  mutate(var = factor(var,levels=c("temp_mean_cor_torpor_", "temp_mean_cor_torpor_temp_change1"),labels=c("TCore","TCore + change")))

wilcox.test(cor~var, data)

p<-ggplot(data, aes(x=var, y=cor))+
  geom_violin()+
  line_pair(aes(group=session_id),size=0.5,alpha=0.6)+
  point_indiv()+
  labs(x=element_blank(),y="Pearson correlation coefficient")+
  ms
p
save_plot("Tcore vs Tcore with change correlation analysis", w=4.5,h=4.5)

data<-lm_df%>%select(temp_mean_cor_EntryArousal_torpor__entry,temp_mean_shuf_cor_EntryArousal_torpor__entry,session_id)%>%
  pivot_longer(cols=c(temp_mean_cor_EntryArousal_torpor__entry,temp_mean_shuf_cor_EntryArousal_torpor__entry), names_to = "var",values_to = "cor")%>%
  mutate(var = factor(var,levels=c("temp_mean_shuf_cor_EntryArousal_torpor__entry","temp_mean_cor_EntryArousal_torpor__entry"),labels=c("Train + test = mixed","Train = entry only, test = arousal only")),
         mouse = sub("_.*", "", session_id))%>%
  filter(!is.na(cor))%>%
  merge(lm_df%>%distinct(session_id,pellet),all.x=T)
anova(lme(data=data, fixed = cor ~ var, random=~1|mouse))

p<-p+aes(y=cor)+data+labs(y="r")
p
facet(p,"pellet")

data<-lm_df%>%select(temp_mean_cor_EntryArousal_torpor__arousal,temp_mean_shuf_cor_EntryArousal_torpor__arousal,session_id)%>%
  pivot_longer(cols=c(temp_mean_cor_EntryArousal_torpor__arousal,temp_mean_shuf_cor_EntryArousal_torpor__arousal), names_to = "var",values_to = "cor")%>%
  mutate(var = factor(var,levels=c("temp_mean_shuf_cor_EntryArousal_torpor__arousal","temp_mean_cor_EntryArousal_torpor__arousal"),labels=c("Train+test = mixed","Train = aroual, test=entry")),
         mouse = sub("_.*", "", session_id))%>%
  filter(!is.na(cor))%>%
  merge(lm_df%>%distinct(session_id,pellet),all.x=T)
anova(lme(data=data, fixed = cor ~ var, random=~1|mouse))

p<-p+data
p
facet(p,"pellet")

##Stdev of F0 signal across stimuli and pellets
data<-sumdf%>%filter(!is.na(sd_f0_bin))%>%group_by(unit_id_id, session_id_type)%>%summarize(mean_sd_f0_bin=mean(sd_f0_bin))%>%
  merge(sumdf%>%ungroup()%>%distinct(unit_id_id,  session_id_type,session_type, mouse, pellet),all.x=T)%>% #add metadata
  filter(!is.na(mean_sd_f0_bin))

anova(lme(data=data%>%filter(pellet!="pre-OVX"), fixed = mean_sd_f0_bin ~ pellet, random=~1|mouse))

for (st in unique(data$session_type)){
  p<-ggplot(data%>%filter(session_type==st), aes(x=pellet, y=mean_sd_f0_bin))+
    geom_violin(aes(fill=pellet))+
    scale_fill_manual(values=pellet_scale)+
    point_indiv(position=position_jitter(width=0.25,height=0,seed=123))+
    point_summary(aes(color=mouse),position=position_jitter(width=0.15,height=0,seed=123),alpha=0.7)+
    labs(title=st)+
    ms+theme(legend.position = "none")
  p
  save_plot(paste("F0 stdev",st),w=6,h=5)
}

##Visualize temporal lag tuning
p<-ggplot(unit_df, aes(x=temp_cor_torpor, y=strongest_temp_lag_torpor))+
  xy_point()+
  ms
p+facet_wrap(vars(pellet))
save_plot("temporal lag tuning",w=6,h=5)

##PCA
for (sid in unique(pca_time$session_id)){
  # if(verbose){print(sid)}
  #telem_ts (each dot is a timepoint, collapsed across all cells)
  # if(verbose){print("telem_ts")}
  d<-pca_time%>%filter(session_id==sid)
  
  set<-list(theme(legend.title = element_text(size=12,face="bold",vjust=0.85),
                  axis.title.x = element_text(margin=margin(t=3,unit="pt")),
                  axis.title.y = element_text(margin=margin(r=3,unit="pt")),
                  legend.justification = "center",
                  legend.key.width = unit(16,"pt"),
                  legend.position = "top",
                  legend.key.height = unit(14,"pt"),
                  legend.title.position = "left",
                  legend.box.margin = margin(t=0,r=0,b=-15,l=-15,unit="pt"),
                  legend.text = element_text(size=12,color="white",face="bold",margin=margin(t=-12,b=0,l=0,r=0,unit = "pt"))))
  
  for (sidt in unique(d$session_id_type)){
    if(grepl("torpor",sidt)){
      # if(verbose){print(sidt)}
      data<-d%>%filter(session_id_type==sidt)
      var_data<-(pca_ls$telem_ts_var)%>%filter(session_id_type==sidt)
      p<-ggplot(data, aes(x=PC1,y=PC2))+
        scale_color_viridis_c(name="T-Core",breaks=seq(25,35,5),limits=c(22,38))+
        xy_point(aes(color=temp),shape=21,size=1.8,stroke=0.9,alpha=0.5)+
        labs(x=paste0("PC1 (",var_data%>%filter(PC==1)%>%pull(var), "%)"), y=paste0("PC2 (",var_data%>%filter(PC==2)%>%pull(var),"%)"))+
        ms+set
      p
      save_plot(paste("PCA by temp",sidt),w=2,h=2)
      
      p<-ggplot(data, aes(x=temp,y=PC1))+
        xy_point2()+
        regression_line(formula=y~x)+
        labs(x="Core temperature (Deg. C)")+
        ms+set
      p
      save_plot(paste("PC1 - temp correlation",sidt),w=6,h=5)
      
      p<-ggplot(data, aes(x=PC1,y=PC2))+
        scale_color_manual(values=colors,name="Torpor status")+
        xy_point(aes(color=torpor_status))+
        ms+set
      p
      save_plot(paste("PCA by torpor status",sidt),w=6,h=5)
    }
    
    if(grepl("cold",sidt)){
      data<-d%>%filter(session_id_type==paste0(sid,"_cold") | session_id_type==paste0(sid,"_heat"))%>%filter(!is.na(ambient_temp_interpolated))
      var_data<-(pca_ls$telem_ts_var)%>%filter(session_id_type==paste0(sid,"_cold") | session_id_type==paste0(sid,"_heat"))
      # if(verbose){print(paste(sid,"ambient"))}
      p<-ggplot(data, aes(x=PC1,y=PC2))+
        scale_color_viridis_c(name="T-Amb")+
        xy_point(aes(color=ambient_temp_interpolated),shape=21,size=1.8,stroke=0.9,alpha=0.5)+
        labs(x=paste0("PC1 (",var_data%>%filter(PC==1)%>%pull(var), "%)"), y=paste0("PC2 (",var_data%>%filter(PC==2)%>%pull(var),"%)"))+
        ms+set+
        theme(legend.box.margin = margin(t=0,b=-15,l=-50,r=0,unit="pt"))
      p
      save_plot(paste("PCA by ambient_temp",sid,"ambient"),w=1.9,h=1.8)
      
      p<-ggplot(data, aes(x=ambient_temp_interpolated,y=PC1))+
        xy_point2()+
        regression_line(formula=y~x)+
        labs(x="Ambient temperature (Deg. C)")+
        ms+set
      p
      save_plot(paste("PC1 - ambient temp correlation",sid, "ambient"),w=6,h=5)
    }
    if(grepl("male_interaction",sidt)){
      data<-d%>%filter(session_id_type==sidt, !is.na(male_interaction))
      # if (verbose){print(paste(sidt))}
      p<-ggplot(data, aes(x=PC1,y=PC2))+
        scale_color_viridis_c(name="Male social stimulus")+
        xy_point(aes(color=male_interaction))+
        ms+set
      p
      save_plot(paste("PCA by male_interaction",sidt),w=6,h=5)
    }
  }
  #unit_id_id (each dot is a cell, collapsed across all timepoints)
  # if(verbose){print("unit_id_id")}
  d<-pca_cell%>%filter(session_id==sid)
  for (sidt in unique(d$session_id_type)){
    if(grepl("torpor",sidt)){
      data<-d%>%filter(session_id_type==sidt)
      p<-ggplot(data, aes(x=PC1,y=PC2))+
        xy_point(aes(color=temp_cor_torpor))+
        scale_color_viridis_c()+
        ms+set
      p
      save_plot(paste("PCA by temp_cor_torpor",sid),w=6,h=5)
    }
    if(grepl("cold",sidt)){
      data<-d%>%filter(session_id_type==paste0(sid,"_cold") | session_id_type==paste0(sid,"_heat"))
      p<-ggplot(data, aes(x=PC1,y=PC2))+
        xy_point(aes(color=ambient_temp_interpolated_cor_ambient))+
        scale_color_viridis_c()+
        ms+set
      p
      save_plot(paste("PCA by ambient_temp_cor_ambient",sid),w=6,h=5)
    }
    if(grepl("male_interaction",sidt)){
      data<-d%>%filter(session_id_type==sidt)
      p<-ggplot(data, aes(x=PC1,y=PC2))+
        xy_point(aes(color=male_interaction_auc_sig))+
        scale_color_viridis_d()+
        ms+set
      p
      save_plot(paste("PCA by male_interaction_auc_sig",sid),w=6,h=5)
    }
  }
}

## Cell locations
d<-A_all%>%merge(unit_df, all.x=T)

for (sid in unique(d$session_id)){
  data<-d%>%filter(session_id==sid)
  data<-data%>%mutate(x0=range(width)%>%mean(),
                      y0=range(height)%>%mean(),
                      dist_from_center=((width-x0)^2 + (height-y0)^2)^0.5)
  
  p<-ggplot(data, aes(x=width,y=height, color=temp_cor_sig_torpor))+
    geom_point(aes(alpha=A),shape=15,size=0.0001)+
    geom_mark_hull(aes(group=unit_id_id), expand=unit(0,"pt"),radius=unit(0,"pt"),concavity=4,color="black",linewidth=0.2)+
    scale_color_manual(values=c(cell_type_scale[2],cell_type_scale[1],cell_type_scale[3]))+
    # geom_point(x = range(data$width)%>%mean(), y = range(data$height)%>%mean(), size = max(range(data$width)%>%diff(), range(data$height)%>%diff())/4.5, shape = 21, fill = NA,color="black")+
    geom_circle(aes(x0=x0,y0=y0,r=max(dist_from_center)*1.05),color="black",linewidth=0.7,n=720)+
    labs(x=element_blank(),y=element_blank())+
    scale_x_continuous(breaks=c())+
    scale_y_continuous(breaks=c())+
    ms+
    theme(legend.position = "none",
          axis.line = element_blank()
          # panel.border = element_rect(colour = "black", fill = NA)
          # axis.title.y=element_text(margin=margin(r=3,unit = "pt")),
          # axis.title.x = element_text(margin=margin(t=3,unit="pt"))
    )
  save_plot(paste0("A by temp_cor_sig_torpor ",sid),plot=p, w=2.1,h=2)
  if (sid=="MT29_2025_05_22_session1"){
    as_ggplot(get_legend(p+theme(legend.position = "right")))
    save_plot("A by temp_cor_sig_torpor legend",w=2,h=2)
  }
}

#test for spatial localization of cell types
spatial_test_obs<-tibble()
spatial_test_null<-tibble()
pie_data<-tibble()
for (col in target_cols_binary){
  print(col)
  test<-spatial_test(d, col,shuf_iters = shuffle_iterations)
  test_obs<-(test$observed)%>%mutate(cell_type=.data[[col]],var=col)%>%select(-all_of(col))%>%merge(sumdf%>%ungroup()%>%distinct(session_id,.keep_all = T), all.x=T)
  test_null<-(test$null)%>%mutate(cell_type=.data[[col]],var=col)%>%select(-all_of(col))%>%merge(sumdf%>%ungroup()%>%distinct(session_id, .keep_all = T), all.x=T)
  spatial_test_obs<-rbind(spatial_test_obs,test_obs)
  spatial_test_null<-rbind(spatial_test_null, test_null)
  pie_data<-rbind(pie_data,
                  transform_data_piegraph(test_obs%>%merge(unit_df,all.x=T),animal_var="pellet",cell_var="sig")%>%mutate(var=col))
}

##plot observed vs null
p<-ggplot(spatial_test_obs%>%filter(session_id=="MT29_2025_05_22_session1",!is.na(obs_dist),var=="temp_cor_sig_torpor"), aes(x=cell_type,y=obs_dist))+
  geom_violin(inherit.aes = F, data=spatial_test_null%>%filter(session_id=="MT29_2025_05_22_session1",var=="temp_cor_sig_torpor",!is.na(perm_dist)), aes(x=cell_type,y=perm_dist),fill=NA)+
  point_summary(aes(color=sig),color="black")+
  coord_cartesian(ylim=c(0,NA))+
  labs(x=element_blank(), y="Distance (pixels)")+
  ms
p
save_plot("example observed vs null spatial distance",w=5,h=3)

##plot piegraph of frequencies
pie<-ggplot(pie_data%>%mutate(var=factor(var,levels=target_cols_binary, labels=c("T-Core","T-Amb","Male"))), aes(x="", y=percent, fill=sig))+ms+theme_pie+
  theme(legend.position = "right",
        legend.title = element_text(hjust=0.5,margin=margin(t=0,b=3,l=0,r=0)),
        legend.title.position = "top",
        legend.text = element_text(size=12,face="bold"))+
  geom_bar(stat="identity", width=1,color="white",position = position_stack(reverse=T)) +
  scale_fill_manual(values=c("grey50","red"), name="Spatial enrichment")+
  coord_polar("y", start=0)+
  geom_text(aes(label = paste0(round(percent,digits=0),"%"),color=sig,x=1.1),position = position_stack(vjust=0.5,reverse = T), size=5, fontface="bold")+
  scale_color_manual(values=c("black","black","black"))+
  guides(color="none")
pie+facet_grid(vars(pellet),vars(var))
save_plot("spatial analysis torpor by pellet and var",w=4,h=4)
pie+(pie$data)%>%filter(pellet=="pre-OVX")+facet_wrap(vars(var),nrow=1)+theme(legend.position = "top",legend.box.margin = margin(t=0,b=-15, l=0, r=0,unit="pt"))
save_plot("spatial analysis torpor intact",w=3,h=2)
pie+(pie$data)%>%filter(pellet!="pre-OVX")+facet_grid(vars(pellet),vars(var),switch = "y")+theme(legend.position = "top",legend.box.margin = margin(t=0,b=-15, l=0, r=0,unit="pt"))
save_plot("spatial analysis torpor ovx",w=3,h=3)

####Plot events
#example
data<-event_df%>%
  filter(event=="fed",session_id=="MT35_2026_02_10_session1",event_aligned_time_seconds%>%between(-60,60))%>%
  mutate(event_time_bin=ifelse(event_aligned_time_seconds>0,"after","before")%>%factor(levels=c("before","after")))
sum_data<-data%>%
  group_by(event_time_bin,unit_id)%>%
  summarize(mean_YrA=mean(YrA))%>%
  pivot_wider(id_cols=unit_id, names_from=event_time_bin, values_from=mean_YrA)%>%
  ungroup()%>%
  mutate(fc=after/before)
data<-data%>%mutate(unit_id = factor(unit_id, levels=sum_data%>%arrange(desc(fc))%>%pull(unit_id)%>%unique()))%>%
  merge(unit_df%>%select(unit_id_id,unit_id,temp_cor_sig_torpor),all.x=T)

set<-list(theme(text=element_text(size=12),
                plot.title = element_text(size=12,margin=margin(t=3,b=3,l=0,r=0,unit="pt")),
                plot.margin = margin(t=0, b=0, l=-3, r=0, "pt"),
                panel.spacing = unit(0.0375,"inches"),
                axis.title.y=element_text(margin=margin(r=-35,unit="pt"))))

rect_label_set<-list(geom_tile(aes(fill=temp_cor_sig_torpor)),
                     scale_x_discrete(expand=c(0,0)),
                     scale_fill_manual(values=c(cell_type_scale[2],cell_type_scale[1],cell_type_scale[3])),
                     theme(axis.line=element_blank(),axis.text = element_blank(),legend.position = "none",axis.ticks = element_blank(),axis.title = element_blank()))

p1<-ggplot(data, aes(x=event_aligned_time_seconds, y=unit_id))+ms+ls+ridge_set+set+scale_y_discrete(breaks=c())+geom_vline(xintercept=0,alpha=0.5)+labs(title="F / max(F)")+theme(axis.line.y = element_blank())
p2<-ggplot(data, aes(x=event_aligned_time_seconds, y=event_time_bin,group=mouse))+ms+ls+fed_set+set+scale_y_discrete(breaks=c("before","after"),labels=c("Fasting","Fed"))
p3<-ggplot(data, aes(x=event_aligned_time_seconds, y=temp))+ms+ls+temp_set+set+scale_y_continuous(expand = c(0.2,0.2),breaks=c(26.8,27.3))+labs(title="Core body temperature",y=element_blank(),x="Time from re-feeding (seconds)")+scale_x_continuous(breaks=seq(-60,60,20),expand=c(0,0))
p4<-ggplot(data%>%ungroup()%>%distinct(unit_id,.keep_all = T),aes(x="",y=unit_id))+ms+set+rect_label_set+labs(x=element_blank(),y=element_blank())
save_plot("example session re-feeding response",plot=p1+p4+p2+plot_spacer()+p3+plot_layout(heights=c(12,1,1),widths=c(20,1),ncol=2),w=3.6,h=6)

#compiled data
data<-event_df%>%
  filter(event_aligned_time_seconds%>%between(-60,60))%>%
  mutate(event_time_bin = cut(event_aligned_time_seconds, seq(-60,60,1), labels=seq(-60,59,1))%>%as.character()%>%as.numeric())%>%
  group_by(unit_id_id, event, event_time_bin)%>%
  summarize(mean_scaled_YrA = mean(scaled_YrA))%>%
  merge(unit_df, all.x=T)%>%
  mutate(temp_cor_sig_torpor=factor(temp_cor_sig_torpor,levels=c("neutral","negative","positive"),labels=c("Neutral","Negative","Positive")))

for (.event in unique(data$event)){
  print(.event)
  d<-data%>%filter(event==.event)
  .title<-case_when(.event=="male_removed" ~ "Male removed",
                    .event=="male_added" ~ "Male added",
                    .event=="fed" ~ "Re-feeding",
                    .event=="first_bite" ~ "First bite",
                    T ~ .event)
  
  p<-ggplot(d, aes(x=event_time_bin,y=mean_scaled_YrA,color=temp_cor_sig_torpor))+ms+
    continuous_errorbar()+
    continuous_line()+
    labs(y="F / max(F)", x="Time (seconds)",title=.title)+
    scale_color_manual(values=cell_type_scale)+
    geom_vline(xintercept=0,linewidth=1,alpha=0.5)+
    theme(plot.title = element_text(size=12,margin=margin(b=3,unit="pt")),
          legend.position = "none")
  p+facet_wrap(vars(pellet),axes="all")
  save_plot(paste(.event, "summary by pellet"),w=5,h=2.5)
  p+(p$data)%>%filter(pellet=="pre-OVX")
  save_plot(paste(.event, "summary intact"),w=3.5,h=2.5)
  p+(p$data)%>%filter(pellet!="pre-OVX")+facet_wrap(vars(pellet),axes="all")
  save_plot(paste(.event, "summary ovx"), w=4.5,h=2.5)
}

####Write final outputs
write_sessioninfo()

