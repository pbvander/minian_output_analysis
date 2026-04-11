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
library(pROC)
library(rstatix)
library(sperrorest)
library(ggnewscale)
library(GGally)
# summarize<-dplyr::summarize
source("C:/Users/General Correa Lab/Documents/GitHub/minian_output_analysis/functions.R")
# filter<-dplyr::filter

##### Things to set manually:
#Paths/directories:
exp_direc<-"C:/Users/General Correa Lab/Box/correalab/Member Folders/Paul Vander/Experiments"
output_dir<-"C:/Users/General Correa Lab/Box/correalab/Member Folders/Paul Vander/Data/Torpor project cross-experiment analyses/Miniscope"
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
session_id_type_to_exclude<-c("MT29_2025_05_23_session1_heat",  "MT34_2025_12_20_session1_cold", "MT35_2026_03_07_session1_cold" # These animals do not have data from the full ambient temperature challenge (some issue arose during heat and/or cold exposure that led to exclusion of one "arm" of the ambient temperature challenge). They have been excluded so that all mice represented in the analysis have the same stimuli exposure
                              )

#Parameters
shuffle_iterations<-100

#Global graph settings:
ms<-list(theme_prism())
theme_pie <- theme(axis.line=element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),axis.ticks = element_blank(),axis.text.x = element_blank(), legend.title = element_blank())
group_gonad_scale<-c("black","#56B4E9","black","#D55E00")
pellet_scale<-c("black","#56B4E9","#D55E00")
cell_type_scale<-c("black","#F0E442","#0072B2")
cell_type_scale2<-c("grey60","#F0E442","#0072B2")
ambient_temp_bin_scale<-c("#0072B2","black","#E69F00")
male_interaction_scale<-c("black","#F0E442")

##### Read in miniscope data:
### Initalize variables
bad_frames<-list()
sumdf<-tibble()
A_all<-tibble()

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
event_mdf<-event_mdf%>%select(event, session_id,start_time, event_ts)%>%pivot_wider(names_from = event,values_from = event_ts)

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

### Read in Miniscope data from each session and process individually
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
            if ("./output/int" %nin% list.dirs()){dir.create("./output/intermediate")} #create directory if it doesn't exist
            write_output_rds(df, direc="./output/intermediate/", name=paste0("df_checkpoint1+",gsub(separator,"+",path)))
            setwd(exp_direc)

            # Add metadata on session type
            df<-df%>%merge(session_mdf, all.x=T)%>%mutate(session_id_type = paste0(session_id,"_",session_type))%>%filter(session_id_type %nin% session_id_type_to_exclude)

            # Add event-based metadata
            df<-df%>%merge(event_mdf, all.x=T)
            df<-df%>%mutate("male_interaction" = case_when(session_type != "male_interaction" ~ NA,
                                                           session_type == "male_interaction" ~ ifelse(miniscope_ts > male_added & miniscope_ts < male_removed, 1, 0)),
                            "first_interaction_aligned_time_seconds" = case_when (session_type != "male_interaction" ~ NA,
                                                                                  session_type == "male_interaction" ~ as.duration(miniscope_ts - first_interaction)%>%as.numeric()),
                            "first_interact" = case_when(session_type != "male_interaction" ~ NA,
                                                         session_type == "male_interaction" ~ ifelse(first_interaction_aligned_time_seconds < 1, T, F)),
                            "injection" = case_when(session_type != "E2_injection" ~ NA,
                                                    session_type == "E2_injection" & miniscope_ts >= E2_injection ~ "E2",
                                                    session_type == "E2_injection" & miniscope_ts >= veh_injection ~ "Veh",
                                                    session_type == "E2_injection" & miniscope_ts >= sham_injection ~ "Sham",
                                                    session_type == "E2_injection" & miniscope_ts < sham_injection ~ "Baseline")%>%factor(levels=c("Baseline","Sham","Veh","E2")),
                            "fed_status" = case_when(session_type != "torpor" ~ NA,
                                                     session_type == "torpor" & miniscope_ts < fed ~ "Fasting",
                                                     session_type == "torpor" & miniscope_ts > first_bite ~ "Eating",
                                                     session_type == "torpor" & miniscope_ts >= fed & miniscope_ts <= first_bite ~ "Fed",
                                                     T ~ "Fasting")%>%factor(levels=c("Fasting","Fed","Eating"),),
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

            ##### Checkpoint 2
            setwd(output_dir)
            write_output_rds(df, direc="./output/intermediate/", name=paste0("df_checkpoint2+",gsub(separator,"+",path)))
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
              ls<-list(scale_x_continuous(expand=c(0,0), breaks=NULL),
                       theme(text=element_text(size=28),plot.title = element_text(size=28)),
                       labs(x="Time (minutes)"),
                       facet_wrap(vars(start_time), nrow=1, scales="free_x", space="free_x"),
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
                                         geom_vline(xintercept = (data%>%filter(first_interact)%>%pull(session_time_minutes))[1], size=1,color="red"),
                                         scale_y_continuous(breaks=c(0,1)))
              fed_set<-list(geom_line(linewidth=1),
                            labs(y="",x=element_blank(),title="Fed/fasted status"))
              injection_set<-list(geom_line(linewidth=1),
                                  labs(y="",title="Injection"))

              # All frames
              p1<-ggplot(data, aes(x=session_time_minutes, y=unit_id))+ms+ls+ridge_set
              p2<-ggplot(data, aes(x=session_time_minutes, y=motion_distance))+ms+ls+motion_set
              p3<-ggplot(data, aes(x=session_time_minutes, y=temp))+ms+ls+temp_set
            
              if (grepl("torpor",id)){
                aligned_tim<-interval(df$fstart[1], ymd_hms(paste0(start_date,"_00_00_01")))%>%time_length("hours")
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
          }
        }
      }
    }
  }
}
###Checkpoint 3
setwd(output_dir)
write_output(sumdf)
sumdf<-read_rds("./output/sumdf.rds")

##### Single-cell analysis
unit_df<-unit_analysis(sumdf%>%filter(!is.na(df_f0_bin)), roc_session_type = c("torpor","heat","cold","male_interaction"), shuf_iters=shuffle_iterations)

##### Population-level analysis
### Linear model
torpor_lm_ls<-lm_analysis(sumdf%>%filter(!is.na(df_f0_bin)), id_col="telem_ts", .session_type = "torpor", response="temp", predictor="df_f0_bin", cv_folds=5, shuf_iters=shuffle_iterations)
torpor_w_tempchange1_lm_ls<-lm_analysis(sumdf%>%filter(!is.na(df_f0_bin)), id_col="telem_ts", .session_type = "torpor", response="temp", predictor="df_f0_bin", additional_x_var = "temp_change1", cv_folds=5, shuf_iters=shuffle_iterations)
ambient_lm_ls<-lm_analysis(sumdf%>%filter(!is.na(df_f0_bin)), id_col="ambient_ts", .session_type = c("heat","cold"), response="ambient_temp_interpolated", predictor="df_f0_bin", cv_folds=5, shuf_iters=shuffle_iterations)

lm_df<-merge(torpor_lm_ls$lm_df, ambient_lm_ls$lm_df,all=T)%>%merge(torpor_w_tempchange1_lm_ls$lm_df,all=T)%>%merge(torpor_w_tempchange1_lm_ls$lm_add_x_var_coef_df,all=T)%>% #combine data
  merge(sumdf%>%ungroup()%>%distinct(session_id,.keep_all = T),all.x=T) #add metadata
lm_predict_df<-merge(torpor_lm_ls$predict_df, ambient_lm_ls$predict_df,all=T)%>%merge(torpor_w_tempchange1_lm_ls$predict_df,all=T)%>% #combine data
  merge(sumdf%>%ungroup()%>%distinct(session_id,.keep_all = T),all.x=T) #add metadata
unit_df<-merge(torpor_lm_ls$lm_coef_df, ambient_lm_ls$lm_coef_df,all=T)%>%merge(torpor_w_tempchange1_lm_ls$lm_coef_df,all=T)%>% #combine data
  merge(unit_df,all=T) #add coefficients from population model to unit_df

### PCA
pca_ls<-pca(sumdf%>%filter(!is.na(df_f0_bin)), predictor = "df_f0_bin", dims=c("unit_id_id","telem_ts"))

##### Checkpoint 3
#Write new
setwd(output_dir)
write_output(lm_df)
write_output(lm_predict_df)
write_output(unit_df)
write_output_rds(pca_ls)

#Read all
setwd(output_dir)
sumdf<-read_rds("./output/sumdf.rds")
lm_df<-read_rds("./output/lm_df.rds")
lm_predict_df<-read_rds("./output/lm_predict_df.rds")
unit_df<-read_rds("./output/unit_df.rds")
pca_ls<-read_rds("./output/pca_ls.rds")

########### Graph ###########
### Number of cells per group/session
counts<-sumdf%>%filter(!is.na(df_f0_bin))%>%ungroup()%>%distinct(unit_id_id,.keep_all = T)%>%group_by(session_id,pellet)%>%count()
t_test(counts%>%ungroup()%>%mutate(pellet=as.character(pellet))%>%filter(pellet!="pre-OVX"), n~pellet)

p<-ggplot(counts, aes(x=pellet,y=n))+
  point_indiv()+
  point_errorbar()+
  point_summary()+ms
p


### dF/F0 - body temperature relationship (single unit analysis)
## By temperature value
# Graph all sessions
for (id in sumdf%>%filter(session_type=="torpor")%>%pull(session_id)%>%unique()){
  print(id)
  data<-sumdf%>%filter(!is.na(df_f0_bin), session_type=="torpor", session_id==id)%>%
    mutate(unit_id = factor(unit_id, levels = unit_df%>%filter(session_id==id)%>%arrange(temp_cor_torpor)%>%pull(unit_id)))
  
  p<-ggplot(data, aes(x=temp, y=df_f0_bin))+
    xy_point2(alpha=0.5)+
    regression_line()+
    ms+
    labs(y="dF/F0", title="Cell ID", x="Core body temperature (Deg. C)")+
    theme(text = element_text(size=24))+
    facet_wrap(vars(unit_id), axes="all",scales="free_y")#+theme(strip.text.x = element_blank())
  p
  save_png_large(paste("df_f0 by body temperature",id), w=40,h=25)
}


# # Graph a single session
# data<-data%>%filter(session_id_type==single_session)%>%
#   mutate(unit_id=factor(unit_id, levels=unit_df%>%filter(session_id_type==single_session)%>%arrange(temp_cor_torpor)%>%pull(unit_id)))
# 
# p<-ggplot(data, aes(x=temp, y=df_f0_bin))+
#   xy_point2(alpha=0.5)+
#   regression_line()+
#   ms+
#   labs(y="dF/F0", title="Cell ID", x="Core body temperature (Deg. C)")+
#   theme(text = element_text(size=24))+
#   facet_wrap(vars(unit_id), axes="all",scales="free_y")#+theme(strip.text.x = element_blank())
# p
# save_png_large(paste("df_f0 by body temperature",single_session),w=25,h=15)

### Number of observations
##All data
#Numbers of rows in data (# of cells x # of timepoints)
p<-ggplot(sumdf%>%filter(session_type=="torpor"), aes(x=temp))+
  geom_histogram(aes(fill=pellet),position = position_dodge(),breaks=seq(min(sumdf$temp)%>%round(digits=0)-1, max(sumdf$temp)%>%round(digits=0)+1, 1))+
  scale_x_continuous(expand=c(0,0),breaks=seq(20,50,1))+
  scale_y_continuous(expand=c(0,0))+
  labs(x="Core body temperature",y="Observations")+
  scale_fill_manual(values=pellet_scale)+
  ms
p
save_png_large("torpor observations by temp and pellet", w=7,h=5)

#Number of timepoints
data<-sumdf%>%filter(session_type == "torpor")%>%ungroup()%>%distinct(telem_ts, mouse, session_id,.keep_all = T)

p<-ggplot(data, aes(x=temp))+
  geom_histogram(aes(fill=pellet),position = position_dodge(),breaks=seq(min(sumdf$temp)%>%round(digits=0)-1, max(sumdf$temp)%>%round(digits=0)+1, 1))+
  scale_x_continuous(expand=c(0,0),breaks=seq(20,50,1))+
  scale_y_continuous(expand=c(0,0))+
  labs(x="Core body temperature",y="Timepoints")+
  scale_fill_manual(values=pellet_scale)+
  ms
p
save_png_large("torpor timepoints by temp and pellet", w=7,h=5)

##Downsampled
data<-sumdf%>%filter(session_type == "torpor")%>%downsample_data_temporal()

p<-ggplot(data, aes(x=temp))+
  geom_histogram(aes(fill=pellet),position = position_dodge(),breaks=seq(min(sumdf$temp)%>%round(digits=0)-1, max(sumdf$temp)%>%round(digits=0)+1, 1))+
  scale_x_continuous(expand=c(0,0),breaks=seq(20,50,1))+
  scale_y_continuous(expand=c(0,0))+
  labs(x="Core body temperature",y="Observations")+
  scale_fill_manual(values=pellet_scale)+
  ms
p
save_png_large("torpor observations by temp and pellet downsampled", w=7,h=5)

#Number of timepoints
data<-data%>%ungroup()%>%distinct(telem_ts, mouse, session_id,.keep_all = T)

p<-ggplot(data, aes(x=temp))+
  geom_histogram(aes(fill=pellet),position = position_dodge(),breaks=seq(min(sumdf$temp)%>%round(digits=0)-1, max(sumdf$temp)%>%round(digits=0)+1, 1))+
  scale_x_continuous(expand=c(0,0),breaks=seq(20,50,1))+
  scale_y_continuous(expand=c(0,0))+
  labs(x="Core body temperature",y="Timepoints")+
  scale_fill_manual(values=pellet_scale)+
  ms
p
save_png_large("torpor timepoints by temp and pellet downsampled", w=7,h=5)

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

# Correlation coefficient by gonad/E2 state or pellet
p<-ggplot(unit_df, aes(x=group_gonad, y=temp_cor_torpor))+
  geom_violin(aes(fill=group_gonad))+
  point_indiv()+
  scale_fill_manual(values=group_gonad_scale)+
  facet_wrap(vars(temp_cor_sig_torpor),axes="all")+
  ms
p
save_plot("torpor temperature correlation coefficient by cell type and group_gonad", w=12,h=8)

t_test(unit_df%>%filter(temp_cor_sig_torpor!="neutral", gonad=="ovx")%>%group_by(mouse,pellet,temp_cor_sig_torpor)%>%summarize(mean_coef=mean(temp_cor_torpor))%>%group_by(temp_cor_sig_torpor), mean_coef ~ pellet)
for (cell_type in unique(unit_df$temp_cor_sig_torpor)){
  if (cell_type=="neutral"){next}
  print(cell_type)
  anov<-anova(lme(data=unit_df%>%filter(gonad=="ovx",temp_cor_sig_torpor==cell_type), fixed=temp_cor_torpor ~ pellet, random=~1|mouse))
  print(anov)
}

p<-ggplot(unit_df%>%filter(temp_cor_sig_torpor!="neutral"), aes(x=pellet, y=temp_cor_torpor))+
  geom_violin(aes(fill=pellet))+
  labs(y="Pearson correlation coefficient")+
  point_summary(aes(color=mouse),position=position_jitter(width=0.05,height=0,seed=123))+
  point_indiv()+
  scale_fill_manual(values=pellet_scale)+
  facet_wrap(vars(temp_cor_sig_torpor),axes="all")+
  ms+
  theme(legend.position = "none")
p
save_plot("torpor temperature correlation coefficient by cell type and pellet", w=12,h=8)

# Correlation type frequencies
#Grouped by pellet
data<-transform_data_piegraph(unit_df, animal_var = "pellet", cell_var = "temp_cor_sig_torpor")%>%
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

pie<-ggplot(data, aes(x="", y=percent, fill=temp_cor_sig_torpor)) +
  theme_prism()+
  theme_pie+
  geom_bar(stat="identity", width=1,color="white",position = position_stack(reverse=T)) +
  scale_fill_manual(values=cell_type_scale)+
  coord_polar("y", start=0)+
  geom_text(aes(label = paste0(round(percent,digits=0),"%"),color=temp_cor_sig_torpor,x=1.1),position = position_stack(vjust=0.5,reverse = T), size=5, fontface="bold")+
  scale_color_manual(values=c("grey90","black","black"))+
  guides(color="none")+
  facet_wrap(vars(pellet))
pie
save_png_large("torpor temperature correlation types by pellet",w=8,h=5)

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
save_png_large("torpor temperature correlation types by mouse and pellet",w=8,h=5)

#Slope by gonad/E2 state
t_test(unit_df%>%filter(temp_cor_sig_torpor!="neutral", gonad=="ovx")%>%group_by(mouse,pellet,temp_cor_sig_torpor)%>%summarize(mean_slope=mean(temp_slope_torpor))%>%group_by(temp_cor_sig_torpor), mean_slope ~ pellet)
for (cell_type in unique(unit_df$temp_cor_sig_torpor)){
  if (cell_type=="neutral"){next}
  print(cell_type)
  anov<-anova(lme(data=unit_df%>%filter(gonad=="ovx",temp_cor_sig_torpor==cell_type), fixed=temp_slope_torpor ~ pellet, random=~1|mouse))
  print(anov)
}

p<-ggplot(unit_df%>%filter(temp_cor_sig_torpor!="neutral"), aes(x=pellet, y=temp_slope_torpor))+
  geom_violin(aes(fill=pellet))+
  point_summary(aes(color=mouse),position=position_jitter(width=0.05,height=0,seed=123))+
  point_indiv()+
  labs(x=element_blank(),y="Slope")+
  scale_fill_manual(values=pellet_scale)+
  facet_wrap(vars(temp_cor_sig_torpor),axes="all",scales="free")+
  ms+
  theme(legend.position = "none")
p
save_plot("torpor temperature slope by cell type and pellet", w=12,h=8)

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
save_plot("df_f0 non-torpor vs deep torpor",w=20,h=15)

# All torpor status bins
p<-ggplot(sumdf%>%filter(!is.na(df_f0_bin)), aes(x=torpor_status,y=df_f0_bin))+ms+set
p
save_plot("df_f0 by torpor status",w=20,h=15)

# ROC analysis
#Cell type frequencies
#By pellet
data<-transform_data_piegraph(unit_df, animal_var="pellet", cell_var = "torpor_auc_sig")%>%
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
save_png_large("torpor roc types by pellet",w=8,h=5)

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
save_png_large("torpor roc types by pellet and mouse",w=8,h=5)

#Fold-change between status bins
p<-ggplot(unit_df%>%filter(!is.na(torpor_auc_sig))%>%mutate(torpor_auc_sig=factor(torpor_auc_sig,levels=c("neutral","activated","suppressed"))),aes(x=torpor_auc,y=abs(log2(torpor_fc))))+
  xy_point3(aes(color=torpor_auc_sig))+
  scale_color_manual(values=cell_type_scale2)+
  ms+
  theme(legend.position = "none")
p
save_plot("torpor status bins fc",w=10,h=7)

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

### dF-F0 - body temperature (population level analysis)
##LM significance by pellet
data<-transform_data_piegraph(lm_df,"pellet","temp_cor_sig_torpor_")

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

##LM accuracy
p<-ggplot(lm_df%>%filter(temp_cor_sig_torpor_=="significant"), aes(x=pellet,y=temp_mean_cor_torpor_))+
  geom_violin(aes(fill=pellet))+
  point_summary(aes(color=mouse),position=position_jitter(width=0.05,height=0,seed=123))+
  point_indiv()+
  labs(x=element_blank(),y="Correlation coefficient")+
  scale_fill_manual(values=pellet_scale)+
  ms+
  coord_cartesian(ylim=c(0,1))+
  theme(legend.position = "none")
p
save_plot("lm correlation coefficient by pellet",w=6,h=6)

##LM predictions
for (sess in sumdf%>%filter(session_type=="torpor")%>%pull(session_id)%>%unique()){
  data<-torpor_lm_ls$predict_df
  data<-data%>%filter(session_id==sess)
  
  p<-ggplot(data, aes(x=predicted,y=true))+
    labs(x="Predicted temperature", y="Observed temperature")+
    geom_abline(slope=1,linetype="dashed",size=2)+
    xy_point2(alpha=0.5)+
    regression_line()+
    ms
  p
  save_plot(paste("predicted temp",sess,"torpor"),w=6,h=6)
}

##PCA
pca_df<-tibble()
for (sid in unique(pca_ls$telem_ts$session_id)){ #Each dot is a timepoint, collapsed across all cells
  data<-(pca_ls$telem_ts)%>%filter(session_id==sid)
  set<-list()
  p<-ggplot(data, aes(x=PC1,y=PC2))+
    scale_color_viridis_c()+
    xy_point(aes(color=temp))+
    ms
  p
  save_plot(paste("PCA by temp",sid),w=6,h=5)
  
  p<-ggplot(data, aes(x=PC1,y=temp))+
    xy_point2()+
    regression_line(formula=y~x)+
    labs(y="Core body temperature (Deg. C)")+
    ms
  p
  save_plot(paste("PC1 - temp correlation",sid),w=6,h=5)
  
  p<-ggplot(data, aes(x=PC1,y=PC2))+
    scale_color_manual(values=colors)+
    xy_point(aes(color=torpor_status))+
    ms
  p
  save_plot(paste("PCA by torpor status",sid),w=6,h=5)
  
  cor<-cor(data$temp, data$PC1, method = "pearson")
  shuf<-c()
  for (i in 1:shuffle_iterations){
    d<-data%>%mutate(temp = sample(temp, length(data$temp)))
    shuf<-c(shuf, cor(d$temp,d$PC1, method="pearson"))
  }
  rank<-(c(cor, shuf)%>%rank())[1]
  cor_sig<-case_when(rank<shuffle_iterations*0.025 ~ "sig",
                     rank>shuffle_iterations*0.975 ~ "sig",
                     T ~ "non-sig")
  pca_d<-tibble("session_id" = sid, "temp_PC1_cor" = cor, "temp_PC1_cor_sig"=cor_sig)
  pca_df<-rbind(pca_df,pca_d)
}
pca_df<-pca_df%>%merge(sumdf%>%ungroup()%>%distinct(session_id,.keep_all = T),all.x=T)

p<-ggplot(pca_df%>%filter(!is.na(temp_PC1_cor),temp_PC1_cor_sig=="sig"), aes(x=pellet,y=abs(temp_PC1_cor)))+
  point_summary(aes(color=pellet))+
  point_errorbar(aes(group=pellet))+
  point_indiv()+
  scale_color_manual(values=pellet_scale)+
  ms+theme(legend.position = "none")
p
save_plot("PC1 - temp correlation by pellet", w=6,h=6)

for (sid in unique(pca_ls$unit_id_id$session_id)){ #Each dot is a cell, collapsed across timepoints
  data<-(pca_ls$unit_id_id)%>%filter(session_id==sid)%>%merge(unit_df,all.x=T)
  p<-ggplot(data, aes(x=PC1,y=PC2))+
    xy_point(aes(color=temp_cor_torpor))+
    scale_color_viridis_c()+
    ms
  p
  save_plot(paste("PCA by temp_cor_torpor",sid),w=6,h=5)
}

### dF/F0 - ambient temperature relationship
## Plot ambient temperature challenge schematic
ggplot(sumdf%>%filter(session_id=="MT30_2025_05_23_session1", session_type %in% c("cold","heat"))%>%mutate(session_type=factor(session_type,levels=c("heat","cold"))),aes(x=session_time_minutes,y=ambient_temp_interpolated))+
  labs(x="Time (minutes)",y="Ambient temperature (Deg. C)")+
  continuous_line()+
  ms+theme(panel.spacing=unit(1,"in"))+
  facet_wrap(vars(session_type),scales="free_x")
save_plot("ambient temperature schematic",w=7,h=5)

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
  save_png_large(paste("df_f0 by ambient temperature",id), w=40,h=25)
  
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

t_test(unit_df%>%filter(temp_cor_sig_torpor!="neutral", gonad=="ovx")%>%group_by(mouse,pellet,temp_cor_sig_torpor)%>%summarize(mean_slope=mean(temp_slope_torpor))%>%group_by(temp_cor_sig_torpor), mean_slope ~ pellet)
for (cell_type in unique(unit_df$temp_cor_sig_torpor)){
  if (cell_type=="neutral"){next}
  print(cell_type)
  anov<-anova(lme(data=unit_df%>%filter(gonad=="ovx",temp_cor_sig_torpor==cell_type), fixed=temp_slope_torpor ~ pellet, random=~1|mouse))
  print(anov)
}

p<-ggplot(unit_df%>%filter(ambient_temp_interpolated_cor_sig_ambient!="neutral"), aes(x=pellet, y=ambient_temp_interpolated_cor_ambient))+
  geom_violin(aes(fill=pellet))+
  point_summary(aes(color=mouse),position=position_jitter(width=0.05,height=0,seed=123))+
  point_indiv()+
  labs(x=element_blank(),y="Slope")+
  scale_fill_manual(values=pellet_scale)+
  facet_wrap(vars(ambient_temp_interpolated_cor_sig_ambient),axes="all",scales="free")+
  ms+
  theme(legend.position = "none")
p
save_plot("ambient temperature slope by cell type and pellet", w=12,h=8)

### dF/F0 - male social stimulus relationship
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

p<-ggplot(unit_df%>%filter(male_interaction_auc_sig!="neutral"), aes(x=pellet, y=male_interaction_auc))+
  geom_violin(aes(fill=pellet))+
  point_summary(aes(color=mouse),position=position_jitter(width=0.05,height=0,seed=123))+
  point_indiv()+
  labs(x=element_blank(),y="Slope")+
  scale_fill_manual(values=pellet_scale)+
  facet_wrap(vars(male_interaction_auc_sig),axes="all",scales="free")+
  ms+
  theme(legend.position = "none")
p
save_plot("male auc by cell type and pellet", w=12,h=8)

p<-ggplot(unit_df%>%filter(male_interaction_auc_sig!="neutral"), aes(x=pellet, y=log2(male_interaction_fc)))+
  geom_violin(aes(fill=pellet))+
  point_summary(aes(color=mouse),position=position_jitter(width=0.05,height=0,seed=123))+
  point_indiv()+
  labs(x=element_blank(),y=expression(bold(log[bold(2)](fold~change))))+
  scale_fill_manual(values=pellet_scale)+
  facet_wrap(vars(male_interaction_auc_sig),axes="all",scales="free")+
  ms+
  theme(legend.position = "none")
p
save_plot("male log2fc auc by cell type and pellet", w=12,h=8)

##Comparison of tuning across stimuli
target_cols<-c("temp_cor_torpor","temp_change1_cor_torpor", 
               "ambient_temp_interpolated_cor_ambient","male_interaction_auc")
target_cols_binary<-c("temp_cor_sig_torpor","temp_change1_cor_sig_torpor",
                      "ambient_temp_interpolated_cor_sig_ambient","male_interaction_auc_sig")
labs<-c("TCore", "TCore_change",
        "TAmb", "Male")

# Heat map
#Values
t=1
for (target in target_cols){
  heatmap_df<-unit_df%>%
    filter(!is.na(ambient_temp_interpolated_cor_ambient))%>%
    pivot_longer(cols=target_cols,names_to = "var",values_to = "val")%>%
    mutate(unit_id_id=factor(unit_id_id, levels=unit_df%>%arrange(!!sym(target_cols[t]))%>%pull(unit_id_id)%>%unique()),
           var=factor(var,levels=target_cols,labels=labs))

  p<-ggplot(heatmap_df, aes(x=var,y=unit_id_id,fill=val))+
    geom_tile(data=heatmap_df%>%filter(var!="Male"))+
    scale_y_discrete(breaks=NULL)+
    # scale_x_discrete(labels = c("TCore during torpor","Ambient temp during heat/cold","Male interaction"))+
    scale_fill_gradient2(high="red",low="blue",mid = "white",midpoint=0)+
    new_scale_fill()+
    geom_tile(data=heatmap_df%>%filter(var=="Male"),aes(fill=val))+
    scale_fill_gradient2(high="red",low="blue",mid="white",midpoint=0.5)+
    scale_color_manual(values=pellet_scale)+
    ms
  p
  save_plot(paste("Tuning heatmap",target),w=8,h=5)
  
  t=t+1
}

#scatter plots
combos<-combinations(length(target_cols),2,target_cols)
combos_binary<-combinations(length(target_cols_binary),2,target_cols_binary)
for (i in 1:length(combos[,1])){
  data<-unit_df%>%
    filter(!is.na(!!sym(combos[i,1])))%>%
    filter(!is.na(!!sym(combos[i,2])))%>%
    filter(!!sym(combos_binary[i,1]) != "neutral" & !!sym(combos_binary[i,2]) != "neutral")%>%
    mutate(r = .data[[combos[i,1]]], pred =  .data[[combos[i,2]]])
  print(combos[i,])
  print(dim(data))
  print(cor(data$r, data$pred, method="pearson"))
  p<-ggplot(unit_df%>%filter(!!sym(combos_binary[i,1]) != "neutral" & !!sym(combos_binary[i,2]) != "neutral"),aes(x=!!sym(combos[i,1]),y=!!sym(combos[i,2])))+
    regression_line()+
    xy_point2()+ms
  p
  save_plot(paste(combos[i,1],"-",combos[i,2]),w=5,h=4)
}

#Binary significance
heatmap_df<-unit_df%>%
  filter(!is.na(ambient_temp_interpolated_cor_ambient))%>%
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
  save_plot(paste("lm correlation significance by pellet",target),w=7,h=5)
  
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

##Stdev of F0 signal across stimuli and pellets
data<-sumdf%>%filter(!is.na(sd_f0_bin))%>%group_by(unit_id_id, session_id_type)%>%summarize(sd_f0=mean(sd_f0_bin))%>%
  merge(sumdf%>%ungroup()%>%distinct(unit_id_id,  session_id_type, .keep_all = T),all.x=T)%>% #add metadata
  filter(!is.na(sd_f0_bin))

for (st in unique(data$session_type)){
  p<-ggplot(data%>%filter(session_type==st), aes(x=pellet, y=sd_f0))+
    geom_violin(aes(fill=pellet))+
    scale_fill_manual(values=pellet_scale)+
    point_indiv()+
    point_summary(aes(color=mouse),position=position_jitter(width=0.05,height=0,seed=123),alpha=0.7)+
    labs(title=st)+
    ms+theme(legend.position = "none")
  p
  save_plot(paste("F0 stdev",st),w=6,h=5)
}

####Write final outputs
write_sessioninfo()

