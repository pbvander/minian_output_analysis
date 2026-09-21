### This code is written by PV. It is designed to take the final outputs from miniscope_telemetry analysis.R and package in a more user-friendly format for data-sharing following publication

library(tidyverse)
source("C:/Users/paulv/Documents/GitHub/minian_output_analysis/functions.R")

output_dir<-"C:/Users/paulv/Box/correalab/Member Folders/Paul Vander/Data/Torpor project cross-experiment analyses/Miniscope/output"
setwd(output_dir)

## Set new column names for improved interpretability (optional, removed renaming step below to keep original names)----
rename_list<-list(
  master_session_id = "session_id",
  cross_day_master_session_id = "cr_session_id",
  master_cell_id = "unit_id_id",
  cross_day_master_cell_id = "cr_unit_id_id",
  session_type = "session_type",
  frame_number = "frame",
  time_ms = "time_ms",
  timestamp = "miniscope_ts",
  treatment_group = "pellet",
  food_status = "fed_status",
  minian_yra = "YrA",
  z_scored_delta_f = "z",
  t_core = "temp",
  locomotor_activity = "act",
  t_amb = "ambient_temp_interpolated",
  male_presence = "male_interaction"
)

rename_list2<-list(
  master_session_id = "session_id",
  cross_day_master_session_id = "cr_session_id",
  master_cell_id = "unit_id_id",
  cross_day_master_cell_id = "cr_unit_id_id",
  session_type = "session_type",
  timestamp = "telem_ts",
  treatment_group = "pellet",
  minian_yra_binned = "YrA_bin",
  z_scored_delta_f_binned = "z_bin",
  t_core = "temp",
  locomotor_activity = "act",
  food_status = "fed_status",
  t_amb = "ambient_temp_interpolated",
  male_presence = "male_interaction"
)


## Read, write, organize raw frame-wise data ----
files<-grep("^2_",list.files("./int"),value = T)
if (!"./for dryad" %in% list.dirs()){
  dir.create("./for dryad")
  dir.create("./for dryad/raw framewise data")
  dir.create("./for dryad/raw framewise data/pre_ovx")
  dir.create("./for dryad/raw framewise data/post_ovx")
  for (session_type in c("torpor","heat_cold","social")){
    dir.create(paste0("./for dryad/raw framewise data/pre_ovx/",session_type))
    dir.create(paste0("./for dryad/raw framewise data/post_ovx/",session_type))
  }
  print("Folders created")
}else{print("Already exists")}

for (file in files){
  gonad=ifelse(grepl("pre-ovx", file), "pre_ovx","post_ovx")
  print(paste(gonad,file))
  d<-read_rds(paste0("./int/",file))%>%
    select(session_id,
           cr_session_id,
           unit_id_id,
           cr_unit_id_id,
           session_type,
           frame,
           time_ms,
           miniscope_ts,
           pellet,
           YrA,
           z,
           temp,
           act,
           fed_status,
           ambient_temp_interpolated,
           male_interaction
    )
  d_torpor<-d%>%filter(session_type=="torpor")%>%select(-ambient_temp_interpolated, -male_interaction)
  d_amb<-d%>%filter(session_type %in% c("heat","cold"))%>%select(-male_interaction, -fed_status)
  d_male<-d%>%filter(session_type=="male_interaction")%>%select(-ambient_temp_interpolated, -fed_status)
  
  for (dataset in list(d_torpor, d_amb, d_male)){
    if(nrow(dataset)==0){next}
    destination<-case_when(dataset$session_type[1] == "torpor" ~ "torpor",
                           dataset$session_type[1] %in% c("heat","cold") ~ "heat_cold",
                           dataset$session_type[1] == "male_interaction" ~ "social")
    mouse<-strsplit((dataset%>%pull(master_cell_id))[1], "_")[[1]][1]
    path<-paste0("./for dryad/raw framewise data/",gonad,"/",destination,"/",mouse,".csv")
    print(path)
    write_csv(dataset, path)
  }
}

## Write summarized data
sumdf<-read_rds("sumdf.rds")%>%
  select(session_id,
         cr_session_id,
         unit_id_id,
         cr_unit_id_id,
         session_type,
         telem_ts,
         pellet,
         YrA_bin,
         z_bin,
         temp,
         act,
         fed_status,
         ambient_temp_interpolated,
         male_interaction
  )
if (!dir.exists("./for dryad/binned data")){dir.create("./for dryad/binned data")}
write_csv(sumdf, "./for dryad/binned data/sumdf.csv")

A_all<-read_rds("A_all.rds")%>%select(
  session_id,
  cr_session_id,
  unit_id_id,
  cr_unit_id_id,
  cr,
  A
)
if (!dir.exists("./for dryad/spatial data")){dir.create("./for dryad/spatial data")}
write_csv(sumdf, "./for dryad/binned data/A_all.csv")
