### This code is written by Paul Vander as part of minian_output_analysis. It is meant to generate a spreadsheet with metadata for each recording session, which manual annotations for "session_type" can be added to
library(tidyverse)

## Things to set manually:
setwd("C:/Users/correalab/Box/correalab/Member Folders/Paul Vander/Experiments")
direcs<-c("250417_circulating_E2_torpor_miniscope/pre-OVX_torpor",
          "250417_circulating_E2_torpor_miniscope/post-OVX_torpor")
output_direc<-c("250417_circulating_E2_torpor_miniscope")
separator<-"/" #set depending on OS

## Create table
annotation_table<-data.frame()
i=1
for (dir in direcs){
  # print(paste0(".",dir))
  for (mouse in list.dirs(dir, full.names=F,recursive=F)){
    # print(paste0("..",mouse))
    for (start_date in list.dirs(paste(dir,mouse,sep=separator),full.names=F,recursive=F)){
      # print(paste0("...",start_date))
      for (session in list.dirs(paste(dir,mouse,start_date,sep=separator),full.names=F,recursive=F)){
        # print(paste0("....",session))
        for (start_time in list.dirs(paste(dir,mouse,start_date,session,sep=separator),full.names=F,recursive=F)){
          if (start_time == "concatenated"){next}
          exp<-strsplit(strsplit(dir,separator)[[1]][1],"_")[[1]][1]
          annotation_table[i,"experiment"] = exp
          annotation_table[i,"mouse"] = mouse
          annotation_table[i,"start_date"] = start_date
          annotation_table[i,"session"] = session
          annotation_table[i,"start_time"] = start_time
          annotation_table[i,"session_type"] = ""
          i=i+1
        }
      }
    }
  }
}

## Write to .csv output
write_csv(annotation_table, paste(output_direc,"session_type_metadata.csv",sep = separator))
