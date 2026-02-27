##Miniscope data
read_csv_minian <- function(file){
  d<-read_csv(file, show_col_types = F)
  d<-d%>%mutate(experiment=exp,
                mouse=mouse,
                start_date=start_date,
                session=session)
  if ("unit_id" %in% colnames(d)){
    d<-d%>%filter(unit_id %nin% bad_cells)%>% #remove cells labeled as bad in CScreener
      mutate(init_unit_id = unit_id, #preserve original unit_id labels
             unit_id = factor(unit_id)%>%as.numeric()%>%factor()) #renumbers unit_id starting from 1
    }
  return(d)
}

read_cell_label <- function(file){
  cell_label<-read_csv(file,show_col_types = F,col_names="label")
  cell_label$unit_id<-0:(nrow(cell_label)-1)
  bad_cells<-cell_label%>%filter(label==0)%>%pull(unit_id)
}

read_timestamps <- function(file){
  d<-read_csv(file, show_col_types = F)
  d<-d%>%mutate(experiment=exp,
                mouse=mouse,
                start_date=start_date,
                session=session)
  if("frame" %nin% colnames(d)){d<-d%>%rename(frame = `Frame Number`, time_ms = `Time Stamp (ms)`, buffer_index = `Buffer Index`)%>%mutate(start_time=start_time)}
  return(d)
}

scale_temporal <- function(d){
  d<-d%>%group_by(unit_id_id)%>%mutate(scaled_C = C/max(C),
                                       scaled_S = S/max(S),
                                       scaled_YrA = YrA/max(YrA))
  return(d)
}

scale_temporal_bin <- function(d){
  d<-d%>%group_by(unit_id_id)%>%mutate(scaled_C_bin = C_bin/max(C_bin),
                                       scaled_S_bin = S_bin/max(S_bin),
                                       scaled_YrA_bin = YrA_bin/max(YrA_bin))
  return(d)
}

read_motion <- function(path){
  pattern<-"motion[0-9]*.csv"
  if (!any(grepl(pattern, list.files(path)))){stop("No motion files detected")}
  files<-grep(pattern, list.files(path),value=T)%>%mixedsort()
  if (length(files) == 1){d<-read_csv_minian(paste0(path,separator,files))}
  if (length(files) > 1){
    d<-tibble()
    f<-0
    for (file in files){
      m<-read_csv_minian(paste0(path,separator,file))%>%filter(!is.na(height) & !is.na(width))
      m$frame <- (0:(nrow(m)-1)) + f
      d<-rbind(d,m)
      f<-f + nrow(m)
    }
  }
  d<-d%>%mutate(motion_distance = sqrt(height^2 + width^2))%>%rename(motion_width = width, motion_height=height)
  return(d)
}

##Single cell analysis
unit_analysis <- function(data, roc_session_type, .predictor="df_f0_bin", shuf_iters=1000){
  ##ROC analysis with binned Y variables
  print("Performing ROC analysis with binned response variables")
  roc_df<-roc_analysis(data, session_type=roc_session_type, predictor=.predictor, shuf_iters=shuf_iters)
  
  ##Correlation analysis with continuous Y variables
  print("Performing correlation analysis with continuous response variables")
  cor_df_torpor<-correlation_analysis(data, .session_type="torpor", response="temp", shuf_iters = shuf_iters)
  cor_df_ambient<-correlation_analysis(data, .session_type=c("cold","heat"), response=c("ambient_temp_interpolated","temp"), shuf_iters = shuf_iters)
  
  ##Combine data and add metadata
  d<-merge(roc_df,cor_df_torpor,all=T)%>%merge(cor_df_ambient,all=T)
  d<-merge(d,data%>%distinct(unit_id_id, .keep_all = T),all.x=T)
  return(d)
}

roc_analysis <- function(data, session_type, predictor="df_f0_bin", shuf_iters=1000){
  ##This function performs ROC analysis as in https://github.com/hongw-lab/Code_for_2024_ZhangM/blob/main/ROC.m (from https://www.nature.com/articles/s41586-023-06973-x#Sec9 paper)
  ii=1
  for (type in session_type){
    #Set up data
    auc_col<-paste0(type,"_auc")
    sig_col<-paste0(type,"_auc_sig")
    roc_df<-tibble(unit_id_id=character(0),"{auc_col}":=numeric(0), "{sig_col}":=character(0))
    if (type=="torpor"){roc_data<-data%>%filter(torpor_status=="non-torpor" | torpor_status=="deep_torpor")%>%mutate(label=ifelse(torpor_status=="non-torpor",0,1))}
    if (type=="heat"){roc_data<-data%>%filter(session_type=="heat")%>%filter(ambient_temp_bin=="22C" | ambient_temp_bin=="37C")%>%mutate(label=ifelse(ambient_temp_bin=="22C",0,1))}
    if (type=="cold"){roc_data<-data%>%filter(session_type=="cold")%>%filter(ambient_temp_bin=="22C" | ambient_temp_bin=="5C")%>%mutate(label=ifelse(ambient_temp_bin=="22C",0,1))}
    if (type=="male_interaction"){roc_data<-data%>%filter(session_type=="male_interaction", !is.na(male_interaction))%>%mutate(label=male_interaction)}
    
    #Run analysis and return results
    for (id in unique(roc_data$unit_id_id)){
      #Real data
      d<-roc_data%>%filter(unit_id_id==id)
      if (length(unique(d$label)) != 2){next} #excludes sessions without both labels (will cause error in roc function)
      roc<-d%>%roc_("label",predictor, direction="<",levels=c(0,1))
      
      #Shuffled data
      shuf_auc<-c()
      for (i in 1:shuf_iters){
        d$shuf_label<-sample(d$label,length(d$label),replace=F)
        shuf_auc<-c(shuf_auc, d%>%roc_("shuf_label",predictor, direction="<",levels=c(0,1))%>%auc())
      }
      #Determine rank and add to roc_df
      rank<-(c(auc(roc),shuf_auc)%>%rank())[1]
      auc_sig<-case_when(rank<shuf_iters*0.025 ~ "suppressed",
                         rank>shuf_iters*0.975 ~ "activated",
                         T ~ "neutral")
      roc_df<-rbind(roc_df, tibble(unit_id_id=id, "{auc_col}":=auc(roc), "{sig_col}":=auc_sig))
    }
    if (ii==1){compiled_df<-roc_df}
    if(ii>1){compiled_df<-merge(compiled_df,roc_df,all=T)}
    ii=ii+1
  }
  return(compiled_df)
}

correlation_analysis <- function(data, response, .session_type, predictor="df_f0_bin", shuf_iters=1000, method="pearson"){
  type<-case_when(.session_type == "torpor" ~ "torpor",
                  "cold" %in% .session_type & "heat" %in% .session_type ~ "ambient")[1]
  print(paste0("Session type = ",type, ", Response = ",paste(response,collapse = ", ")))
  data<-data%>%filter(session_type %in% .session_type)
  ii=1
  for (resp in response){
    #Set up data
    cor_col<-paste0(resp,"_cor_",type)
    sig_col<-paste0(resp,"_cor_sig_",type)
    slope_col<-paste0(resp,"_slope_",type)
    cor_df<-tibble(unit_id_id=character(0),"{cor_col}":=numeric(0), "{sig_col}":=character(0), "{slope_col}":=numeric(0))
    data<-data%>%mutate(r = .data[[resp]], pred =  .data[[predictor]])
    
    #Calculate
    for (id in unique(data$unit_id_id)){
      #Calculate real correlation
      d<-data%>%filter(unit_id_id==id)
      cor<-cor(d$pred, d$r, method=method)
      slope<-coef(lm(d$r ~ d$pred))[[2]]
      
      #Calculate shuffled correlations
      shuf_cor<-c()
      for (i in 1:shuf_iters){
        d$r_shuf<-sample(d$r, length(d$r), replace=F)
        shuf_cor<-c(shuf_cor, cor(d$r_shuf, d$pred, method=method))
      }
      #Determine rank and add to cor_df
      rank<-(c(cor,shuf_cor)%>%rank())[1]
      cor_sig<-case_when(rank<shuf_iters*0.025 ~ "negative",
                         rank>shuf_iters*0.975 ~ "positive",
                         T ~ "neutral")
      cor_df<-rbind(cor_df, tibble(unit_id_id=id, "{cor_col}":=cor, "{sig_col}":=cor_sig, "{slope_col}":=slope))
    }
    if (ii==1){compiled_df<-cor_df}
    if (ii>1){compiled_df<-merge(compiled_df,cor_df,all=T)}
    ii=ii+1
  }
  return(compiled_df)
}

##Telemetry
read_telemetry_data <- function(file, metadata_file, format = "starr-lifesci", idinfo = c("mouse", "misc", "measure", "misc2"), round = F){
  #Read in telemetry data
  d<-read_csv(file, show_col_types = F)
  if (format == "starr-lifesci"){ ###For STARR Life Sciences VitalView software output
    idinfo<-c("mouse", "misc", "measure", "misc2")
    colnames(d)<-gsub(" ", "_", colnames(d))#remove spaces from column names and add underscores
    d<-d%>%
      pivot_longer(3:length(colnames(d)),names_to="ID",values_to="data")%>%
      separate(ID, sep="_",idinfo)%>% ##Ignore warning about additional pieces
      select(!starts_with("misc"))%>%
      mutate(measure=case_when(measure=="Deg."~"temp",measure=="Cnts"~"act",T~"unknown"))%>%
      mutate(telem_ts = mdy_hms(paste(date,time)))
    if (round){d<-d%>%mutate(telem_ts = telem_ts%>%round_date(telem_ts, unit = "minute"))} #round data to nearest minute (only necessary if timestamps are not aligned on the minute)
    if (any(grepl(d$measure, pattern="unknown"))){stop("Unknown measure in telemetry data")}
    d<-d%>%spread(measure,data)
    
  } else {
    stop("This function currently only supports 'starr-lifesci' as a format")
  }
  
  #Read in telemetry metadata and add to telemetry data
  md<-read_csv(metadata_file, show_col_types = F)%>%mutate(fstart = mdy_hms(paste(fstart, "10:00:00"))) #read metadata and convert fstart to date-time
  md<-md%>%mutate(group_gonad = factor(paste0(group,"_",gonad),levels=c("veh_intact","veh_ovx","e2_intact","e2_ovx")), #useful combination column
                  pellet = factor(pellet, levels=c("pre-OVX","OVX+Veh","OVX+E2"))) 
  trial1_end <- md%>%filter(trial=="trial1")%>%pull(fstart)%>%unique() + days(4) #calculate end of trial1
  d<-d%>%
    mutate(trial=ifelse(telem_ts < trial1_end, "trial1","trial2"))%>% #add trial metadata to telem data
    merge(md)%>% #add metadata
    mutate(aligned_time  = difftime(telem_ts, fstart, units = "hours")%>%as.numeric()%>%round(digits=2))
  
  #Process telemetry data
  #Calculate lagged data
  d<-d%>%
    mutate(telem_ts = telem_ts + minutes(1))%>%
    select(telem_ts,temp,mouse)%>%rename(temp_lag1 = temp)%>%
    merge(d, all.y=T)%>%
    mutate(temp_change1 = temp-temp_lag1)
  
  #Define torpor states
  d<-d%>%mutate(torpor_status = case_when(temp<31 ~ "deep_torpor",
                                          temp<34 & temp_change1 < -0.05 ~ "entry",
                                          temp<34 & temp_change1 > 0.1 ~ "arousal",
                                          temp<34 ~ "shallow_torpor",
                                          temp>=34 ~ "non-torpor"),
                torpor_status=factor(torpor_status, levels=c("non-torpor","shallow_torpor","entry","arousal","deep_torpor")))

  return(d)
}

##Ambient temeprature
read_ambient_data <- function(file, upsample_interval_seconds = 60){
  ## Read in ambient temeprature data file
  d<-read_csv(file, show_col_types=F)%>%
    mutate(ambient_ts = ymd_hms(paste(start_date,time)))
  
  ## Upsample data 
  # Get target times for upsampling and add to d
  start_stop_df<-d%>%group_by(start_date,session_type,mouse)%>%summarize(start = min(ambient_ts), stop = max(ambient_ts))
  times_df<-tibble()
  for (i in 1:nrow(start_stop_df)){
    times_df<-rbind(times_df,
                    tibble("ambient_ts" = (seq(start_stop_df[[i,"start"]], start_stop_df[[i,"stop"]], seconds(60))),
                           "mouse"=start_stop_df[[i,"mouse"]],
                           "start_date"=start_stop_df[[i,"start_date"]],
                           "session_type"=start_stop_df[[i,"session_type"]]))
  }
  d<-merge(d, times_df, all=T)%>%
    mutate(timepoint = factor(ambient_ts)%>%as.numeric())
  
  # Upsample via linear interpolation
  interp<-approx(x = d$timepoint, y=d$ambient_temp, xout=unique(d$timepoint))
  interp_df<-tibble("timepoint"=interp$x, "ambient_temp_interpolated"=(interp$y)%>%round(digits=2))
  d<-merge(d,interp_df)%>%
    select(!c(timepoint,time))%>%
    mutate(ambient_temp_bin=case_when(ambient_temp_interpolated < 5.5 ~ "5C",
                                     ambient_temp_interpolated > 36.5 ~ "37C",
                                     ambient_temp_interpolated%>%between(21,23) ~ "22C",
                                     T ~ NA),
           ambient_temp_bin=factor(ambient_temp_bin,levels=c("5C","22C","37C")))
  return(d)
}

## Event metadata
read_event_metadata <- function(file){
  ## Read in data and get path to timeStamp files required
  d<-read_csv(file, show_col_types=F)%>%
    filter(!is.na(webcam_frame))%>%
    mutate(across(everything(),as.character),
           path=paste(paste0(experiment,"_circulating_E2_torpor_miniscope"),timepoint,mouse,start_date,session,start_time,"My_WebCam","timeStamps.csv",sep=separator))
  ts_d<-tibble()
  for (pth in unique(d$path)){
    d2<-d%>%filter(path==pth)
    t<-read_csv(pth, show_col_types = F)%>%
      rename(webcam_frame = `Frame Number`, webcam_time_ms = `Time Stamp (ms)`, buffer_index = `Buffer Index`)%>%
      mutate(exp=unique(d2$experiment),
             mouse=unique(d2$mouse),
             start_date=unique(d2$start_date),
             session=unique(d2$session),
             start_time=unique(d2$start_time))
    ts_d<-rbind(ts_d, t)
  }
  d<-merge(d,ts_d,all.x=T)%>%mutate(session_id = paste(mouse,start_date,session,sep="_"))
  for (i in 1:nrow(d)){
    frm<-(df%>%
            filter(session_id==d[[i,"session_id"]], start_time==d[[i,"start_time"]])%>%
            mutate(ts_diff = abs(time_ms-d[[i,"webcam_time_ms"]]))%>%
            filter(ts_diff == min(ts_diff))%>%
            pull(frame)%>%
            unique())[1] # [1] is important in case there is a "tie" between two frames for closest to the time stamp of interest
    d[i,"frame"] <- frm
  }
  d<-d%>%filter(!is.na(frame))
  return(d)
}

##Saving plots/data
save_plot<- function(name, plot=last_plot(), direc="./output/",w=NA,h=NA,units="in", ...){
  # print("Saving pdf...")
  ggsave(filename = paste0(direc,name,".pdf"), plot = plot, width=w, height=h, units=units, ...)
  # print("Saving svg...")
  ggsave(filename = paste0(direc,name,".svg"), plot = plot, width=w, height=h, units=units, ...)
  # print("Saving png...")
  ggsave(filename = paste0(direc,name,".png"), plot = plot, width=w, height=h, units=units, ...)
  # print("Saving RDS...")
  # saveRDS(object=last_plot(), file = paste0(direc,name,".rds"))
}

save_plot_lowdpi <- function(name, plot=last_plot(), direc="./output/",w=NA,h=NA,units="in"){
  # print("Saving pdf...")
  ggsave(filename = paste0(direc,name,".pdf"), plot = plot, width=w, height=h, units=units,dpi=72)
  # print("Saving svg...")
  ggsave(filename = paste0(direc,name,".svg"), plot = plot, width=w, height=h, units=units)
  # print("Saving png...")
  ggsave(filename = paste0(direc,name,".png"), plot = plot, width=w, height=h, units=units)
  # print("Saving RDS...")
  # saveRDS(object=last_plot(), file = paste0(direc,name,".rds"))
}

save_plot_facet <- function(name,facet_name=NULL, plot=last_plot(), by="mouse", direc="./output/",w=NA,h=NA,fw=NA,fh=NA,units="in"){
  #Define facet name
  if (is.null(facet_name)){facet_name<-paste(name, "by", by)}
  
  # print("Saving pdf...")
  ggsave(filename = paste0(direc,name,".pdf"), plot = plot, width=w, height=h, units=units)
  # print("Saving svg...")
  ggsave(filename = paste0(direc,name,".svg"), plot = plot, width=w, height=h, units=units)
  # print("Saving png...")
  ggsave(filename = paste0(direc,name,".png"), plot = plot, width=w, height=h, units=units)
  # print("Saving pdf...")
  ggsave(filename = paste0(direc,facet_name,".pdf"), plot = facet(plot,by), width=fw, height=fh, units=units)
  # print("Saving svg...")
  ggsave(filename = paste0(direc,facet_name,".svg"), plot = facet(plot,by), width=fw, height=fh, units=units)
  # print("Saving png...")
  ggsave(filename = paste0(direc,facet_name,".png"), plot = facet(plot,by), width=fw, height=fh, units=units)
  # print("Saving RDS...")
  # saveRDS(object=last_plot(), file = paste0(direc,name,".rds"))
}

save_plot_facet_exp <- function(name,facet_name=NULL, plot=last_plot(), by="exp", direc="./output/",w=NA,h=NA,fw=NA,fh=NA,units="in"){
  #Define facet name
  if (is.null(facet_name)){facet_name<-paste(name, "by", by)}
  
  # print("Saving pdf...")
  ggsave(filename = paste0(direc,name,".pdf"), plot = plot, width=w, height=h, units=units)
  # print("Saving svg...")
  ggsave(filename = paste0(direc,name,".svg"), plot = plot, width=w, height=h, units=units)
  # print("Saving png...")
  ggsave(filename = paste0(direc,name,".png"), plot = plot, width=w, height=h, units=units)
  # print("Saving pdf...")
  ggsave(filename = paste0(direc,facet_name,".pdf"), plot = facet(plot,by), width=fw, height=fh, units=units)
  # print("Saving svg...")
  ggsave(filename = paste0(direc,facet_name,".svg"), plot = facet(plot,by), width=fw, height=fh, units=units)
  # print("Saving png...")
  ggsave(filename = paste0(direc,facet_name,".png"), plot = facet(plot,by), width=fw, height=fh, units=units)
  # print("Saving RDS...")
  # saveRDS(object=last_plot(), file = paste0(direc,name,".rds"))
}

save_pdf <- function(name, plot=last_plot(), direc="./output/",w=NA,h=NA,units="in"){
  # print("Saving pdf...")
  ggsave(filename = paste0(direc,name,".pdf"),plot = plot, width=w, height=h, units=units)
}

save_png_large <- function(name, plot=last_plot(), direc="./output/",w=NA,h=NA,units="in",dpi=300){ #to be used with very large plots (raster plots is what I wrote this for). Using ggsave seems to cause issues with color
  png(filename = paste0(direc,name,".png"), width=w, height=h, res=dpi, units=units)
  print(plot)
  dev.off()
}

stat_save<- function(data, name=NA, direc="./output/"){
  if ("groups" %in% colnames(data)){
    data$groups<-as.character(data$groups)}
  write_csv(x=data, file=paste0(direc,ifelse(is.na(name),deparse(substitute(data)),name),".csv"))}

write_output<-function(data, direc="./output/", name = NA){
  saveRDS(object=data, file=paste0(direc,ifelse(is.na(name),deparse(substitute(data)),name),".rds"))
  write_csv(x=data, file=paste0(direc,ifelse(is.na(name),deparse(substitute(data)),name),".csv"))}

write_output_rds<-function(data, direc="./output/", name = NA){
  saveRDS(object=data, file=paste0(direc,ifelse(is.na(name),deparse(substitute(data)),name),".rds"))}

##ggplot themes/aesthetics/settings
point_summary <- function(...,size=8){
  geom_point(size=size,stat="summary",...)
}

point_indiv <- function(...,fill="grey60",size=3,seed=123,position=position_jitter(width=0.05,height=0,seed=seed),shape=21,color="grey20",alpha=0.7){
  geom_point(fill=fill,size=size,position=position,shape=shape,color=color,alpha=alpha,...)
}

point_errorbar <- function(...,width=0.35,size=1,color="grey15"){
  geom_errorbar(stat="summary",width=width,size=size,color=color,...)
}

continuous_line <- function(...,stat="summary",size=1.5){
  geom_line(size=size,stat=stat,...)
}

continuous_errorbar<- function(...,width=0,alpha=0.2){
  geom_errorbar(stat="summary",width=width,alpha=alpha,...)
}

xy_point <- function(...,size=4){
  geom_point(size=size,...)
}

xy_point2 <- function(...,size=3,shape=21,color="grey20",alpha=0.6, stroke=1.2){
  geom_point(size=size,shape=shape,color=color,alpha=alpha, stroke=stroke,...)
}

regression_line <- function(..., size=1.5, method="lm", se=F){
  geom_smooth(size=size,method = method, se = se,...)
}
  
line_error <- function(..., width=0,alpha=0.3){
  geom_errorbar(width=width,alpha=alpha,stat="summary",...)
}

line_pair <- function(..., color="grey",seed=123,position=position_jitter(width=0.05,height=0,seed=seed),size=1,alpha=0.8){
  geom_path(color=color,position=position,size=size,alpha=alpha,...)
}

draw_pvalue <- function(..., data,label.size=6,bracket.size=1){
  stat_pvalue_manual(data=data,label.size=label.size,bracket.size=bracket.size,...)
}

annotate_pvalue <- function(..., geom="text",hjust=0,x=0,y=0,p,size=4.7,fontface="bold"){ #set p equal to cell in t_test where p-value is located!
  annotate(geom=geom,hjust=hjust,x=x,y=y,label=paste0("p=",p),size=size,fontface=fontface,...)
} 

annotate_text <- function(..., geom="text",hjust=0,x=0,y=0,label="Lorem Ipsum",size=4.7,fontface="bold"){
  annotate(geom=geom,hjust=hjust,x=x,y=y,label=label,size=size,fontface=fontface,...)
}

ridgeline <- function(..., scale=0.9,linewidth=0.35,fill=NA){
  geom_ridgeline(scale=scale,linewidth=linewidth,fill=fill, ...)
}

ridgeline_guide <- function(..., var = "unit_id",color="grey"){
  geom_hline(aes(yintercept=!!sym(var)),color=color, ...)
}

colors<-c("#000000", "#E69F00", "#56B4E9", "#009E73","#F0E442", "#0072B2", "#D55E00", "#CC79A7") #http://bconnelly.net/posts/creating_colorblind-friendly_figures/

##Misc
`%nin%` <- Negate(`%in%`)

pca_matrix<- function(df, row_names_from = 1){
  df<-as.data.frame(df)
  rownames(df)<-df[,row_names_from]
  df<-df%>%select(!row_names_from)
  return(df)
}

write_sessioninfo<- function(){
  writeLines(capture.output(sessionInfo()), "./output/sessionInfo.txt")
}
