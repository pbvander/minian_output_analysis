##Miniscope data
read_csv_minian <- function(file){
  d<-read_csv(file, show_col_types = F)
  d<-d%>%mutate(experiment=exp,
                mouse=mouse,
                start_date=start_date,
                session=session)
  if ("unit_id" %in% colnames(d)){d<-d%>%filter(unit_id %nin% bad_cells)}
  return(d)
}

read_cell_label <- function(file){
  cell_label<-read_csv(file,show_col_types = F)%>%rename(label = "1")
  cell_label$unit_id<-0:(nrow(cell_label)-1)
  bad_cells<-cell_label%>%filter(label==0)%>%pull(unit_id)
}

read_timestamps <- function(file){
  d<-read_csv(file, show_col_types = F)
  if("frame" %nin% colnames(d)){d<-d%>%rename(frame = `Frame Number`, time_ms = `Time Stamp (ms)`, buffer_index = `Buffer Index`, start_time = start_time)}
  return(d)
}

scale_temporal <- function(d){
  d<-d%>%group_by(unit_id)%>%mutate(scaled_C = C/max(C),
                                    scaled_S = S/max(S),
                                    scaled_YrA = YrA/max(YrA))
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

##Telemetry
read_telemetry_data <- function(file, format = "starr-lifesci", idinfo = c("mouse", "misc", "measure", "misc2"), round = F, round_unit = "00:01:00"){
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
    if (round){d<-d%>%mutate(telem_ts = telem_ts%>%round_date(telem_ts, unit = "5 mins"))}
    if (any(grepl(d$measure, pattern="unknown"))){stop("Unknown measure in telemetry data")}
    d<-d%>%spread(measure,data)
    
  } else {
    stop("This function currently only supports 'starr-lifesci' as a format")
  }
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
