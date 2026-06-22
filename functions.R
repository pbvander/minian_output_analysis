##Reading and processing Miniscope data
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

equalize_data_temporal <- function(data, response="temp", grouping_variable="pellet", temporal_variable = "telem_ts", subject_variable="mouse", bin_width=1,verbose=T){
  if(verbose){print("Warning! This function could undersample data to equalize if the unique number of subjects across groups exceeds 2")}
  #Set up column names and group size
  bin_col<-paste0(response,"_bin",bin_width)
  total_groups<-data%>%pull({{grouping_variable}})%>%unique()%>%length()
  subject_counts<-data%>%group_by(across(all_of(c(grouping_variable,subject_variable))))%>%count()%>%group_by(across(all_of(grouping_variable)))%>%summarize(n_subjects=n())%>%ungroup()%>%mutate(scale=n_subjects/min(n_subjects))

  #Create bin_col by using cut on response
  d<-data%>%mutate({{bin_col}}:=cut(.data[[response]], seq(min(data%>%pull({{response}}))%>%floor(), max(data%>%pull({{response}}))%>%ceiling(), bin_width)))
  
  #Get the number of timepoints per mouse in each group/bin combination
  timepoint_counts<-d%>%ungroup()%>%distinct(.data[[temporal_variable]], .data[[subject_variable]], .keep_all = T)%>%group_by(across(all_of(c(bin_col,grouping_variable))))%>%summarize(timepoints_per_pellet=n())
  
  #Calculate minimum number of timepoints present in each bin across all groups (set to 0 if not all groups are represented)
  target_counts<-timepoint_counts%>%
    group_by(across(all_of(bin_col)))%>%
    mutate(n_distinct=n_distinct(.data[[grouping_variable]]),
           target_timepoints_per_pellet=ifelse(n_distinct==total_groups, min(timepoints_per_pellet), 0))%>%
    merge(subject_counts,all=T)%>%
    mutate(target_timepoints_per_pellet_scaled = target_timepoints_per_pellet * scale)%>%
    group_by(across(all_of(bin_col)))%>%
    mutate(target_timepoints_per_pellet_scaled = case_when(any(target_timepoints_per_pellet_scaled > timepoints_per_pellet) ~ target_timepoints_per_pellet_scaled / max(scale),
                                                           T ~target_timepoints_per_pellet_scaled))
  
  #Filter ONLY temporally (collapse across other variables (cells), then filter timepoints)
  d2<-d%>%
    merge(target_counts,all=T)%>%
    distinct(.data[[temporal_variable]], .data[[subject_variable]],.keep_all = T)%>%
    group_by(across(all_of(c(bin_col,grouping_variable))))%>%
    group_modify(~ slice_sample(.x, n=first(.x$target_timepoints_per_pellet_scaled)%>%round(digits=0)))
  # d2%>%group_by(temp_bin1,pellet)%>%count()  #Check that number of timepoints in each bin are the same across groups
  
  #Apply filtering of timepoints to original dataset
  timepoint_subjects_to_keep<-d2%>%ungroup()%>%select(all_of(c(temporal_variable,subject_variable, bin_col)))
  ds_data<-merge(timepoint_subjects_to_keep,data)
}

##Single cell analysis
unit_analysis <- function(data, roc_session_type, .predictor="df_f0_bin", lag_window=seq(-5,5,1), shuf_iters=1000, verbose=T){
  ##ROC analysis with binned Y variables
  if (verbose){print("Performing ROC analysis with binned response variables")}
  roc_df<-roc_analysis(data, session_type=roc_session_type, predictor=.predictor, shuf_iters=shuf_iters)
  
  ##Correlation analysis with continuous Y variables
  if (verbose){print("Performing correlation analysis with continuous response variables")}
  cor_df_torpor<-correlation_analysis(data, .session_type="torpor", response=c("temp", "temp_change1"), shuf_iters = shuf_iters, verbose=verbose)
  cor_df_torpor_only<-correlation_analysis(data%>%filter(temp<34), .session_type="torpor", response=c("temp", "temp_change1"), shuf_iters = shuf_iters, verbose=verbose)%>%
    rename_with(~ paste0(.x, "_TempBelow34"), contains("temp"))
  cor_df_ambient<-correlation_analysis(data, .session_type=c("cold","heat"), response=c("ambient_temp_interpolated","temp"), shuf_iters = shuf_iters, verbose=verbose)
  
  ##Lag correlation analysis
  if (verbose){print("Performing lag analysis during torpor")}
  max_cor_df_torpor<-unit_lag_analysis(telem_data = t_df, miniscope_data = data, lag_session_type = "torpor", response="temp", window=lag_window, shuf_iters = shuf_iters, verbose=verbose)
  
  ##Combine data and add metadata
  d<-merge(roc_df,cor_df_torpor,all=T)%>%merge(cor_df_torpor_only,all=T)%>%merge(cor_df_ambient,all=T)%>%merge(max_cor_df_torpor, all=T)
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
    fc_col<-paste0(type,"_fc")
    roc_df<-tibble(unit_id_id=character(0),"{auc_col}":=numeric(0), "{sig_col}":=character(0))
    if (type=="torpor"){roc_data<-data%>%filter(torpor_status=="non-torpor" | torpor_status=="deep_torpor")%>%mutate(label=ifelse(torpor_status=="non-torpor",0,1))}
    if (type=="heat"){roc_data<-data%>%filter(session_type=="heat")%>%filter(ambient_temp_bin=="22C" | ambient_temp_bin=="37C")%>%mutate(label=ifelse(ambient_temp_bin=="22C",0,1))}
    if (type=="cold"){roc_data<-data%>%filter(session_type=="cold")%>%filter(ambient_temp_bin=="22C" | ambient_temp_bin=="5C")%>%mutate(label=ifelse(ambient_temp_bin=="22C",0,1))}
    if (type=="male_interaction"){
      roc_data<-data%>%filter(session_type=="male_interaction", !is.na(male_interaction))%>%mutate(label=male_interaction)
      male_removal_minutes<-roc_data%>%filter(male_interaction==1)%>%group_by(session_id)%>%summarize(male_removal_minutes=max(session_time_minutes))
      roc_data<-roc_data%>%merge(male_removal_minutes,all.x=T)%>%filter(session_time_minutes < male_removal_minutes-1) #remove post-male time from ROC analysis
      }
    
    #Run analysis and return results
    for (id in unique(roc_data$unit_id_id)){
      #Real data
      d<-roc_data%>%filter(unit_id_id==id)
      if (length(unique(d$label)) != 2){next} #excludes sessions without both labels (will cause error in roc function)
      roc<-d%>%roc_("label",predictor, direction="<",levels=c(0,1))
      sum_roc<-d%>%group_by(label)%>%summarize(mean=mean(.data[[predictor]]))
      fc<-1+
          (sum_roc%>%filter(label==1)%>%pull(mean) - 
           sum_roc%>%filter(label==0)%>%pull(mean))
      
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
      roc_df<-rbind(roc_df, tibble(unit_id_id=id, "{auc_col}":=auc(roc), "{sig_col}":=auc_sig, "{fc_col}":=fc))
    }
    if (ii==1){compiled_df<-roc_df}
    if(ii>1){compiled_df<-merge(compiled_df,roc_df,all=T)}
    ii=ii+1
  }
  return(compiled_df)
}

correlation_analysis <- function(data, response, .session_type, predictor="df_f0_bin", shuf_iters=1000, method="pearson",verbose=T){
  type<-case_when(.session_type == "torpor" ~ "torpor",
                  "cold" %in% .session_type & "heat" %in% .session_type ~ "ambient")[1]
  if(verbose){print(paste0("Session type = ",type, ", Response = ",paste(response,collapse = ", ")))}
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

unit_lag_analysis <- function(telem_data, miniscope_data, window, response,lag_session_type="torpor", predictor="df_f0_bin", animal_var="mouse", timepoint_var="telem_ts", cor_method="pearson", verbose=T, shuf_iters=1000){
  ##Generate lagged telemetry data
  if(verbose){print("Creating lagged telemetry data")}
  i=1
  for (lag in window){
    new_name<-paste0(response,"_",lag)
    d<-telem_data%>%ungroup()%>%
      mutate("{timepoint_var}":= .data[[timepoint_var]] + minutes(lag))%>%
      select(all_of(c(timepoint_var,animal_var,response)))%>%
      rename("{new_name}":=.data[[response]])
    if(i==1){telem_data_w_lag<-d}
    if(i>1){telem_data_w_lag<-merge(telem_data_w_lag,d,all=T)}
    i=i+1
  }
  telem_data_w_lag<-telem_data_w_lag%>%drop_na()
  
  ##Add to miniscope data
  miniscope_d<-merge(miniscope_data%>%filter(session_type==lag_session_type), telem_data_w_lag)
  
  ##Calculate correlation for each unit at each lag value
  if(verbose){print("Computing correlations")}
  cor_df<-tibble()
  for (lag in window){
    col<-paste0(response,"_",lag)
    cor_col<-paste0(col,"_cor_",lag_session_type)
    cor_d<-correlation_analysis(miniscope_d, response=col, .session_type = lag_session_type, shuf_iters = shuf_iters, verbose=verbose)%>%
      mutate(lag=lag)%>%
      rename(cor = .data[[paste0(col,"_cor_",lag_session_type)]], cor_sig = .data[[paste0(col,"_cor_sig_",lag_session_type)]],slope=.data[[paste0(col,"_slope_",lag_session_type)]])%>%
      mutate(abs_cor=abs(cor))%>%
      select(!contains(col))
    cor_df<-rbind(cor_df,cor_d)
  }
  
  ##Find max correlation value for each unit at each lag value (filter to only significant correlations)
  strongest_cor_col<-paste0("strongest_",response,"_cor_",lag_session_type)
  strongest_lag_col<-paste0("strongest_",response,"_lag_",lag_session_type)
  sig_col<-paste0(response,"_lag_cor_sig",lag_session_type)
  
  max_cor_df<-cor_df%>%group_by(unit_id_id)%>%summarize("{strongest_cor_col}":=cor[which.max(abs_cor)], 
                                                        "{strongest_lag_col}":=lag[which.max(abs_cor)],
                                                        "{sig_col}":=cor_sig[which.max(abs_cor)])
  return(max_cor_df)
}

population_lag_analysis <- function(telem_data, miniscope_data, window, response,lag_session_type="torpor", predictor="df_f0_bin", animal_var="mouse", timepoint_var="telem_ts", cor_method="pearson", verbose=T, shuf_iters=1000){
  ##Generate lagged telemetry data
  if(verbose){print("Creating lagged telemetry data")}
  i=1
  for (lag in window){
    new_name<-paste0(response,"_",lag)
    d<-telem_data%>%ungroup()%>%
      mutate("{timepoint_var}":= .data[[timepoint_var]] + minutes(lag))%>%
      select(all_of(c(timepoint_var,animal_var,response)))%>%
      rename("{new_name}":=.data[[response]])
    if(i==1){telem_data_w_lag<-d}
    if(i>1){telem_data_w_lag<-merge(telem_data_w_lag,d,all=T)}
    i=i+1
  }
  telem_data_w_lag<-telem_data_w_lag%>%drop_na()
  
  ##Add to miniscope data
  miniscope_d<-merge(miniscope_data%>%filter(session_type==lag_session_type), telem_data_w_lag)
  
  ##Do lm analysis at each lag value
  if(verbose){print("Performing lm analysis")}
  cor_df<-tibble()
  coef_data<-tibble()
  for (lag in window){
    col<-paste0(response,"_",lag)
    cor_col<-paste0(col,"_cor_",lag_session_type)
    cor_d<-lm_analysis(miniscope_d, id_col = "telem_ts", response=col, .session_type = lag_session_type, shuf_iters = shuf_iters, verbose=verbose)
    lm_d<-(cor_d$lm_df)%>%
      mutate(lag=lag)%>%
      rename(cor = .data[[paste0(col,"_mean_cor_",lag_session_type,"_")]], cor_sig = .data[[paste0(col,"_cor_sig_",lag_session_type,"_")]])%>%
      select(!contains(col), !contains("shuf"))
    coef_d<-(cor_d$lm_coef_df)%>%
      mutate(lag=lag)%>%
      rename(weight = .data[[paste0(col,"_MeanLmCoef_",lag_session_type)]])%>%
      mutate(abs_weight=abs(weight))
    cor_df<-rbind(cor_df,lm_d)
    coef_data<-rbind(coef_data, coef_d)
  }
  ###########FINISH THIS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  ##Find max correlation value for each session at each lag value (filter to only significant correlations)
  strongest_cor_col<-paste0("strongest_",response,"_cor_",lag_session_type)
  strongest_lag_col<-paste0("strongest_",response,"_lag_",lag_session_type)
  sig_col<-paste0(response,"_lag_cor_sig",lag_session_type)
  
  max_cor_df<-cor_df%>%group_by(session_id)%>%summarize("{strongest_cor_col}":=cor[which.max(cor)], 
                                                        "{strongest_lag_col}":=lag[which.max(cor)],
                                                        "{sig_col}":=cor_sig[which.max(cor)])
  
  ##Find max correlation value for each unit at each lag value
  heaviest_weight_col<-paste0("heaviest_weight_",response,"_cor_",lag_session_type)
  heaviest_lag_col<-paste0("heaviest_lag_",response,"_cor_",lag_session_type)
  
  max_cor_unit_df<-coef_data%>%group_by(unit_id_id)%>%summarize("{heaviest_weight_col}":= abs_weight[which.max(abs_weight)],
                                                                "{heaviest_lag_col}":= lag[which.max(abs_weight)])
  
  return(list("temporal_session_df"=max_cor_df, "temporal_unit_df"=max_cor_unit_df))
}

transform_data_piegraph <- function(data, animal_var, cell_var){
  d<-data%>%filter(!is.na(.data[[cell_var]]))%>%
    group_by(across(all_of(c(animal_var,cell_var))))%>%
    count()%>%
    ungroup()%>%
    group_by(across(all_of(animal_var)))%>%
    mutate(percent=round((n/sum(n))*100, digits=1))
  return(d)
}

##Popoulation-level analysis
partition_data <- function(d, id_col, partition_type, partition_col, response, cv_folds=5, verbose=T){
  if (partition_type=="standard"){
    #Find number of bins to use for stratification of response variable (there must be more observations than cv_folds in each bin)
    bins<-1
    for (f in 2:cv_folds){
      min<-d%>%pull({{response}})%>%cut(f)%>%table()%>%min()
      if(min<cv_folds){break} #if bin size has gotten too small, stop here (and don't update bins variable)
      bins<-f
    }
    if (bins>1){pcs<-partition_cv_strat(d, coords=colnames(d),strat = d%>%pull({{response}})%>%cut(bins)%>%as.numeric()%>%as.factor(),nfold=cv_folds)} #use stratified partition when data has enough observations in each bin
    if (bins==1){ #use a non-stratified partition when data doesn't have enough observations
      pcs<-partition_cv(d, coords=colnames(d),nfold=cv_folds)
      if (verbose){print("Using non-stratified partition here")}
    } 
  }
  
  if (partition_type=="entry_arousal"){
    pcs<-partition_factor_cv(d, coords=colnames(d),fac=d%>%pull({{partition_col}})%>%droplevels(), nfold=2)
  }
  return(pcs)
}

safe_var <- function(x) paste0("`", x, "`")

lm_analysis <- function(data, .session_type, id_col, predictor = "df_f0_bin", partition_col=NULL, partition_type="standard", additional_x_var=NULL, response, cv_folds=5, shuf_iters=1000, verbose=T){
  if (partition_type=="standard"){folds<-1:cv_folds}
  if (partition_type=="entry_arousal"){folds<-cv_folds}
  lm_df<-tibble()
  lm_coef_df<-tibble()
  lm_add_x_var_coef_df<-tibble()
  predict_df<-tibble()
  data<-data%>%filter(session_type %in% .session_type)%>%ungroup()
  for (sid in unique(data$session_id)){
    ##Set up column names based on session type
    type<-case_when(.session_type=="torpor" ~ "torpor",
                    "heat" %in% .session_type & "cold" %in% .session_type ~ "ambient")[1]
    if(partition_type=="standard"){
      cor_col=paste0(response,"_mean_cor_",type,"_",additional_x_var)
      cor_shuf_col=paste0(response,"_mean_cor_shuf_",type,"_",additional_x_var)
      sig_col=paste0(response,"_cor_sig_",type,"_",additional_x_var)
      coef_col<-paste0(response,"_MeanLmCoef_",type,additional_x_var)
    }
    if(partition_type=="entry_arousal"){
      cor_col=paste0(response,"_mean_cor_EntryArousal_",type,"_",additional_x_var)
      sig_col=paste0(response,"_cor_sig_EntryArousal_",type,"_",additional_x_var)
      coef_col<-paste0(response,"_MeanLmCoef_EntryArousal_",type,additional_x_var)
      shuf_cor_col=paste0(response,"_mean_shuf_cor_EntryArousal_",type,"_",additional_x_var)
    }
    if(verbose){print(paste(sid, type, additional_x_var))}
    
    ##Format and pivot data
    d<-data%>%filter(session_id==sid)%>%
      select(unit_id_id, {{id_col}}, {{response}}, {{predictor}},{{additional_x_var}},{{partition_col}})%>%
      pivot_wider(values_from = {{predictor}}, names_from = unit_id_id)
    if (sum(is.na(d))>0)(warning(paste("NAs present in dataset",sid)))
    
    ##Partition data for cross-validation
    if (partition_type %nin% c("standard","entry_arousal")){stop("Partition type not recognized/supported")}
    pcs<-partition_data(d, id_col = id_col, response=response, partition_type=partition_type, partition_col=partition_col, cv_folds=cv_folds, verbose=verbose)
    
    ##Make models with observed data for each cv fold
    cor_df<-tibble()
    coef_df<-tibble()
    add_x_var_coef_df<-tibble()
    formula<-as.formula(paste0(safe_var(response),"~ . - ",id_col))
    for (fold in folds){
      #Set up data
      if(partition_type=="standard"){
        idx = pcs[["1"]][[fold]]
        train = d[idx$train,]
        test = d[idx$test,]}
      if(partition_type=="entry_arousal"){
        train=d%>%filter(.data[[partition_col]] == fold)%>%select(-{{partition_col}})
        test=d%>%filter(.data[[partition_col]] != fold)%>%select(-{{partition_col}})}
      
      #Train and test model
      model<-lm(formula, train)
      predict<-predict(model, newdata=test)
      cor<-cor(test%>%pull({{response}}), predict, method="pearson")
      cor_df<-rbind(cor_df, tibble("fold"=fold,"cor"=cor))
      coef_d<-tibble("unit_id_id"=names(model$coefficients),"coefficient"=unlist(model$coefficients),"fold"=fold)%>%filter(unit_id_id %nin% c("(Intercept)",additional_x_var))
      coef_df<-rbind(coef_df,coef_d)
      add_x_var_coef_d<-tibble("add_x_var"=names(model$coefficients),"coefficient"=unlist(model$coefficients),"fold"=fold,"session_id"=sid)%>%filter(add_x_var %in% additional_x_var)
      add_x_var_coef_df<-rbind(add_x_var_coef_df,add_x_var_coef_d)
      predict_df<-rbind(predict_df,tibble("predicted"=predict,"true"=test%>%pull({{response}}),"id"=test%>%pull({{id_col}}), "session_id"=sid,"fold"=fold))
    }
    
    ##Make models with shuffled data using same partitions for each cv fold
    shuf_cor_df<-tibble()
    shuf_add_x_var_coef_df<-tibble()
    for (i in 1:shuf_iters){
      shuf_cor_d<-tibble()
      for (fold in folds){
        #Set up data
        if(partition_type=="standard"){
          idx = pcs[["1"]][[fold]]
          train = (d[idx$train,])
          shuf_train = train%>%mutate("{response}":=sample(train%>%pull({{response}}),length(train%>%pull({{response}}))))
          test = d[idx$test,]}
        if(partition_type=="entry_arousal"){
          shuf_d = d%>%mutate("{partition_col}":=sample(d%>%pull({{partition_col}}), length(d%>%pull({{partition_col}}))))
          shuf_train = shuf_d%>%filter(.data[[partition_col]] == fold)%>%select(-{{partition_col}})
          test = shuf_d%>%filter(.data[[partition_col]] != fold)%>%select(-{{partition_col}})
        }
        
        #Train and test model
        model<-lm(formula, shuf_train)
        predict<-predict(model, newdata=test)
        cor<-cor(test%>%pull({{response}}), predict, method="pearson")
        shuf_cor_d<-rbind(shuf_cor_d,tibble("fold"=fold,"cor"=cor,"iteration"=i))
        shuf_add_x_var_coef_d<-tibble("add_x_var"=names(model$coefficients),"shuf_coefficient"=unlist(model$coefficients),"fold"=fold,"session_id"=sid,iter=i)%>%filter(add_x_var %in% additional_x_var)
        shuf_add_x_var_coef_df<-rbind(shuf_add_x_var_coef_df, shuf_add_x_var_coef_d)
      }
      if (partition_type=="standard"){shuf_cor_df<-rbind(shuf_cor_df, tibble("cor"=mean(shuf_cor_d$cor)))}
      if (partition_type=="entry_arousal"){shuf_cor_df<-rbind(shuf_cor_df, shuf_cor_d)}
    }
    
    if (partition_type=="standard"){
      ##Get rank of observed data within shuffled data
      mean_cor<-mean(cor_df$cor)
      rank<-(c(mean_cor,shuf_cor_df$cor)%>%rank())[1]
      cor_sig<-ifelse(rank>shuf_iters*0.95, "significant","non-significant")
      
      ##Compile data
      lm_d<-tibble("{cor_col}":= mean_cor, "{sig_col}":=cor_sig, "{cor_shuf_col}":=mean(shuf_cor_d$cor), "session_id"=sid)
      lm_df<-rbind(lm_df,lm_d)
      lm_coef_d<-coef_df%>%group_by(unit_id_id)%>%summarize("{coef_col}":=mean(coefficient))
      lm_coef_df<-rbind(lm_coef_df,lm_coef_d)
    }
    if (partition_type=="entry_arousal"){
      for (.fold in folds){
        ##Get rank of observed data within shuffled data
        mean_cor<-cor_df%>%filter(fold==.fold)%>%pull(cor)
        rank<-(c(mean_cor,shuf_cor_df%>%filter(fold==.fold)%>%pull(cor))%>%rank())[1]
        cor_sig<-case_when(rank>shuf_iters*0.975 ~ "improved",
                           rank<shuf_iters*0.025 ~ "worse",
                           T ~ "non-significant")
        
        ##Compile data
        lm_d<-tibble("train" = .fold,"{cor_col}":= mean_cor, "{shuf_cor_col}":=shuf_cor_df%>%filter(fold==.fold)%>%pull(cor)%>%mean(), "rank"=rank/shuf_iters, "{sig_col}":=cor_sig, "session_id"=sid)
        lm_df<-rbind(lm_df,lm_d)
      }
      lm_coef_d<-coef_df%>%group_by(unit_id_id,fold)%>%summarize("{coef_col}":=mean(coefficient))
      lm_coef_df<-rbind(lm_coef_df,lm_coef_d)
    }
    
    for (xvar in additional_x_var){
      #Get rank
      mean_obs_coef<-add_x_var_coef_df%>%filter(add_x_var==xvar)%>%group_by(add_x_var)%>%summarize(mean=mean(coefficient))%>%pull(mean)
      shuf_coef<-shuf_add_x_var_coef_df%>%group_by(iter)%>%summarize(shuf_mean=mean(shuf_coefficient))%>%pull(shuf_mean)
      coef_rank<-(c(mean_obs_coef,shuf_coef)%>%rank())[1]
      coef_sig<-case_when(coef_rank>shuf_iters*0.975 ~ "positive",
                          coef_rank<shuf_iters*0.025 ~ "negative",
                          T ~ "netural")
      #Compile
      xvar_coef_col<-paste0(xvar,"_MeanLmCoef_",type)
      xvar_sig_col<-paste0(xvar,"_sig_",type)    
      lm_add_x_var_coef_d<-add_x_var_coef_df%>%group_by(add_x_var,session_id)%>%summarize("{xvar_coef_col}":=mean(coefficient))%>%mutate("{xvar_sig_col}":=coef_sig)
      lm_add_x_var_coef_df<-rbind(lm_add_x_var_coef_df,lm_add_x_var_coef_d)
    }
  }
  if (partition_type=="entry_arousal"){
    lm_df<-lm_df%>%pivot_wider(id_cols=session_id, names_from=train,values_from=c(cor_col,shuf_cor_col,sig_col,rank))
    lm_coef_df<-lm_coef_df%>%pivot_wider(id_cols=unit_id_id,names_from=fold,values_from=coef_col)
    }
  
  ls<-list("lm_df" = lm_df, "lm_coef_df"=lm_coef_df, "predict_df"=predict_df, "lm_add_x_var_coef_df"=lm_add_x_var_coef_df)
  return(ls)
}

format_data_pca <- function(data, predictor, dim_to_reduce, id_dim){
  d<-data%>%
    select({{predictor}}, {{dim_to_reduce}}, {{id_dim}})%>%
    pivot_wider(values_from = {{predictor}}, names_from = {{dim_to_reduce}})%>%
    as.data.frame()
  rownames(d)<-as.character(d%>%pull({{id_dim}}))
  d<-d%>%select(!{{id_dim}})
  return(d)
}

plot_pca_var <- function(pca, .name){
  #Get pca_var
  pca.var <- pca$sdev^2
  pca.var.per <- round(pca.var/sum(pca.var)*100, 1)
  pcvardf<-tibble(PC = 1:length(pca.var.per),var=pca.var.per)
  
  #Plot
  p<-ggplot(pcvardf,aes(x=PC,y=var))+
    scale_x_continuous(breaks=seq(1,20,1),expand=c(0.01,0.01))+
    labs(y="% variance explained")+
    scale_y_continuous(expand=c(0,0))+
    geom_col()+
    ms
  save_plot(name = .name, plot=p, w=5,h=5)
  return(pcvardf)
}

pca <- function(data, predictor = "df_f0_bin", dims){
  if (length(dims)>2){stop("Length of dims exceeds 2")}
  pca_ls<-list()
  for (i in 1:length(dims)){
    pca_df<-tibble()
    reduce<-dims[i]
    id<-dims[-i]
    print(paste0("Reducing over ",reduce,", ID = ",id))
    pcvardf<-tibble()
    for (sidt in unique(data$session_id_type)){
      sidt_data<-data%>%filter(session_id_type==sidt)
      sid<-unique(sidt_data$session_id)
      if (grepl("cold",sidt)){sidt_data<-data%>%filter(session_id_type==paste0(sid,"_heat") | session_id_type==paste0(sid,"_cold"))}
      if (grepl("heat",sidt)){next} #prevents doubling up on ambient temperature data
      input_d<-format_data_pca(sidt_data, predictor = predictor, dim_to_reduce = reduce, id_dim=id)
      pca<-prcomp(input_d, center=F, scale=F)
      pcvardf<-rbind(pcvardf,plot_pca_var(pca, paste("PCA variance",id,sidt))%>%mutate(session_id_type=sidt))
      pca_data<-as.data.frame(pca$x)
      pca_data<-pca_data%>%mutate("{id}":= rownames(pca_data),session_id_type=sidt,session_id=sid)
      # pca_data<-merge(pca_data, sidt_data%>%ungroup()%>%mutate("{id}":=as.character(.data[[id]]))%>%distinct(!!sym(id), .keep_all = T),all.x=T)
      pca_df<-bind_rows(pca_df,pca_data)
    }
    pca_ls[[id]]<-pca_df
    pca_ls[[paste0(id,"_var")]]<-pcvardf
  }
  return(pca_ls)
}

plot_pca <- function(pca_ls, unit_df, verbose=T){

}

spatial_test <- function(data, var, shuf_iters=1000){
  cell_centroids <- data %>%
    filter(!is.na(.data[[var]])) %>%
    group_by(across(all_of(c(var,"unit_id_id","session_id")))) %>%
    summarize(x = sum(width * A) / sum(A), y = sum(height * A) / sum(A))
  observed<-cell_centroids%>%
    group_by(across(all_of(c(var,"session_id"))))%>%
    summarize(obs_dist = mean(nndist(x,y)))
  
  null_distribution <- map_dfr(1:shuf_iters, function(i) {
    cell_centroids %>%
      ungroup()%>%
      mutate("{var}":= sample(cell_centroids%>%pull({{var}}), length(cell_centroids%>%pull({{var}})))) %>% # Shuffle labels
      group_by(across(all_of(c(var,"session_id")))) %>%
      summarize(perm_dist = mean(nndist(x, y)))
  })%>%ungroup()
  null_distribution_sum<-null_distribution%>%group_by(across(all_of(c(var,"session_id"))))%>%summarize(null_dist=list(perm_dist))
  
  observed<-observed%>%
    merge(null_distribution_sum,all=T)%>%
    rowwise()%>%
    mutate(rank = rank(c(obs_dist,null_dist))[1],
           p = rank/shuf_iters,
           p.adj=p*length(unique(observed%>%pull({{var}}))),
           p.adj=ifelse(p.adj>1, 1, p.adj),
           sig=ifelse(p.adj<0.05, "p<0.05","ns"))%>%
    select(-null_dist)%>%
    ungroup()
  out<-list("observed"=observed, "null"=null_distribution)
  return(out)
}

##Telemetry
read_telemetry_data <- function(file, metadata_file, format = "starr-lifesci", idinfo = c("mouse", "misc", "measure", "misc2"), round = T, round_units = "minute", fstart_time="10:00:00", trial_length_days = 4){
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
    if (round){ #round data to nearest minute (only necessary if timestamps are not aligned on the minute)
      print(paste("Rounding telem_ts to nearest",round_units))
      d<-d%>%mutate(telem_ts_og = telem_ts,
                    telem_ts = round_date(telem_ts, round_units))
      } 
    if (any(grepl(d$measure, pattern="unknown"))){stop("Unknown measure in telemetry data")}
    d<-d%>%spread(measure,data)
    
  } else {
    stop("This function currently only supports 'starr-lifesci' as a format")
  }
  
  #Read in telemetry metadata and add to telemetry data
  md<-read_csv(metadata_file, show_col_types = F)%>%mutate(fstart = mdy_hms(paste(fstart, fstart_time))) #read metadata and convert fstart to date-time
  md<-md%>%mutate(group_gonad = factor(paste0(group,"_",gonad),levels=c("veh_intact","veh_ovx","e2_intact","e2_ovx")), #useful combination column
                  pellet = factor(pellet, levels=c("pre-OVX","OVX+Veh","OVX+E2"))) 
  trial1_end <- md%>%filter(trial=="trial1")%>%pull(fstart)%>%unique() + days(trial_length_days) #calculate end of trial1
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
  ## Read in events metadata csv and get path to directories for reading timeStamp files
  d<-read_csv(file, show_col_types=F)%>%
    filter(!is.na(webcam_frame))%>%
    mutate(across(everything(),as.character),
           path=paste(paste0(experiment,"_circulating_E2_torpor_miniscope"),timepoint,mouse,start_date,session,start_time,sep=separator))
  
  ## Get timestamps for WebCam and Miniscope
  webcam_ts_d<-tibble()
  miniscope_ts_d<-tibble()
  for (pth in unique(d$path)){
    d2<-d%>%filter(path==pth)
    
    #WebCam
    w_pth<-paste(pth,"My_WebCam","timeStamps.csv",sep=separator)
    t<-read_csv(w_pth, show_col_types = F)%>%
      rename(webcam_frame = `Frame Number`, webcam_time_ms = `Time Stamp (ms)`, buffer_index = `Buffer Index`)%>%
      mutate(exp=unique(d2$experiment),
             mouse=unique(d2$mouse),
             start_date=unique(d2$start_date),
             session=unique(d2$session),
             start_time=unique(d2$start_time))
    webcam_ts_d<-rbind(webcam_ts_d, t)
    
    #Miniscope
    m_pth<-paste(pth,"My_V4_Miniscope","timeStamps.csv",sep=separator)
    t<-read_csv(m_pth, show_col_types = F)%>%
      rename(frame = `Frame Number`, time_ms = `Time Stamp (ms)`, buffer_index = `Buffer Index`)%>%
      mutate(exp=unique(d2$experiment),
             mouse=unique(d2$mouse),
             start_date=unique(d2$start_date),
             session=unique(d2$session),
             start_time=unique(d2$start_time))
    miniscope_ts_d<-rbind(miniscope_ts_d, t)
  }
  
  ## Add webcam timestamp info back to d and add session_id as metadata to d and miniscope data
  d<-merge(d,webcam_ts_d,all.x=T)%>%mutate(session_id = paste(mouse,start_date,session,sep="_"), event_ts = ymd_hms(paste0(start_date,"_",start_time)) + milliseconds(webcam_time_ms))
  miniscope_ts_d<-miniscope_ts_d%>%mutate(session_id = paste(mouse,start_date,session,sep="_"))
  
  ## Find frame in miniscope data that is closest in timestamp to each event (deprecated, since it led to ambiguity with concatenated session timeStamps.csv)
  # for (i in 1:nrow(d)){
  #   frm<-(miniscope_ts_d%>%
  #           filter(session_id==d[[i,"session_id"]], start_time==d[[i,"start_time"]])%>%
  #           mutate(ts_diff = abs(time_ms-d[[i,"webcam_time_ms"]]))%>%
  #           filter(ts_diff == min(ts_diff))%>%
  #           pull(frame)%>%
  #           unique())[1] # [1] is important in case there is a "tie" between two frames for closest to the time stamp of interest
  #   d[i,"miniscope_frame"] <- frm
  # }
  if (NA %in% d$event_ts){warning("Not all event registered to a miniscope frame")}
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

point_cell <- function(...,alpha=0.3){
  point_indiv(..., size=1.5, position = position_jitter(width=0.25,height=0,seed=123),alpha=alpha)
}

point_mouse <- function(...){
  point_summary(...,aes(color=mouse),position=position_jitter(width=0.15,height=0,seed=123),size=3,stroke=1,shape=21,fill=NA,alpha=0.8)
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

xy_point3 <- function(...,size=3,shape=21,alpha=0.6, stroke=1.2){
  geom_point(size=size,shape=shape,alpha=alpha, stroke=stroke,...)
}

regression_line <- function(..., linewidth=1, method="lm", se=F){
  geom_line(stat="smooth",linewidth=linewidth,method = method, se = se,...)
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

ridgeline <- function(..., scale=0.9,linewidth=0.25,fill=NA){
  geom_ridgeline(scale=scale,linewidth=linewidth,fill=fill, ...)
}

ridgeline_guide <- function(..., var = "unit_id",color="grey"){
  geom_hline(aes(yintercept=!!sym(var)),color=color, ...)
}

get_torpor_y <- function(data){
  torpor_y_increment<-case_when(range(data$temp)%>%diff() < 1.5 ~ 0.5,
                                range(data$temp)%>%diff() < 3 ~ 1,
                                range(data$temp)%>%diff() < 6 ~ 2,
                                range(data$temp)%>%diff() < 9 ~ 3,
                                range(data$temp)%>%diff() < 12 ~ 4,
                                T ~ 5)
  torpor_y<-if(torpor_y_increment %in% c(1,2,3,4,6)){scale_y_continuous(breaks=seq(19,43,torpor_y_increment),expand=c(0.1,0.1))}else{scale_y_continuous(breaks=seq(16,46,torpor_y_increment),expand=c(0.1,0.1))}
  return(torpor_y)
}


moving_avg <- function(formula, data, window = 1, ...) {
  structure(list(data = data, window = window), class = "moving_avg")
}

predict.moving_avg <- function(object, newdata, se.fit = FALSE, ...) {
  x <- object$data$x
  y <- object$data$y
  w <- object$window
  
  results <- lapply(newdata$x, function(xi) {
    in_window <- x >= (xi - w/2) & x <= (xi + w/2)
    yvals <- y[in_window]
    n <- sum(!is.na(yvals))
    m <- mean(yvals, na.rm = TRUE)
    s <- sd(yvals, na.rm = TRUE)
    sem <- if (n > 1) s / sqrt(n) else NA_real_
    list(fit = m, sem = sem)
  })
  
  fit <- sapply(results, `[[`, "fit")
  sem <- sapply(results, `[[`, "sem")
  
  if (se.fit) {
    fit_matrix <- cbind(
      fit = fit,
      lwr = fit - sem,   # mean - 1 SEM
      upr = fit + sem    # mean + 1 SEM
    )
    list(fit = fit_matrix, se.fit = sem)
  } else {
    fit
  }
}

moving_avg_trail <- function(formula, data, window = 1, ...) {
  structure(list(data = data, window = window), class = "moving_avg_trail")
}

predict.moving_avg_trail <- function(object, newdata, se.fit = FALSE, ...) {
  x <- object$data$x
  y <- object$data$y
  w <- object$window
  
  results <- lapply(newdata$x, function(xi) {
    in_window <- x >= (xi - w) & x <= (xi)
    yvals <- y[in_window]
    n <- sum(!is.na(yvals))
    m <- mean(yvals, na.rm = TRUE)
    s <- sd(yvals, na.rm = TRUE)
    sem <- if (n > 1) s / sqrt(n) else NA_real_
    list(fit = m, sem = sem)
  })
  
  fit <- sapply(results, `[[`, "fit")
  sem <- sapply(results, `[[`, "sem")
  
  if (se.fit) {
    fit_matrix <- cbind(
      fit = fit,
      lwr = fit - sem,   # mean - 1 SEM
      upr = fit + sem    # mean + 1 SEM
    )
    list(fit = fit_matrix, se.fit = sem)
  } else {
    fit
  }
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

get_mode <- function(x) {
  ux <- unique(na.omit(x))
  ux[which.max(tabulate(match(x, ux)))]
}
