predict_combined_original = function(object, newdata2, process = "event",Tstart,
                                     times, return_newdata = TRUE){
  
  id_var <- object$model_info$var_names$idVar
  obs_time_var <- object$model_info$var_names$time_var
  stat_time_var <- object$model_info$var_names$Time_var
  event_var <- object$model_info$var_names$event_var
  
  
  newdata2 <- newdata2[newdata2[[obs_time_var]] <= Tstart, ]
  newdata2[[stat_time_var]] <- Tstart
  newdata2[[event_var]] <- 0
  all_ID = unique(newdata2[[id_var]])
  
  
  pred_res = lapply(1:length(all_ID),function(kk){
    ND <- newdata2[newdata2[[id_var]]==all_ID[kk], ]
    ND = ND %>% arrange(!!sym(obs_time_var))
    if(nrow(ND)<1){
      return(NULL)
    }else{
      predSurv <- predict(object, newdata = ND, process = "event",
                          times = times,
                          return_newdata = TRUE)
      return(predSurv)
    }
  }) 
  
  predSurv_res = lapply(pred_res,function(uu){
    tmp = uu %>% dplyr::select(c(pred_CIF,low_CIF,upp_CIF),any_of(c(id_var,stat_time_var)))
    return(tmp)
  }) %>% do.call(rbind,.)
  
  return(predSurv_res)  
}
