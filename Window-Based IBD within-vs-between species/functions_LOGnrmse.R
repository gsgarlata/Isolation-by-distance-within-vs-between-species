#This script contains functions that are used for computing genetic distances.


#'LOGnrme_regression' is used to compute the Normalised Root Mean Squared Error (NMRSE) distribution for the between-species individual pairwise comparisons,
#using the 'within-species' IBD model to obtain the predicted genetic distances, which are compared to the observed between-species individual genetic distances.
#'df': a dataframe containing individual pairwise comparisons within- and between-species of a given candidate species-pair.
#'list_pairs': within- and between-species comparisons for a given candidate species-pair (i.e., ).
LOGnrme_regression<-function(df,list_pairs){
  
  #It splits the strings containing the species pairs labels. 
  #'list_pairs' is a three elements vector such as: c("Mdan-Mdan" "Mrav-Mrav" "Mdan-Mrav").
  pairs_sp = do.call(rbind,strsplit(list_pairs,'-'))
  
  #It identifies the within-species groups (e.g., c("Mdan-Mdan" "Mrav-Mrav")).
  coord_with = which(pairs_sp[,1]==pairs_sp[,2])
  #It identifies the between-species group (e.g., c("Mdan-Mrav")).
  coord_betw = which(!pairs_sp[,1]==pairs_sp[,2])
  
  final_within = NULL
  
  #It loops over each "within-species" group (e.g., c("Mdan-Mdan" "Mrav-Mrav")).
  for(s in 1:length(coord_with)){
    
    #It subsets one of the two within-species groups.
    sp_within_id = list_pairs[coord_with[s]]
    #It subsets the 'df' dataframe for a specific within-species group.
    sp_within_df = subset(df,combo==sp_within_id)
    
    #It performs linear regression between genetic distance and the logarithm of geographical distance
    #for the within-species group pairiwse comparisons.
    lm.within = lm(value ~ log_geo_dist,data=sp_within_df)
    
    #It subsets the between-species group.
    sp_between_id = list_pairs[coord_betw]
    #It subsets the 'df' dataframe for the between-species group.
    sp_between_df = subset(df,combo==sp_between_id)
    #It creates a vector of predicted genetic distances for 'between-candidate species' individual pairwise comparisons using 
    #the 'within-species' IBD model parameters (fitted parameters in the linear regression 'lm.within').
    yy<-as.vector(predict.lm(lm.within, sp_between_df))
    #It computes the Normalised Root Mean Squared Error (NMRSE) between observed and 
    #predicted 'between-species' individual pairwise genetic distances.
    nrmse<-sqrt(mean((yy - sp_between_df$value)^2,na.rm=TRUE)) / diff(range(sp_between_df$value,na.rm = TRUE))
    
    tmp_within<-data.frame(within=sp_within_id,between=sp_between_id,nrmse=nrmse)
    final_within<-rbind(final_within,tmp_within)
  }
  
  return(final_within)
}

#'LOGnrme_regressionWithinSP' is equivalent to 'LOGnrme_regression' but for the null-model species ('M. tavaratra')
#for which we treat 'fragmented' and 'continuous' populations as they were different species.
LOGnrme_regressionWithinSP<-function(df){
  
  coord_with = which(df$type=='Within')
  coord_betw = which(df$type=='Between')
  
  sp_within_df = df[coord_with,]
  lm.within = lm(value ~ log_geo_dist,data=sp_within_df)
  
  sp_between_df = df[coord_betw,]
  
  yy<-as.vector(predict.lm(lm.within, sp_between_df))
  nrmse<-sqrt(mean((yy - sp_between_df$value)^2,na.rm=TRUE)) / diff(range(sp_between_df$value,na.rm = TRUE))
  
  final_within<-data.frame(within=unique(df$combo),between=paste0(unique(df$sp1),'C-',unique(df$sp1),'F'),nrmse=nrmse)
  
  return(final_within)
}
