#This script is used to compute the Normalised Root Mean Squared Error (NMRSE) distribution for the between-species individual pairwise comparisons.
#In particular, 'observed' between-species individual genetic distances are compared to the corresponding 'predicted' individual genetic distances,
#obtained using the 'between-species' geographic distances and the within-species IBD model parameters.


require(ggplot2)

path_input = args[1]
path_output = args[2]
window_size = args[3]

source(paste0(path_bin,'/functions_LOGnrmse.R'))

species_list = list(c('MmurC','MmurN'),c('Mgan','MmurN'),c('Mgan','MmurC'),
                    c('Mman','MmurN'), c('Mman','MmurC'), c('Mtav'),
                    c('Marn','Msp1'), c('Mger','Mjol'), c('Mmac','Mjon'),
                    c('Mmit','Mleh'),c('Mman','Mmur'), c('Mgan','Mmur'),
                    c('Mgan','Mman'),c('Mmyo','Mruf'),c('Mmyo','Mber'),
                    c('Mber','Mruf'), c('Mmarg','Msam'), c('Mmam','Msam'),
                    c('Mmam','Mmarg'),c('Msim','Mbor'),c('Mbon','Mdan'),
                    c('Mrav','Mbon'),c('Mrav','Mdan'))

#It loops over all candidate species-pairs. 
for(r in 1:length(species_list)){
  
  species_abbr = species_list[[r]]
  
  if(length(species_abbr)==1){
  filename = paste0(path_input,'/',species_abbr,'_versus_',species_abbr,'_allwindows')
  }else{
  filename = paste0(path_input,'/',species_abbr[1],'_versus_',species_abbr[2],'_allwindows')
  }
  
  #It reads the dataframe with the average individual pairwise genetic distances for each genomic window. 
  df_data = read.table(file=paste0(filename,'_dataPIXY.txt'))
  
  #It create a vector of window indeces.
  wind_list = unique(df_data$window)
  
  final_lognrmse = NULL
  
  #It loops over each genomic window.
  for(wind_index in wind_list){
    
    #It subset a specific genomic window.
    sub_df = subset(df_data,window==wind_index)
    #It identify missing data in the geographic distance dataframe and remove them.
    coordNA = which(is.na(sub_df$geo_dist))
    if(length(coordNA)==0){}else{
    sub_df = sub_df[-which(is.na(sub_df$geo_dist)),]
    }
    #It lists all individual pairwise comparisons 
    list_pairs = unique(sub_df$combo)
    
    #It checks whether there are individual pairwise geographic distances equal to 0, and if there are, it adds a 'fictitious 0.01' in order
    #to get a definite number once the logarithm of geographical distance is calculated (i.e., you cannot do the log of 0). 
    c_epsilon = as.numeric(quantile(sub_df$geo_dist,na.rm = T)[2])
   if(c_epsilon==0){
     c_epsilon = 0.01
   }else{}
   
   #It computes the logarithm of geographic distance.
   sub_df$log_geo_dist = log(sub_df$geo_dist + c_epsilon)
   
   if(length(species_abbr)==1){
   lognrmse = LOGnrme_regressionWithinSP(sub_df)
   }else{
   lognrmse = LOGnrme_regression(sub_df,list_pairs)
   }
   
   lognrmse$window = wind_index
   final_lognrmse = rbind(final_lognrmse,lognrmse)
  }
  
  if(length(species_abbr)==1){
  filename = paste0(path_output,'/',species_abbr,'_versus_',species_abbr,'_allwindows_windsize',window_size)
  }else{
  filename = paste0(path_output,'/',species_abbr[1],'_versus_',species_abbr[2],'_allwindows_windsize',window_size)
  }
  write.table(file=paste0(filename,'_LOGnrmsePIXY.txt'),as.matrix(final_lognrmse),quote=FALSE)

  }


