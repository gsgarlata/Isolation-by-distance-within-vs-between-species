#This script contains functions that are used for computing genetic distances.


#'countDiff' is used to count genetic differences for one SNP, given the genotypes of two individuals.
#Diploid genotypes are coded as 0 (no copies of the alternative allele), 1 (one reference and one alternative allele) or 2 (no copies of the reference allele).
#Pairs of individual genotypes are coded, for instance, as '1-0', where one individual is heterozygote ('1') and the other individual is homozygote for the reference allele ('0').
#See Korunes and Samuk, (2021) to understand how the number of pairwise differences are computed between genotypes.
countDiff<-function(x){
  
  #One of the two individuals is heterozygote. The total number of differences is 2.
  if(x=='1-0' | x=='0-1'){res = 2}
  
  #Both individuals are homozygote for the reference allele. The total number of differences is 0.
  if(x=='0-0'){res = 0}
  
  #One individual is homozygote for the alternative allele ('2') and the other individual is homozygote for the reference allele.
  #The total number of differences is 4.
  if(x=='2-0' | x=='0-2'){res = 4}
  
  #Both individuals are homozygote for the alternative allele. The total number of differences is 0.
  if(x=='2-2'){res = 0}
  
  #One individual is homozygote for the alternative allele ('2') and the other individual is heterozygote ('1').
  #The total number of differences is 2.
  if(x=='1-2' | x=='2-1'){res = 2}
  
  #Both individuals are heterozygotes ('1'). The total number of differences is 2.
  if(x=='1-1'){res = 2}
  
  #If any of the two individuals has a missing genotype for the SNP, retuns 'NA'.
  if(length(grep('NA',x))>0){res = 'NA'}
  
  return(res)
}

#'makeCountDIff' is used to compute the expected number of differences for a randomly sampled variant site.
#Specifically, it computes the pairwise differences between genotypes, across all sites and then estimate the expected number of differences 
#by dividing over the total number of pairwise comparisons across all sites.
makeCountDIff<-function(pair_df){
  
  #It creates a matrix of size
  count_mat = matrix(NA,nrow = nrow(pair_df),ncol = 2,dimnames = list(1:nrow(pair_df),c('combo','count')))
  pair_df = cbind(pair_df,count_mat)
  
  pair_df[,3] = paste0(pair_df[,1],'-',pair_df[,2])
  
  #It computes the number of pairwise differences for each genotype across all site.
  pair_df[,4] = as.numeric(sapply(pair_df[,3],countDiff))
  
  pair_df = as.data.frame(pair_df)
  pair_df$count = as.numeric(pair_df$count)
  
  #It removes the sites where one of the two individuals has missing data.
  coordNA = which(is.na(pair_df$count))
  if(length(coordNA)==0){}else{
  pair_df=pair_df[-coordNA,]
  }
  
  #It sum the pairwise differences across all sites
  sum_diff = sum(pair_df$count,na.rm = T)
  #Since we have removed all the sites for which at least one of the two individuals has missing data,
  #we can consider that the total number of comparisons between two diploid genotypes is 6 (see Korunes and Samuk, 2021).
  pair_df$comparisons = 6
  #It sums the total number of pairwise comparisons across sites.
  tot_comparisons = sum(pair_df$comparisons)
  #It compute the percentage of pairwise differences between genotypes across all sites. It provides the epxected number of differences
  #for a variant site (not-missing for any of the two individuals).
  perc_diff = sum_diff/tot_comparisons
  return(data.frame(perc_diff=perc_diff,sum_diff=sum_diff,tot_comparisons=tot_comparisons))
}

#'windowIBDcomputationPIXY' is used to compute individual pairwise genetic distances within and between-taxa.

#'snps': it is the R object including all the 'splitted' genomic windows.
#'window_index': the index of the focal genomic window. (i.e. 1 specifies the first genomic window).
#'species_abbr': the name of the focal candidate species pair.
#'IDs': a vector of individual IDs for all the individuals (across all species) included in the analyses.
#'geo_df': it is the dataframe containing the geographic distances of each individual pairwise comparison in the dataset.
#'popfile': the population file specifies the taxon of each individual.
windowIBDcomputationPIXY<-function(snps,window_index,species_abbr,IDs,geo_df,popfile){
  
  #It reads the 'population' file.
  popdata = read.csv(popfile)
  #It removes missing data.
  popdata = popdata[complete.cases(popdata),]
  
  #It retrieves the 'window_index' genomic window from the 'snps' R object.
  w2_df = snps(window_index)
  #It assign the individuals IDs to the columns of the 'w2_df' R object.
  colnames(w2_df) = IDs
  #It subset the candidate species pair from the population file 
  sub_popdata = popdata[popdata$cluster %in% species_abbr,]
  #It extract the individuals of a given candidate species pair from 'w2_df'
  w2_sps = w2_df[,colnames(w2_df) %in% sub_popdata$individuals]
  
  ###START - It computes average number of genotype pairwise differences among individuals within species.
  tot_within = NULL
  
  #It loops over each taxon of a given candidate species pair, to compute the average number of pairwise differences
  #across all individual pairwise comparisons.
  for(i in 1:length(species_abbr)){
    
    #Extract the 'population' information of individuals for a given taxon in the candidate species pair. 
    coord_id = grep(species_abbr[i],sub_popdata$cluster)
    within_sp = w2_sps[,colnames(w2_sps) %in% sub_popdata$individuals[coord_id]]
    
    #It lists all individual pairwise comparisons within a given taxon of candidate species pair.
    combo_pairs = t(combn(colnames(within_sp),2))
    
    species_df = NULL
    
    #It loops over all pairwise individual comparison within a taxon.
    for(j in 1:nrow(combo_pairs)){
      #It extract an individual pair from the 'snp' R object.
      pair_df = within_sp[,colnames(within_sp) %in% combo_pairs[j,]]
      #It computes the average number of pairwise differences of the genotypes of a pair of individuals.
      perc_df = makeCountDIff(pair_df)
      #It store all the relevant information in a temporary dataframe.
      tmp_df = data.frame(id1=combo_pairs[j,1],id2=combo_pairs[j,2], sp1=species_abbr[i],sp2=species_abbr[i],
                          value=perc_df$perc_diff,count_diffs=perc_df$sum_diff,comparisons=perc_df$tot_comparisons)
      species_df = rbind(tmp_df,species_df)
      rm(tmp_df)
    }
    
    #It classify these pairwise comparisons as 'within-taxon'.
    species_df$type = 'Within'
    species_df$combo = paste0(species_abbr[i],'-',species_abbr[i])
    tot_within = rbind(species_df,tot_within)
  }
  ###END - It computes average number of genotype pairwise differences among individuals within species.
  
  ###START - It computes average number of genotype pairwise differences among individuals between species.
  #It lists all individual pairwise comparisons between taxa of a candidate species pair.
  combo_pairs = t(combn(colnames(w2_sps),2))
  combo_pairs_df = data.frame(id1=combo_pairs[,1],id2=combo_pairs[,2])
  
  #It creates two columns, one for each taxon of the candidate species pair.
  combo_pairs_df$sp1=sub_popdata$cluster[match(combo_pairs[,1],sub_popdata$individuals)]
  combo_pairs_df$sp2=sub_popdata$cluster[match(combo_pairs[,2],sub_popdata$individuals)]
  #It subsets only pairwise individuals from different taxa.
  combo_pairs_df = combo_pairs_df[!combo_pairs_df$sp1==combo_pairs_df$sp2,]
  
  combo_pairs_df$value = NA
  combo_pairs_df$count_diffs = NA
  combo_pairs_df$comparisons = NA
  
  #It loops over all pairwise individual comparisons between taxa of the candidate species pair.
  for(j in 1:nrow(combo_pairs_df)){
    #It extract an individual pair from the 'snp' R object.
    pair_df = w2_sps[,colnames(w2_sps) %in% combo_pairs_df[j,]]
    #It computes the average number of pairwise differences of the genotypes of a pair of individuals.
    perc_diff = makeCountDIff(pair_df)
    combo_pairs_df$value[j] = perc_diff$perc_diff
    combo_pairs_df$count_diffs[j] = perc_diff$sum_diff
    combo_pairs_df$comparisons[j] = perc_diff$tot_comparisons
  }
  
  #It classify these pairwise comparisons as 'between-taxa'.
  combo_pairs_df$type = 'Between'
  tot_between = combo_pairs_df
  tot_between$combo = NA
  
  #It order the species names by alphabetic order.
  for(k in 1:nrow(tot_between)){
    ord_sps=sort(c(tot_between$sp1[k],tot_between$sp2[k]))
    tot_between$combo[k] = paste0(ord_sps[1],'-',ord_sps[2])
  }
  ###END - It computes average number of genotype pairwise differences among individuals between species.
  
  tot_final = rbind(tot_within,tot_between)
  tot_final$id1 = as.character(tot_final$id1)
  tot_final$id2 = as.character(tot_final$id2)
  
  tot_final$id_combo = NA
  
  for(k in 1:nrow(tot_final)){
    ord_ids = order(c(tot_final$id1[k],tot_final$id2[k]))
    tot_final$id_combo[k] = paste0(tot_final[k,ord_ids[1]],'-',tot_final[k,ord_ids[2]])
  }
  
  #It extract the geographical distance among each individual pairwise comparison (between or within-taxa).
  tot_final$geo_dist = geo_df[match(tot_final$id_combo,geo_df[,'id_combo']),'value']
  
  return(tot_final)
}

#'windowIBDcomputationPIXYwithinBetwenSP' is equivalent to 'windowIBDcomputationPIXY' but only for M. tavaratra (the species used as null model for spatial population structure).
windowIBDcomputationPIXYwithinBetwenSP<-function(snps,window_index,species_abbr,IDs,geo_df,popfile){
  
  popdata = read.csv(popfile)
  ####
  w2_df = snps(window_index)
  colnames(w2_df) = IDs
  
  coord_total = sort(grep(species_abbr,colnames(w2_df)))
  w2_sps = w2_df[,coord_total]
  
  ###START - within species
  coord_id = grep(species_abbr,colnames(w2_sps))
  within_sp = w2_sps[,coord_id]
  
  combo_pairs = t(combn(colnames(within_sp),2))
  
  tot_within = NULL
  
  for(j in 1:nrow(combo_pairs)){
    pair_df = within_sp[,colnames(within_sp) %in% combo_pairs[j,]]
    perc_df = makeCountDIff(pair_df)
    tmp_df = data.frame(id1=combo_pairs[j,1],id2=combo_pairs[j,2], sp1=species_abbr,sp2=species_abbr,
                        value=perc_df$perc_diff,count_diffs=perc_df$sum_diff,comparisons=perc_df$tot_comparisons)
    tot_within = rbind(tmp_df,tot_within)
  }
  
  popdata$id_combo = NA
  for(a in 1:nrow(popdata)){
    popdata$id_combo[a] = paste(sort(c(popdata$id1[a],popdata$id2[a])),collapse='-')
  }
  
  tot_within$id_combo = NA
  
  for(b in 1:nrow(tot_within)){
    tot_within$id_combo[b] = paste(sort(c(tot_within$id1[b],tot_within$id2[b])),collapse='-')
  }
  
  tot_within$type = popdata$type[match(tot_within$id_combo,popdata$id_combo)]
  tot_within$combo = paste0(species_abbr,'-',species_abbr)
  ###END - within species
  
  tot_within$id1 = as.character(tot_within$id1)
  tot_within$id2 = as.character(tot_within$id2)
  
  tot_within$geo_dist = geo_df[match(tot_within$id_combo,geo_df[,'id_combo']),'value']
  
  return(tot_within)
}


nrme_regression<-function(df,list_pairs){
  pairs_sp = do.call(rbind,strsplit(list_pairs,'-'))
  
  coord_with = which(pairs_sp[,1]==pairs_sp[,2])
  coord_betw = which(!pairs_sp[,1]==pairs_sp[,2])
  
  final_within = NULL
  ###within species, separately
  for(s in 1:length(coord_with)){
    
    sp_within_id = list_pairs[coord_with[s]]
    sp_within_df = subset(df,combo==sp_within_id)
    lm.within = lm(value ~ geo_dist,data=sp_within_df)
    
    sp_between_id = list_pairs[coord_betw]
    sp_between_df = subset(df,combo==sp_between_id)
    
    
    yy<-as.vector(predict.lm(lm.within, sp_between_df))
    nrmse<-sqrt(mean((yy - sp_between_df$value)^2,na.rm=TRUE)) / diff(range(sp_between_df$value,na.rm = TRUE))
    
    tmp_within<-data.frame(within=sp_within_id,between=sp_between_id,nrmse=nrmse)
    final_within<-rbind(final_within,tmp_within)
  }
  
  return(final_within)
}


