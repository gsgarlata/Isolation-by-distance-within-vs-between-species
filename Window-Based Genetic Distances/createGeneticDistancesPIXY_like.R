##This script is used to compute 'Isolation-by-distance' (IBD) between candidate species-pairs for each genomic window, separately.
#It relies on the "vcf_windower" function of the "locustr" R package by Li and Ralph, 2019. Please, if you use this scripts,
#cite both our study and Li and Ralph, (2019).


require(lostruct)
require(ggplot2)

###FUNCTIONS###
matrixConvertLower<-function (triMatrix, colname = c("sp1", "sp2", "dist")) {
  m <- as.matrix(triMatrix)
  m2 <- reshape2::melt(m)[reshape2::melt(lower.tri(m))$value, 
  ]
  names(m2) <- colname
  invisible(m2)
}
###FUNCTIONS###

#Input arguments
path_input = args[1]
path_info = args[2]
path_out = args[3]
path_bin = args[4]

geo_file = 'geo_dist_matrix'
pop_name = 'popfile.csv'
genomic_fle = 'mydata.bcf'
window_size = 1000

source(paste0(path_bin,'/functions_regressionPIXYlike.R'))

#it defines the population/species classification for each individual
popfile = paste0(path_info,'/',pop_name)
#it defines the genomic file (.bcf format)
input_file = paste0(path_input,"/",genomic_fle)

##START: Get GEOGRAPHIC DISTANCES matrix##
#read a matrix of pairwise geographic distances between all mouse lemurs individuals.
geo_df = read.table(paste0(path_input,'/',geo_file,'.txt'))
#convert the matrix in a three-columns dataframe.
dat.geo<-matrixConvertLower(geo_df,colname=c('id1','id2','value'))
dat.geo$id1 = as.character(dat.geo$id1)
dat.geo$id2 = as.character(dat.geo$id2)

#it order the names' pairs of the individuals, following alphabetic order.
dat.geo$id_combo = NA
for(q in 1:nrow(dat.geo)){
  ord_ids = order(c(dat.geo$id1[q],dat.geo$id2[q]))
  dat.geo$id_combo[q] = paste0(dat.geo[q,ord_ids[1]],'-',dat.geo[q,ord_ids[2]])
}
##END: Get GEOGRAPHIC DISTANCES matrix##

#it reads the genomic data input file (.bcf format) using the 'vcf_windower' function from Li and Ralph (2019), 
#which splits the data in genomic windows of fixed size 'window_size' (in terms of number of SNPs and not physical position).
#For instance, a 'window_size' of 1000 means that each window will contain 1000 SNPs. 
snps <- vcf_windower(file,size=window_size,type='snp')
#extract the ID of the individuals.
IDs = attr(snps,"samples")
#extract the number of splitted genomic windows.
n_windows = attr(snps,"max.n")

#List all the pairwise taxon comparisons.
species_list = list(c('MmurC','MmurN'),c('Mgan','MmurN'),c('Mgan','MmurC'),
                    c('Mman','MmurN'), c('Mman','MmurC'), c('Mtav'),
                    c('Marn','Msp1'), c('Mger','Mjol'), c('Mmac','Mjon'),
                    c('Mmit','Mleh'),c('Mman','Mmur'), c('Mgan','Mmur'),
                    c('Mgan','Mman'),c('Mmyo','Mruf'),c('Mmyo','Mber'),
                    c('Mber','Mruf'), c('Mmarg','Msam'), c('Mmam','Msam'),
                    c('Mmam','Mmarg'),c('Msim','Mbor'),c('Mbon','Mdan'),
                    c('Mrav','Mbon'),c('Mrav','Mdan'))


for(r in 1:length(species_list)){
  
  species_abbr = species_list[[r]]
  
  ####ONLY LINEAR REGRESSION WITH NRMSE
  final_window = NULL
  
  #loop over each genomic window
  for(j in 1:n_windows){
  
  #If conditioning is used only for the Mtav case, for which 'fragmented' and 'continuous' populations are distringuished for the purposes
  #of the null model of spatial population structure.
  if(length(species_abbr)==1){
    #It computes IBD between and within taxon, separately, for each genomic window. This function is used specifically for 'Mtav'.
    tmpDF = windowIBDcomputationPIXYwithinBetwenSP(snps,window_index=j,species_abbr=species_abbr,IDs=IDs,geo_df=dat.geo,popfile=popfile)
      }else{
    #It computes IBD between and within taxon, separately, for each genomic window.
    tmpDF = windowIBDcomputationPIXY(snps,window_index=j,species_abbr=species_abbr,IDs=IDs,geo_df=dat.geo,popfile=popfile)
      }
    #It add the window index in the columnn 'window'.
    tmpDF$window = j

  final_window = rbind(final_window,tmpDF)
  rm(tmpDF)
  }
  
  if(length(species_abbr)==1){
  filename = paste0(path_out,'/',species_abbr,'_versus_',species_abbr,'_allwindows')
  }else{
  filename = paste0(path_out,'/',species_abbr[1],'_versus_',species_abbr[2],'_allwindows')
  }
  
write.table(file=paste0(filename,'_dataPIXY.txt'),as.matrix(final_window),quote=FALSE)

}

      
      
      
  