#This script allows to compute the signal/error ratio for a specific window size. It follows the procedure described
#in (Li and Ralph, 2019), please check their article for more details.


require(lostruct)

path_bin = args[1]
path_data = '/Users/gabrielemariasgarlata/Documents/Projects/Ongoing_projects/Microcebus_Phylogenomics_Consortium/Microcebus_phylogeny/analysis/IBD/output/signal2error/wind1000' #args[2]
win_size = as.numeric(args[3])
jack_block = as.numeric(args[4])
filename = args[4]
path_out = args[5]

#It loads scripts needed for compute 'signal' and 'error' statistics.
source(paste0(path_bin,'/functions_signal2error_PCs.R'))

#It lists the files (one for each chromosome) containing the PCA output of the genomic windows laying in that chromosome.
list_pc_files = list.files(path=path_data,pattern='.txt')
#list_pc_files = list_pc_files[-c(1,2)]
###START: create dataset of PCs across all genomic windows###
final_df = NULL
final_snp = NULL

#It loops over each PC file and R object containing the respective genomic data.
for(i in 1:length(list_pc_files)){
  
  pre_file = gsub('\\_pcs.txt','',list_pc_files[i])
  #It reads the PC file for a specific chromosome
  txt_file = paste0(path_data,'/',pre_file,'_pcs.txt')
  tmp_txt_data = read.csv(txt_file)
  
  #It reads the genomic data for the respective chromosome
  snps_file = paste0(path_data,'/',pre_file,'_snps.RData')
  tmp_snp_data = readRDS(snps_file)

  #It combines the SNP data of each chromosome (including multiple genomic windows) in a unique R object ('final_snp').
  final_snp = c(final_snp,tmp_snp_data)
  
  #This 'if conditioning' is used to reformat the datasets.
  if(ncol(tmp_txt_data) > 2){
    tmp_txt_data[,1] = NULL
  }else{
    rownames(tmp_txt_data) = tmp_txt_data[,1]
    tmp_txt_data[,1] = NULL
    tmp_txt_data = t(tmp_txt_data)
  }
  
  #It binds the reformatted datasets in a unique R object ('final_df').
  final_df = rbind(final_df,tmp_txt_data)
  rm(tmp_txt_data, tmp_snp_data)
}

rownames(final_df) = NULL
final_df = as.matrix(final_df)
###END: create dataset of PCs###

#It uses a function included in the 'functions_signal2error_PCs.R' script to compute the SIGNAL variance in the dataset,
#following Li and Ralph, (2019).
signalDataframe = computeVarSignal(data_df=final_df)

n_wind = nrow(final_df)

#It uses a function included in the 'functions_signal2error_PCs.R' script to compute the ERROR variance in the dataset,
#following Li and Ralph, (2019), that is using the block jackknife approach.
errorDataframe = NULL
for(i in 1:n_wind){
  
  data_snps = final_snp[[i]]

  tmp_b = computeJackknifeVarError(data_snps,jack_block=jack_block)
  errorDataframe = rbind(errorDataframe,tmp_b)
}


write.csv(file=paste0(path_out,'/',filename,'_jackBLOCK_',jack_block,'_window_size',win_size,'_signal.txt'),signalDataframe,quote = FALSE)
write.csv(file=paste0(path_out,'/',filename,'_jackBLOCK_',jack_block,'_window_size',win_size,'_error.txt'),errorDataframe,quote = FALSE)








