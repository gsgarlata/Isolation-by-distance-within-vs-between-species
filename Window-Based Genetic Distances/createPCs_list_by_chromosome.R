
require(lostruct)

path_data = args[1]
chrom_bcf = args[2]
wind_size = as.numeric(args[3])
outgroup_list = args[4]
chrom_ID = args[5]
filename = args[6]

file_chrm = paste0(path_data,'/',chrom_bcf)

snps_chrm <- vcf_windower(file_chrm,size=wind_size,type='snp')
pcs_chrm <- eigen_windows(snps_chrm,k=2)

pc1_chrm_coord = grep('PC_1_',colnames(pcs_chrm))
pc1_final = pcs_chrm[,pc1_chrm_coord]
pc1_final = pc1_final[,-c(1:3)]

if(is.null(outgroup_list)){}else{
  pc1_final = pc1_final[,-match(outgroup_list,colnames(pc1_final))]
}

write.csv(file=paste0(path_out,'/',filename,'_chromo',chrom_ID,'_window_size',win_size,'_pcs.txt'),pc1_final,quote = FALSE)

