#This script is used to plot the difference between SIGNAL and ERROR variances across various window sizes.
#In our work, we decided to select the minimum window size reaching a plateau value in the difference between
#SIGNAL and ERROR variances.

require(ggplot2)

path_input = args[1]
pre_filename = args[2]
first_size = as.numeric(args[3])
int_size = as.numeric(args[4])
last_size = as.numeric(args[5])
jack_block = as.numeric(args[6])

#It lists the tested window sizes.
all_sizes = seq(first_size,last_size,int_size)

final_tot = NULL

#It loops over all window sizes.
for(win_size in all_sizes){
  
  #It reads the ERROR dataframe.
  tmp_error = read.csv(paste0(path_input,'/',dataset,'_jackBLOCK_',jack_block,'_window_size',win_size,'_error.txt'),header = TRUE)[,2]
  
  #It reads the SIGNAL dataframe.
  tmp_signal = read.csv(paste0(path_input,'/',dataset,'_jackBLOCK_',jack_block,'_window_size',win_size,'_signal.txt'),header = TRUE)[,2]
  
  #It creates a dataframe with the results of each window size.
  df_tot = data.frame(win_size=win_size,signal=mean(tmp_signal,na.rm = T),error=mean(tmp_error,na.rm = T))
  
  final_tot =  rbind(final_tot,df_tot)
  rm(tmp_error,tmp_signal,df_tot)
}

final_tot$delta = final_tot$signal - final_tot$error

deltaPLOT = ggplot(final_tot,aes(x=win_size,y=delta))+
  geom_point()+
  geom_smooth()+
  xlab('Window size (n° of SNPs)')+
  ylab('Signal - error')+
  theme_bw()


out_file = paste0(path_input,'/',pre_filename)

ggsave(file=paste0(out_file,"_deltaSignal2Error_new.png"),plot=deltaPLOT,device="png",width=4,height=4)


