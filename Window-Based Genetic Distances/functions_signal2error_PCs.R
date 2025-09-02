

####Functions to compute variance SIGNAL####
varfunction<-function(data,nval){
  var=sum((data-mean(data))^2)/(nval)
  return(var)
}

computeVarSignal<-function(data_df){
  
  #It creates a matrix of the same size as 'data_df'.
  A_final=matrix(0,nrow=nrow(data_df),ncol=ncol(data_df))
  #It creates a vector of the size of the number of genomic windows included in 'data_df'.
  S_vec=rep(0,nrow(A_final))
  
  #It loops over each genomic window.
  for(i in 1:nrow(A_final)){
    #It chooses the sign that best match the relative positioning of each individual in PC1 
    # across genomic windows. The 'if conditioning' considers that if ||PC1 - PC1j|| > |PC1 + PCj||
    #then PC1j = -PC1j. 'PC1j' is equivalent to 'u' in Li and Ralph, (2019). 
    if (sum((data_df[1,]-data_df[i,])^2)<sum((data_df[1,]+data_df[i,])^2))
    {
      S_vec[i]=1
      A_final[i,]=data_df[i,]	
    }
    else 
    {
      S_vec[i]=-1
      A_final[i,]=-data_df[i,]	
    }
  }
  #It computes the variance SIGNAL#
 b_final=apply(A_final,2,varfunction,nval=nrow(A_final))
  
  return(b_final)
}
####Functions to compute variance SIGNAL####

####Functions to compute variance ERROR####
#It computes the squared standard error#
varfunctionJackknifeError<-function(data,jack_block){
  var=(jack_block-1)/jack_block*sum((data-mean(data))^2)
  return(var)
}

#This function performs the block jackknife to estimate the standard error
#for each PC1ij where, 'i' indicates a specific individual and 'j' the genomic window under consideration.
get.eigenvector <- function(x, d,jack_block) {
  step <- round(nrow(d)/jack_block)
  chunk <- d[-(((x-1)*step + 1):(x*step)), ]
  temp<-chunk
  temp<-data.matrix(temp)
  data=temp
  M=rowMeans(data,na.rm=TRUE)
  M=rep(M,times=ncol(data))
  M=matrix(M,nrow=nrow(data),ncol=ncol(data),byrow=FALSE)
  data=data-M
  cov=cov(data,use="pairwise")
  if(sum(is.na(cov))>0) {return(rep(NA,nrow(cov)))}
  PCA=eigen(cov)
  Vec=PCA$vectors
  lam=PCA$values
  PC1=Vec[,1]
  return(PC1)
}

computeJackknifeVarError<-function(df,jack_block=10){ 
  
  PC1s <- sapply(1:jack_block, get.eigenvector, d=df,jack_block=jack_block)
  
  if(sum(is.na(PC1s))>0) {return(NA)}
  
  a=as.matrix(PC1s)
  a=t(a)
  A=matrix(0,nrow=nrow(a),ncol=ncol(a))
  S=rep(0,nrow(A))
  for(i in 1:nrow(A))
  {
    if (sum((a[1,]-a[i,])^2)<sum((a[1,]+a[i,])^2))
    {
      S[i]=1
      A[i,]=a[i,]     
    }
    else 
    {
      S[i]=-1
      A[i,]=-a[i,]    
    }
  }
  b=apply(A,2,varfunctionJackknifeError,jack_block=jack_block)
  return(mean(b))
}
####Functions to compute variance ERROR####
