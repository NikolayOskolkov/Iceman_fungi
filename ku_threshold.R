
################################################ GROUND TRUTH MAIN AMETA ANALYSIS ####################################################

library("pheatmap")
gt_path<-"/home/nikolay/WABI/A_Gotherstrom/gargammel/gargammel/data_and_lists_for_main_aMeta_analysis/"
#ground_truth_matrix<-read.delim(paste0(gt_path,"ground_truth_all_microbes_mod.txt"),row.names=1,header=TRUE,sep="\t")
ground_truth_matrix<-read.delim(paste0(gt_path,"ground_truth_number_of_reads.txt"),row.names=1,header=TRUE,sep="\t")
rownames(ground_truth_matrix)<-gsub("\\.fa","",rownames(ground_truth_matrix))
ground_truth_matrix<-ground_truth_matrix[order(rownames(ground_truth_matrix)),]
colnames(ground_truth_matrix)<-paste0("Sample",seq(1,10,1))
pheatmap(ground_truth_matrix,display_numbers=TRUE,fontsize=12,cluster_cols=FALSE,cluster_rows=FALSE,
         main="Ground truth: microbial species",number_format="%i")
ground_truth_matrix[1:5,]
colSums(ground_truth_matrix)
seq_depth<-colSums(ground_truth_matrix)

ground_truth_matrix_binary<-ground_truth_matrix
ground_truth_matrix_binary[ground_truth_matrix_binary>0]<-1
ground_truth_matrix_binary[ground_truth_matrix_binary<=0]<-0
pheatmap(ground_truth_matrix_binary,display_numbers=FALSE,fontsize=12,cluster_cols=FALSE,cluster_rows=FALSE,legend_breaks = c(0,1),
         main="Binary ground truth: microbial species")
ground_truth_matrix_binary[1:5,1:5]
colSums(ground_truth_matrix_binary)


################################################# KRAKENUNIQ MAIN AMETA ANALYSIS ####################################################

ku_path<-"/home/nikolay/WABI/A_Gotherstrom/KrakenUniq/krakenuniq_aMeta_simulations_unfiltered/"
library("matrixStats")
max_kmers<-10000; step_kmers<-500
max_taxReads<-300; step_taxReads<-50
kmers_vector<-seq(from=0,to=max_kmers,by=step_kmers)
taxReads_vector<-seq(from=0,to=max_taxReads,by=step_taxReads)
sample_vector<-seq(from=1,to=10,by=1)
IoU_array<-array(rep(NA, length(sample_vector)*length(taxReads_vector)*length(kmers_vector)),
                 c(length(sample_vector), length(taxReads_vector), length(kmers_vector)))
for(s in 1:length(sample_vector))
{
print(paste0("Working with Sample",sample_vector[s]))
IoU_matrix<-matrix(NA,nrow=length(kmers_vector),ncol=length(taxReads_vector))
krakenuniq<-read.delim(paste0(ku_path,"krakenuniq.output_sample",sample_vector[s]),comment.char="#",check.names=FALSE,sep="\t")
krakenuniq$taxName<-trimws(as.character(krakenuniq$taxName))
krakenuniq<-krakenuniq[as.character(krakenuniq$rank)=="species",]
for(j in 1:length(taxReads_vector))
{
krakenuniq<-krakenuniq[krakenuniq$taxReads>taxReads_vector[j],]
#krakenuniq<-krakenuniq[krakenuniq$reads>taxReads_vector[j],]
for(i in 1:length(kmers_vector))
{
krakenuniq<-krakenuniq[krakenuniq$kmers>kmers_vector[i],]

query_list<-krakenuniq$taxName
true_list<-rownames(ground_truth_matrix_binary)[ground_truth_matrix_binary[,paste0("Sample",sample_vector[s])]==1]
length(intersect(true_list,query_list))
IoU_matrix[i,j]<-length(intersect(true_list,query_list))/length(union(true_list,query_list))

IoU_array[s,j,i]<-length(intersect(true_list,query_list))/length(union(true_list,query_list))
}
}
colnames(IoU_matrix)<-taxReads_vector
rownames(IoU_matrix)<-kmers_vector
print(head(IoU_matrix))

filled.contour(IoU_matrix,plot.axes= {
  axis(2,at=as.numeric(colnames(IoU_matrix))/max_taxReads,labels=as.numeric(colnames(IoU_matrix)))
  axis(1,at=as.numeric(rownames(IoU_matrix))/max_kmers,las=2,labels=as.numeric(rownames(IoU_matrix)))},
  nlevels=20,color.palette=terrain.colors,main=paste0("Sample",sample_vector[s]))

mtext(paste0("Optimal kmers = ",rownames(which(IoU_matrix==max(IoU_matrix), arr.ind=T))[1],", IoU_max = ",max(IoU_matrix)))


print(which(IoU_matrix==max(IoU_matrix), arr.ind=T))
}
mean_across_samples<-apply(IoU_array,c(2,3),mean)
sd_across_samples<-apply(IoU_array,c(2,3),sd)
rownames(mean_across_samples)<-taxReads_vector
colnames(mean_across_samples)<-kmers_vector
rownames(sd_across_samples)<-taxReads_vector
colnames(sd_across_samples)<-kmers_vector

mean_across_samples
max(mean_across_samples)
which(mean_across_samples==max(mean_across_samples), arr.ind=T)

sd_across_samples

mean_across_samples[1,3]
sd_across_samples[1,3]

filled.contour(mean_across_samples,plot.axes= {
  axis(2,at=as.numeric(colnames(mean_across_samples))/max_kmers,labels=as.numeric(colnames(mean_across_samples)))
  axis(1,at=as.numeric(rownames(mean_across_samples))/max_taxReads,las=2,labels=as.numeric(rownames(mean_across_samples)))},
  nlevels=30,color.palette=terrain.colors,main="Jaccard similarity averaged across samples: regular microbes dataset",
  xlab="Number of taxReads",ylab="Number of unique k-mers")

mtext(paste0("kmers_max = ",colnames(mean_across_samples)[which(mean_across_samples==max(mean_across_samples), arr.ind=T)[,"col"]],
             ", taxReads_max = ",rownames(mean_across_samples)[which(mean_across_samples==max(mean_across_samples), arr.ind=T)[,"row"]],
             ", Jaccard_max = ",max(mean_across_samples)))

optimal_kmers<-c(2000, 2000, 500, 500, 2000, 1000, 1000, 500, 1000, 1000)
plot(optimal_kmers~seq_depth,pch=19)

optimal_kmers_bins<-c(mean(c(2000, 2000, 500)), mean(c(500, 2000, 1000)), 1000, mean(c(500, 1000, 1000)))
seq_depth_bins<-c(700000, 500000, 400000, 300000)
plot(optimal_kmers_bins~seq_depth_bins,pch=19)

boxplot(c(2000, 2000, 500), c(500, 2000, 1000), 1000, c(500, 1000, 1000),names=c("7e+5","5e+5","4e+5","3e+5"),
        ylab="Optimal number of unique kmers")


#boxplot(c(2000, 2000, 500), c(500, 2000, 1000), 1000, c(500, 1000, 1000, 2000, 1000, 1000), c(1500, 6000, 6000), 
#        c(1000, 500, 1000, 500),
#        names=c("7e+5","5e+5","4e+5","3e+5","2e+5","1e+5"),
#        ylab="Optimal number of unique kmers")
boxplot(c(2000, 2000, 500), c(500, 2000, 1000), 1000, c(500, 1000, 1000, 2000, 1000, 1000), c(1000, 500, 1000, 500),
        names=c("7e+5","5e+5","4e+5","3e+5","1e+5"),
        ylab="Optimal number of unique kmers")


################################################# KRAKENUNIQ FILTER K ####################################################

ku_path<-"/home/nikolay/WABI/A_Gotherstrom/KrakenUniq/krakenuniq_aMeta_simulations_unfiltered/"
library("matrixStats")
par(mfrow=c(2,3))

max_kmers<-10000; step_kmers<-100
kmers_vector<-seq(from=0,to=max_kmers,by=step_kmers)
sample_vector<-seq(from=1,to=10,by=1)

IoU_matrix<-matrix(NA,ncol=length(sample_vector),nrow=length(kmers_vector))

for(s in 1:length(sample_vector))
{
print(paste0("Working with Sample",sample_vector[s]))

krakenuniq<-read.delim(paste0(ku_path,"krakenuniq.output_sample",sample_vector[s]),comment.char="#",check.names=FALSE,sep="\t")
krakenuniq$taxName<-trimws(as.character(krakenuniq$taxName))
krakenuniq<-krakenuniq[as.character(krakenuniq$rank)=="species",]

for(i in 1:length(kmers_vector))
{
  krakenuniq<-krakenuniq[krakenuniq$taxReads>0,]
  krakenuniq<-krakenuniq[krakenuniq$kmers>kmers_vector[i],]
      
  query_list<-krakenuniq$taxName
  true_list<-rownames(ground_truth_matrix_binary)[ground_truth_matrix_binary[,paste0("Sample",sample_vector[s])]==1]
  IoU_matrix[i,s]<-length(intersect(true_list,query_list))/length(union(true_list,query_list))
}
colnames(IoU_matrix)<-sample_vector
rownames(IoU_matrix)<-kmers_vector
print(head(IoU_matrix))
}
plot(rowMeans(IoU_matrix)~kmers_vector,type="o",xlab="K",
     ylab="Intersection Over Union with ground truth",col="darkblue",ylim=c(0,1),pch=19,main="KrakenUniq K filter")
arrows(x0=kmers_vector, y0=rowMeans(IoU_matrix)-rowSds(IoU_matrix), x1=kmers_vector, y1=rowMeans(IoU_matrix)+rowSds(IoU_matrix), 
       code=3, angle=90, length=0.03, col="darkblue")

mtext(paste0("kmers_max = ",names(rowMeans(IoU_matrix))[rowMeans(IoU_matrix)==max(rowMeans(IoU_matrix))][1],
             ", IoU_max = ",max(rowMeans(IoU_matrix))))


################################################# KRAKENUNIQ FILTER C ####################################################

ku_path<-"/home/nikolay/WABI/A_Gotherstrom/KrakenUniq/krakenuniq_aMeta_simulations_unfiltered/"
library("matrixStats")

kmers_vector<-seq(from=0,to=0.005,by=0.00001)
sample_vector<-seq(from=1,to=10,by=1)

IoU_matrix<-matrix(NA,ncol=length(sample_vector),nrow=length(kmers_vector))

for(s in 1:length(sample_vector))
{
  print(paste0("Working with Sample",sample_vector[s]))
  
  krakenuniq<-read.delim(paste0(ku_path,"krakenuniq.output_sample",sample_vector[s]),comment.char="#",check.names=FALSE,sep="\t")
  krakenuniq$taxName<-trimws(as.character(krakenuniq$taxName))
  krakenuniq<-krakenuniq[as.character(krakenuniq$rank)=="species",]
  
  for(i in 1:length(kmers_vector))
  {
    krakenuniq<-krakenuniq[krakenuniq$taxReads>0,]
    krakenuniq<-krakenuniq[krakenuniq$cov>kmers_vector[i],]
    
    query_list<-krakenuniq$taxName
    true_list<-rownames(ground_truth_matrix_binary)[ground_truth_matrix_binary[,paste0("Sample",sample_vector[s])]==1]
    IoU_matrix[i,s]<-length(intersect(true_list,query_list))/length(union(true_list,query_list))
  }
  colnames(IoU_matrix)<-sample_vector
  rownames(IoU_matrix)<-kmers_vector
  print(head(IoU_matrix))
}
plot(rowMeans(IoU_matrix)~kmers_vector,type="o",xlab="C",
     ylab="Intersection Over Union with ground truth",col="magenta",ylim=c(0,1),pch=19,main="KrakenUniq C filter")
arrows(x0=kmers_vector, y0=rowMeans(IoU_matrix)-rowSds(IoU_matrix), x1=kmers_vector, y1=rowMeans(IoU_matrix)+rowSds(IoU_matrix), 
       code=3, angle=90, length=0.03, col="magenta")

mtext(paste0("C_max = ",names(rowMeans(IoU_matrix))[rowMeans(IoU_matrix)==max(rowMeans(IoU_matrix))][1],
             ", IoU_max = ",max(rowMeans(IoU_matrix))))



################################################# KRAKENUNIQ FILTER R ####################################################

ku_path<-"/home/nikolay/WABI/A_Gotherstrom/KrakenUniq/krakenuniq_aMeta_simulations_unfiltered/"
library("matrixStats")

kmers_vector<-seq(from=0,to=1000,by=10)
sample_vector<-seq(from=1,to=10,by=1)

IoU_matrix<-matrix(NA,ncol=length(sample_vector),nrow=length(kmers_vector))

for(s in 1:length(sample_vector))
{
  print(paste0("Working with Sample",sample_vector[s]))
  
  krakenuniq<-read.delim(paste0(ku_path,"krakenuniq.output_sample",sample_vector[s]),comment.char="#",check.names=FALSE,sep="\t")
  krakenuniq$taxName<-trimws(as.character(krakenuniq$taxName))
  krakenuniq<-krakenuniq[as.character(krakenuniq$rank)=="species",]
  
  for(i in 1:length(kmers_vector))
  {
    krakenuniq<-krakenuniq[krakenuniq$taxReads>0,]
    krakenuniq<-krakenuniq[krakenuniq$taxReads>kmers_vector[i],]
    
    query_list<-krakenuniq$taxName
    true_list<-rownames(ground_truth_matrix_binary)[ground_truth_matrix_binary[,paste0("Sample",sample_vector[s])]==1]
    IoU_matrix[i,s]<-length(intersect(true_list,query_list))/length(union(true_list,query_list))
  }
  colnames(IoU_matrix)<-sample_vector
  rownames(IoU_matrix)<-kmers_vector
  print(head(IoU_matrix))
}
plot(rowMeans(IoU_matrix)~kmers_vector,type="o",xlab="R",
     ylab="Intersection Over Union with ground truth",col="cyan",ylim=c(0,1),pch=19,main="KrakenUniq R filter")
arrows(x0=kmers_vector, y0=rowMeans(IoU_matrix)-rowSds(IoU_matrix), x1=kmers_vector, y1=rowMeans(IoU_matrix)+rowSds(IoU_matrix), 
       code=3, angle=90, length=0.03, col="cyan")

mtext(paste0("R_max = ",names(rowMeans(IoU_matrix))[rowMeans(IoU_matrix)==max(rowMeans(IoU_matrix))][1],
             ", IoU_max = ",max(rowMeans(IoU_matrix))))



################################################# KRAKENUNIQ FILTER K/R ####################################################

ku_path<-"/home/nikolay/WABI/A_Gotherstrom/KrakenUniq/krakenuniq_aMeta_simulations_unfiltered/"
library("matrixStats")

ratio_vector<-seq(from=0,to=100,by=0.5)
sample_vector<-seq(from=1,to=10,by=1)

IoU_matrix<-matrix(NA,ncol=length(sample_vector),nrow=length(ratio_vector))

for(s in 1:length(sample_vector))
{
  print(paste0("Working with Sample",sample_vector[s]))
  
  krakenuniq<-read.delim(paste0(ku_path,"krakenuniq.output_sample",sample_vector[s]),comment.char="#",check.names=FALSE,sep="\t")
  krakenuniq$taxName<-trimws(as.character(krakenuniq$taxName))
  krakenuniq<-krakenuniq[as.character(krakenuniq$rank)=="species",]
  
  for(i in 1:length(ratio_vector))
  {
    krakenuniq<-krakenuniq[krakenuniq$taxReads>0,]
    krakenuniq<-krakenuniq[(krakenuniq$kmers/krakenuniq$taxReads)>ratio_vector[i],]
    
    query_list<-krakenuniq$taxName
    true_list<-rownames(ground_truth_matrix_binary)[ground_truth_matrix_binary[,paste0("Sample",sample_vector[s])]==1]
    IoU_matrix[i,s]<-length(intersect(true_list,query_list))/length(union(true_list,query_list))
  }
  colnames(IoU_matrix)<-sample_vector
  rownames(IoU_matrix)<-ratio_vector
  print(head(IoU_matrix))
}
plot(rowMeans(IoU_matrix)~ratio_vector,type="o",xlab="K / R",
     ylab="Intersection Over Union with ground truth",col="darkred",ylim=c(0,1),pch=19,main="KrakenUniq K / R filter")
arrows(x0=ratio_vector, y0=rowMeans(IoU_matrix)-rowSds(IoU_matrix), x1=ratio_vector, y1=rowMeans(IoU_matrix)+rowSds(IoU_matrix), 
       code=3, angle=90, length=0.03, col="darkred")

mtext(paste0("(K / R)_max = ",names(rowMeans(IoU_matrix))[rowMeans(IoU_matrix)==max(rowMeans(IoU_matrix))][1],
             ", IoU_max = ",max(rowMeans(IoU_matrix))))



################################################# KRAKENUNIQ FILTER (K/R)*C ####################################################

ku_path<-"/home/nikolay/WABI/A_Gotherstrom/KrakenUniq/krakenuniq_aMeta_simulations_unfiltered/"
library("matrixStats")

ratio_vector<-seq(from=0,to=0.1,by=0.0005)
sample_vector<-seq(from=1,to=10,by=1)

IoU_matrix<-matrix(NA,ncol=length(sample_vector),nrow=length(ratio_vector))

for(s in 1:length(sample_vector))
{
  print(paste0("Working with Sample",sample_vector[s]))
  
  krakenuniq<-read.delim(paste0(ku_path,"krakenuniq.output_sample",sample_vector[s]),comment.char="#",check.names=FALSE,sep="\t")
  krakenuniq$taxName<-trimws(as.character(krakenuniq$taxName))
  krakenuniq<-krakenuniq[as.character(krakenuniq$rank)=="species",]
  
  for(i in 1:length(ratio_vector))
  {
    krakenuniq<-krakenuniq[krakenuniq$taxReads>0,]
    krakenuniq<-krakenuniq[((krakenuniq$kmers/krakenuniq$taxReads)*krakenuniq$cov)>ratio_vector[i],]
    
    query_list<-krakenuniq$taxName
    true_list<-rownames(ground_truth_matrix_binary)[ground_truth_matrix_binary[,paste0("Sample",sample_vector[s])]==1]
    IoU_matrix[i,s]<-length(intersect(true_list,query_list))/length(union(true_list,query_list))
  }
  colnames(IoU_matrix)<-sample_vector
  rownames(IoU_matrix)<-ratio_vector
  print(head(IoU_matrix))
}
plot(rowMeans(IoU_matrix)~ratio_vector,type="o",xlab="(K / R) * C",
     ylab="Intersection Over Union with ground truth",col="darkgreen",ylim=c(0,1),pch=19,main="KrakenUniq (K / R) * C filter")
arrows(x0=ratio_vector, y0=rowMeans(IoU_matrix)-rowSds(IoU_matrix), x1=ratio_vector, y1=rowMeans(IoU_matrix)+rowSds(IoU_matrix), 
       code=3, angle=90, length=0.03, col="darkgreen")

mtext(paste0("((K / R) * C)_max = ",names(rowMeans(IoU_matrix))[rowMeans(IoU_matrix)==max(rowMeans(IoU_matrix))][1],
             ", IoU_max = ",max(rowMeans(IoU_matrix))))



################################################# KRAKENUNIQ FILTER (K/R)*dexp(C) ##################################################

ku_path<-"/home/nikolay/WABI/A_Gotherstrom/KrakenUniq/krakenuniq_aMeta_simulations_unfiltered/"
library("matrixStats")

ratio_vector<-seq(from=0,to=100,by=1)
sample_vector<-seq(from=1,to=10,by=1)

IoU_matrix<-matrix(NA,ncol=length(sample_vector),nrow=length(ratio_vector))

for(s in 1:length(sample_vector))
{
  print(paste0("Working with Sample",sample_vector[s]))
  
  krakenuniq<-read.delim(paste0(ku_path,"krakenuniq.output_sample",sample_vector[s]),comment.char="#",check.names=FALSE,sep="\t")
  krakenuniq$taxName<-trimws(as.character(krakenuniq$taxName))
  krakenuniq<-krakenuniq[as.character(krakenuniq$rank)=="species",]
  
  for(i in 1:length(ratio_vector))
  {
    krakenuniq<-krakenuniq[krakenuniq$taxReads>0,]
    krakenuniq<-krakenuniq[((krakenuniq$kmers/krakenuniq$taxReads)*(1.3^(18^krakenuniq$cov)))>ratio_vector[i],]
    
    query_list<-krakenuniq$taxName
    true_list<-rownames(ground_truth_matrix_binary)[ground_truth_matrix_binary[,paste0("Sample",sample_vector[s])]==1]
    IoU_matrix[i,s]<-length(intersect(true_list,query_list))/length(union(true_list,query_list))
  }
  colnames(IoU_matrix)<-sample_vector
  rownames(IoU_matrix)<-ratio_vector
  print(head(IoU_matrix))
}
plot(rowMeans(IoU_matrix)~ratio_vector,type="o",xlab="(K / R) * dexp(C)",
     ylab="Intersection Over Union with ground truth",col=j,ylim=c(0,1),pch=19,main="KrakenUniq (K / R) * dexp(C) filter")
arrows(x0=ratio_vector, y0=rowMeans(IoU_matrix)-rowSds(IoU_matrix), x1=ratio_vector, y1=rowMeans(IoU_matrix)+rowSds(IoU_matrix), 
       code=3, angle=90, length=0.03, col=j)

mtext(paste0("((K / R) * dexp(C))_max = ",names(rowMeans(IoU_matrix))[rowMeans(IoU_matrix)==max(rowMeans(IoU_matrix))][1],
             ", IoU_max = ",max(rowMeans(IoU_matrix))))


################################################ CORRELATION BETWEEN QUALITY FILTERS #################################################

ku_path<-"/home/nikolay/WABI/A_Gotherstrom/KrakenUniq/krakenuniq_aMeta_simulations_unfiltered/"
sample_vector<-seq(from=1,to=10,by=1)
for(s in 1:3)#length(sample_vector)
{
  print(paste0("Working with Sample",sample_vector[s]))
  
  krakenuniq<-read.delim(paste0(ku_path,"krakenuniq.output_sample",sample_vector[s]),comment.char="#",check.names=FALSE,sep="\t")
  krakenuniq$taxName<-trimws(as.character(krakenuniq$taxName))
  krakenuniq<-krakenuniq[as.character(krakenuniq$rank)=="species",]
  krakenuniq<-krakenuniq[krakenuniq$taxReads>0,]
}
pheatmap(cor(subset(krakenuniq,select=c("%","reads","taxReads","kmers","cov","dup")),method="spearman"),display_numbers=TRUE)
