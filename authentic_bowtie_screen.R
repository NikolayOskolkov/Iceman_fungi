#Script for quick checking PMDtools outputs from Bowtie2 alignments for many species in metagenomic sample
#Run this script as: for i in $(cat Swab.txt); do echo ${i}; cd Swab/${i}; Rscript ancient_bowtie_screen.R; done


path<-getwd()
my_dirs<-list.dirs()[grepl("valid*",list.dirs())]
my_dirs<-gsub("./","",my_dirs)
my_taxa<-gsub("valid_","",my_dirs)

if(length(my_dirs)==0){q("no")}

CT_CpG_5_vector<-vector()
GA_CpG_3_vector<-vector()
for(i in 1:length(my_dirs))
{
if(file.exists(paste0(path,"/",my_dirs[i],"/PMD_temp.txt")) & file.info(paste0(path,"/",my_dirs[i],"/PMD_temp.txt"))$size!=0)
{
df<-read.delim(paste0(path,"/",my_dirs[i],"/PMD_temp.txt"),header=TRUE,sep="\t")
CT_CpG_5_vector<-append(CT_CpG_5_vector,df$CT_CpG_5[1])
GA_CpG_3_vector<-append(GA_CpG_3_vector,df$GA_CpG_3[1])
}
else
{
CT_CpG_5_vector<-append(CT_CpG_5_vector,NA)
GA_CpG_3_vector<-append(GA_CpG_3_vector,NA)
}
}
report<-data.frame(TaxID=my_taxa,CT_CpG_5=CT_CpG_5_vector,GA_CpG_3=GA_CpG_3_vector)
report$Status<-ifelse(report$CT_CpG_5>=0.05 & report$GA_CpG_3>=0.05,"YES","NO")
print(report)


