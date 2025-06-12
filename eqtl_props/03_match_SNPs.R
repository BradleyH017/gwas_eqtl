
## Script for obtaining control SNPs matched for MAF, LD score and gene density

args <- commandArgs(TRUE) 
group <- gsub("\\_", "\\ ", args[1]) # Do this within group:  group = "Cell type"
print(paste0("group=", group))
snplist <- "temp/output_hg19_final.tsv" #list of SNPs for which matches are found, i.e. all eQTL SNPs TESTING
info_file="snp_annotations/filter_snps.txt"
outfile=paste0("eqtl_props_out/", gsub("\\ ", "\\-", group), ".match_snps.txt")

#====

library(tidyverse)
library(data.table)
library(plyr)

#====

d_info=fread(info_file)

l2_interval=sd(d_info$L2)/10

#====

dg=fread(snplist,header=T) %>% 
  filter(annotation_type == group) %>% 
  rowwise() %>% 
  mutate(
    SNP = b37_variant
  ) %>% 
  select(SNP)


dg=dg[dg$SNP %in% d_info$SNP,]
d_data=left_join(dg,d_info,by="SNP")
d_data = d_data[complete.cases(d_data),]

# Save the per-snp results
write.table(
  d_data %>% mutate(group = group) %>% select(-c(d_tss, gene)),
  file=gsub(".match_snps.txt", ".snps.txt", outfile)
  ,quote = FALSE,sep=",",row.names = F,col.names=T)

#====
NBOOT=1000  #1000 instances of matching
x=paste0("V",as.character(1:NBOOT))
d_match=t(data.frame(x))
colnames(d_match)=x
d_match=d_match[-1,]


#iterate over SNPs
for (i in 1:nrow(d_data)){

  #focal SNP props
  SNP_temp=d_data$SNP[i]
  MAF_temp=d_data$MAF[i]
  TSSD_temp=d_data$TSSD[i]
  gene_temp=d_data$gene[i]
  l2_temp=d_data$L2[i]

  #matching
  min_maf=(MAF_temp)-0.02
  max_maf=(MAF_temp)+0.02
  if (MAF_temp<0.03){max_maf=0.03;min_maf=0.01}

  d_temp=d_info[(d_info$MAF)<=max_maf & (d_info$MAF)>=min_maf & (d_info$TSSD)==TSSD_temp,]
  d_temp$l2_diff=abs(d_temp$L2-l2_temp)
  d_temp=d_temp[d_temp$l2_diff<l2_interval,]

  d_pick=t(data.frame(d_temp[sample(1:nrow(d_temp),NBOOT,replace=T),]$SNP))
  rownames(d_pick)=i; colnames(d_pick)=x
  d_match=rbind(d_match,d_pick)
  
  print(paste0(as.character(i),": ",as.character(nrow(d_temp)),": ",as.character(length(unique(d_temp$gene)))))
 

}


#====

write.table(d_match,file=outfile,quote = FALSE,sep=",",row.names = F,col.names=F)

