
## Script for obtaining bootstrapped samples of eQTLs by resampling tissues and LD blocks

args <- commandArgs(TRUE)
iter <- as.numeric(args[1]) #bootstrapping iteration; performed 1000 times: iter=1
group <- gsub("\\_", "\\ ", args[2]) # Do this within group:  group = "Cell type"
print(paste0("iter=", as.character(iter)))
print(paste0("group=", group))
eqtlfile="temp/output_hg19_final.tsv" # TESTING
blockfile="snp_annotations/snps.LD_blocks.txt"
outfile=paste0("eqtl_props_out/", gsub("\\ ", "-", group), ".boot_",as.character(iter),".txt")
if(!file.exists("eqtl_props_out")){
  dir.create("eqtl_props_out")
}


library(tidyverse)
library(data.table)
set.seed(as.numeric(iter))

print("..Loading eQTLs")
d_eqtl=fread(eqtlfile) %>% 
  filter(annotation_type == !!(group)) %>% 
  rowwise() %>% 
  mutate(
    SNP = b37_variant
  )

print("..Loading blocks")
d_blocks=fread(blockfile,header = F)
colnames(d_blocks)=c("SNP","Block")

print("..Filtering bad snps and joining")
bad_snps=d_blocks[duplicated(d_blocks$SNP),]$SNP #remove a few SNPs assigned to multiple blocks
d_blocks=d_blocks[!(d_blocks$SNP %in% bad_snps),]

df=d_eqtl[(d_eqtl$SNP %in% d_blocks$SNP),]
df=left_join(df,d_blocks,by="SNP")

#####

eqtl_blocks=unique(df$Block)
blocks_temp=sample(eqtl_blocks,size = length(eqtl_blocks),replace = T) # Sampling of LD blocks with replacement

print("..Sampling")
#bootstrap LD blocks
d_boot2=df[df$annotation_type == "XXX",]
tx=aggregate(blocks_temp, list(blocks_temp), length); colnames(tx)=c("block","count")
for (j in 1:nrow(tx)){
  dx=do.call("rbind", replicate(tx$count[j], df[df$Block==tx$block[j],], simplify = FALSE))
  d_boot2=rbind(d_boot2,dx)
}
  
#####

print("..Saving")
d_out=d_boot2 %>% select(-c(Block))
write.table(d_out,file=outfile,quote=F,sep="\t",row.names=F)
