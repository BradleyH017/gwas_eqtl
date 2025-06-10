
## Script to compute promoter-enhancer enrichments of eQTL SNPs.

# This script, as provided, is run on eQTL SNPs. 
# The same script was run using 1000 bootstrapped set of SNPs, as well as 1000 instances of matched SNPs.

set.seed(0)

library(data.table)
library(tidyverse)
library(dplyr)

args <- commandArgs(TRUE) 
group <- args[1] # group = "Cell type"
hit_or_match = args[2] # hit_or_match="hit"
eqtlfile="temp/sampled_ready.txt" # TESTING
print(paste0("hit_or_match=", as.character(hit_or_match)))
print(paste0("group=", group))
infofile="snp_annotations/promoter_enhancer.counts"
outfile=paste0("eqtl_props_out/promoter_enhancer.",gsub("\\ ", "\\-", group), ".", hit_or_match, ".enrichments")

d_info=fread(infofile)

# Construct eqtl file from hits or matches
if(hit_or_match == "hit"){
  hitfs = list.files("eqtl_props_out", pattern = paste0(gsub("\\ ", "\\-", group), ".boot*"), full=TRUE)
  d_eqtl_list = lapply(hitfs, function(x){
    read.delim(x) %>% 
      select(SNP)
  })
} else {
  hitf = paste0("eqtl_props_out/", gsub("\\ ", "\\-", group), ".match_snps.txt")
  d_match = read.csv(hitf, header=F)
  d_eqtl_list = vector("list", length=1000)
  for(i in 1:1000){
    d_eqtl_list[[i]] = data.frame("SNP" = d_match[,i])
  }
}

enrichments = vector("list", length = length(d_eqtl_list))
p0=sapply((d_info %>% select(-SNP)), function(x) length(x[x>0])/length(x) ) # background enrichment of all SNPs in promoters/enhancers
for(i in 1:length(d_eqtl_list)){
  print(i)
  d_eqtl_annots=left_join(d_eqtl_list[[i]],d_info,by="SNP")
  p1=sapply(d_eqtl_annots %>% select(-SNP), function(x) length(x[x>0])/length(x) ) 
  enrichments[[i]]=data.frame(t(p1/p0))
}
enrichments = do.call(rbind, enrichments)

write.table(enrichments,file=outfile,quote=F,sep=",",row.names=F)

# Also compute this for the actual enrichment
if(hit_or_match == "hit"){
 d_eqtl=fread(eqtlfile) %>% 
   filter(annotation_type == !!(group)) %>% 
   rowwise() %>% 
   mutate(
     SNP = unlist(strsplit(phenotype_clump_index, "\\-"))[c(F,T)]
   )
 d_eqtl_annots=left_join(d_eqtl %>% select(SNP),d_info,by="SNP")
 d_eqtl_annots = d_eqtl_annots[complete.cases(d_eqtl_annots), ]
 p1=sapply(d_eqtl_annots %>% select(-SNP), function(x) length(x[x>0])/length(x) ) # Actual hit mean
 enrichments=data.frame(t(p1/p0))
 actual_out=gsub("hit", "actual", outfile)
 write.table(enrichments,file=actual_out,quote=F,sep=",",row.names=F)
}
