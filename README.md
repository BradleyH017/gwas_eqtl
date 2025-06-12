# Promoter/enhancer enrichment to match that of Mostafavi et al., 2023

### Bradley May 2025

## Code

Code was obtained from author's github (<https://github.com/hakha-most/gwas_eqtl.git>)

## Data

snp annotations were obtained from the author's zenodo `wget https://zenodo.org/records/6618073/files/snp_annotations.zip?download=1 -O snp_annotations.zip`

## Overview of approach

1.  Annotate each eQTL by the minimum resolution at which an effect has been obtained (done during plotting scripts)
2.  Convert these positions --\> GRCh37
3.  Sample 1000 hits per resolution, using the LD blocks previously defined (Berisa and Pickrell 2016), sample 1000 other variants, per hit (bootstapped hits)
4.  For each lead variant, calculate 1000 matched snps using the author's code (control)
5.  Quantify the overlap of bootstrapped hits and control snps with the FANTOM and ENCODE genomic annotations

### 1. Annotating eQTLs by min resolution

This was carried out during plotting of eQTLs. See code below:

```         
save = cond_all # Conditional eQTLs mapped to phenotype clump indexes

df.count.sig.egenes <- save %>% 
  mutate(Gene_class = factor(case_when(biotype == 'protein_coding' ~ 'Protein coding',
                                       biotype == 'lncRNA' ~ 'lncRNA',
                                       TRUE ~ 'Other'),
                             levels = c('Other', 'lncRNA', 'Protein coding'))) %>% 
  arrange(Level) %>% 
  distinct(!!sym(distinction), .keep_all = TRUE) 

level_seperation = df.count.sig.egenes %>% distinct(!!sym(distinction), annotation, annotation_type)
  level_seperation = level_seperation %>% 
    rowwise() %>% 
    mutate(index = unlist(strsplit(phenotype_clump_index, "\\-"))[c(F,T)]) %>% 
    select(phenotype_clump_index, annotation_type) %>% 
    filter(!is.na(index))

write.table(level_seperation, paste0("eqtl_out/eqtl_seperation_for_promoter_enhancer_enrichment_.txt"), row.names=F, quote=F, sep = "\t")
```

### 2. Converting eQTLs to GRCh37

The same chain file as used by the author's was downloaded from <https://zenodo.org/records/6618073/files/auxiliary_files.zip?download=1>

```         
# Get the chain file
wget https://zenodo.org/records/6618073/files/auxiliary_files.zip?download=1 -O auxiliary_files.zip
unzip auxiliary_files.zip
rm auxiliary_files.zip

# Convert the input file to bed format (keep the annotation_type column as this is required for division of leads in the subsequent scripts)
awk -F'\t' 'NR > 1 {
    split($2, v, ":")
    chrom = v[1]
    pos = v[2] - 1
    end = v[2]
    ref = v[3]
    alt = v[4]
    print chrom"\t"pos"\t"end"\t"ref"/"alt"\t0\t.\t1\t"$3"\t"$4
}' temp/gene_level_seperation.txt > temp/temp_output.bed
# Assuming 3rd and 4th columns are 'annotation_type' and 'phenotype_clump_index'

# Convert using liftover
lo=/software/team152/bh18/liftOver
$lo -bedPlus=7 temp/temp_output.bed auxiliary_files/hg38ToHg19.over.chain.gz temp/output_hg19.bed temp/unmapped.bed

# Add new variant id col and simplify
awk -F'\t' '{
    sub(/^chr/, "", $1)           # Remove "chr" prefix
    gsub("/", ":", $4)            # Replace "/" with ":"
    b37_variant = $1":"$3":"$4
    print b37_variant"\t"$8" "$9"\t"$10
}' temp/output_hg19.bed > temp/output_hg19_clean.tsv
echo -e "b37_variant\tannotation_type\tphenotype_clump_index" | cat - temp/output_hg19_clean.tsv > temp/output_hg19_final.tsv
# Final format is b37_variant, annotation_type, phenotype_clump_index
```

### 3. Sampling and bootstrapping hits

For this, we are sampling variants from the LD blocks where we also find hits for each level

```  
mkdir -p logs       
MEM=3000
levels=("Cell_type" "Major_population" "All_Cells")
for level in "${levels[@]}"; do
  for i in {1..1000}; do
    bsub -J "sample_hits-${level}-${i}" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G team152 \
      -e logs/sample_hits-${level}-${i}-stderr \
      -o logs/sample_hits-${level}-${i}-stdout \
      "Rscript eqtl_props/02_bootstrap_hits.R $i $level> \
      logs/process_sample_hits-${level}-${i}.Rout"
  done
done
```

### 4. Getting a set of matched snps

For each group, select the variants and identify 1000 matched variants based on their LDscore, MAF and gene density (which may vary for each set - so saving this too)

```         
levels=("Cell_type" "Major_population" "All_Cells")
for level in "${levels[@]}"; do
   bsub -J "get_matched-${level}" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G team152 \
      -e logs/get_matched-${level}-stderr \
      -o logs/get_matched-${level}-stdout \
      "Rscript eqtl_props/03_match_SNPs.R $level> \
      logs/get_matched-${level}.Rout"
done
```

### 5. Calculating enrichment scores

For each group, and each of hits and matched snps, calculate the enrichment score

```         
levels=("Cell type" "Major population" "All Cells")
hit_or_match=("hit" "match")
for level in "${levels[@]}"; do
  for hm in "${hit_or_match[@]}"; do
    bash .... 'Rscript gwas_eqtl/eqtl_props/06_calc_promoter_enhancer_props.R "${level}" "${hm}"'
  done
done
```

### 6. Calculating confidence intervals and results

Gather the results and compute confidence intervals

```         
levels=("Cell type" "Major population" "All Cells")
for level in "${levels[@]}"; do
  bash .... 'Rscript gwas_eqtl/eqtl_props/06_extract_promoter_enhancer_enrichments_CI.R "${level}"'
done
```

##### TEMP

Get a list of test variants for the analysis by sampling the snp

```         
bash getting_test_variants.sh #Already in GRCh37
```
