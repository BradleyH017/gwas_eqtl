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
wget https://zenodo.org/records/6618073/files/auxiliary_files.zip?download=1

# Get liftover software

# Convert the variant IDs to GRCh37 positions
```

### 3. Sampling and bootstrapping hits

For this, we are sampling variants from the LD blocks where we also find hits for each level

```         
levels=("Cell type" "Major population" "All Cells")
for level in "${levels[@]}"; do
  for i in {1..1000}; do
    bash .... 'Rscript gwas_eqtl/eqtl_props/02_bootstrap_hits.R ${i} "${level}"'
  done
done
```

### 4. Getting a set of matched snps

For each group, select the variants and identify 1000 matched variants based on their LDscore, MAF and gene density (which may vary for each set - so saving this too)

```         
levels=("Cell type" "Major population" "All Cells")
for level in "${levels[@]}"; do
  bash .... 'Rscript gwas_eqtl/eqtl_props/03_match_SNPs.R "${level}"'
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
