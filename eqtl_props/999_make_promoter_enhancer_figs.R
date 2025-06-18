### Bradley June 2025
library(tidyverse)
library(ggplot2)
library(ggpubr)

######################
# Define options, load functions
#####################
out="eqtl_plot_out"
together = TRUE # Plot both ENCODE and FANTOM annotations together?
group_encode_enhancer = FALSE # Group both the encode enhancers together?
annot.class.palette <- c('Tissue'='#edf8b1',
                         'All Cells'='#edf8b1',
                         'Major category'='#7fcdbb',
                         'Major population'='#7fcdbb',
                         'Cell type'='#2c7fb8',
                         'All Cells-matched'="#EFF5D0",
                         'Major population-matched'='#AFD1CAF3',
                         'Cell type-matched'='#7CA3BF')

if(!file.exists(out)){
    dir.create(out)
}

calc_CI <- function(d_data,d_boot){
  

  d1= apply(d_boot,2,function(x) {unname(quantile(x,probs = seq(0,1,0.001))[26])})
  d2= apply(d_boot,2,function(x) {unname(quantile(x,probs = seq(0,1,0.001))[976])})

  d1=data.frame(d1);colnames(d1)="lower"
  d2=data.frame(d2);colnames(d2)="upper"
  d_data=data.frame(t(d_data)); colnames(d_data)="mean"
  
  d_all=cbind(d_data,d1,d2)

  return(d_all)
  
}

bar_error_plot = function(x){
    if(length(unique(x$annot))>1){
        ggplot(x, aes(x=annot, y=mean, fill=cat)) + 
            geom_col(position = position_dodge(width = 0.8), width = 0.7) + 
            theme_classic() +
            scale_fill_manual(
                values=annot.class.palette
            ) + 
            facet_grid(. ~ annot_group, scales = "free_x", space = "free_x") +
            theme(
                legend.position = "bottom",
                strip.background = element_blank()
            ) + 
            labs(
                y="Relative enrichment against promoters \n(Compared to random SNPs)",
                x="",
                fill = NULL
            ) + 
            geom_errorbar(
                aes(ymin = lower, ymax = upper),
                position = position_dodge(width = 0.8),
                width = 0.2)  
    } else {
        ggplot(x, aes(x=cat, y=mean, fill=cat)) + 
            geom_bar(stat="identity") + 
            theme_classic() +
            scale_fill_manual(
                values=annot.class.palette
            ) + 
            theme(
                legend.position = "none",
                strip.background = element_blank(),
                axis.text.x = element_text(angle=45, hjust=1)
            ) + 
            labs(
                y="Relative enrichment in enhancers over promoters \n(Compared to random SNPs)",
                x=unique(x$annot),
                fill = NULL,
                title=paste0(unique(x$annot_group), " - Enhancers")
            ) + 
            geom_errorbar(
                aes(ymin = lower, ymax = upper),
                position = position_dodge(width = 0.8),
                width = 0.2)  
    }
}

######################
# Using the bootstrapped values to calculate t-test of fold changes
#####################

# Load the enrichment results for the bootstrapped tests
inputf <- list.files("eqtl_props_out", pattern = "*(hit|actual).enrichments", full.names = TRUE)
res = do.call(rbind, lapply(inputf, function(x){
    read.delim(x, sep=",") %>% # read
        select(matches("FANTOM|ENCODE")) %>%
        select(-matches("LOEUF")) %>% 
        mutate(
            level = gsub("eqtl_props_out/promoter_enhancer.", "", x),
            level = gsub(".enrichments", "", level)
        )
}))

# Calculate the fold change of enhancers against the promoters for each dataset manually
res$FANTOM_enhancer = res$FANTOM_enhancer / res$FANTOM_promoter
res$ENCODE_pELS = res$ENCODE_pELS / res$ENCODE_PLS
res$ENCODE_dELS = res$ENCODE_dELS / res$ENCODE_PLS
res$ENCODE_enhancer_any = res$ENCODE_enhancer_any / res$ENCODE_PLS
res = res %>%
    select(level, FANTOM_enhancer, ENCODE_pELS, ENCODE_dELS, ENCODE_enhancer_any) %>% 
    rowwise() %>% 
    mutate(cat = unlist(strsplit(level, "\\."))[c(T,F)])


# Construct res dataframe
normres = do.call(rbind, lapply(unique(res$cat), function(x){
    temp = res %>% filter(cat == !!x) %>% select(-cat)
    d_eqtl = temp %>% filter(grepl("actual", level)) %>% select(-level)
    d_boot = temp %>% filter(grepl("hit", level)) %>% select(-level)
    d_all = calc_CI(d_eqtl, d_boot)
    d_all$cat = x
    return(d_all)
})) %>%
    rownames_to_column(var = "annot") %>% 
    rowwise() %>% 
    mutate(
        annot_group = ifelse(grepl("FANTOM", annot), "FANTOM", "ENCODE"), # Divide groups
        annot_group = factor(annot_group, levels = c("FANTOM", "ENCODE")),
        cat = factor(gsub("\\-", " ", cat), levels = c("All Cells", "Major population", "Cell type")),
        annot = gsub("1", "", annot),
        annot = gsub("2", "", annot), # Remove suffix from when these had to be non-unique rownames
        annot = recode(annot,
                        "FANTOM_enhancer" = "Enhancer",
                        "ENCODE_dELS" = "Enhancer-like\n(TSS distal)",
                        "ENCODE_pELS" = "Enhancer-like\n(TSS proximal)",
                        "ENCODE_enhancer_any" = "Enhancer-like"
                        ),
        annot = factor(annot, levels = c("Enhancer", "Enhancer-like\n(TSS proximal)", "Enhancer-like\n(TSS distal)", "Enhancer-like"))                                      
    )

# Decide whether to group encode
if(group_encode_enhancer){
    plotres = normres %>% 
        filter(!(annot %in% c("Enhancer-like\n(TSS proximal)", "Enhancer-like\n(TSS distal)")))
} else {
    plotres = normres %>% 
        filter(annot != "Enhancer-like")
}

# Create long res for t-test
raw_data_long <- res %>%
    filter(grepl("hit", level)) %>% 
    select(-level) %>% 
    pivot_longer(cols = colnames(res %>% select(-c(cat, level))),
                names_to = "annot", 
                values_to = "value")

# Calculate t-test
ttest_results <- raw_data_long %>%
  group_by(annot) %>%
  nest() %>%
  mutate(
    ttest = map(data, ~{
      cats <- unique(.x$cat)
      
      combinations <- combn(cats, 2, simplify = FALSE)
      
      map_dfr(combinations, function(pair) {
        group1_data <- .x$value[.x$cat == pair[1]]
        group2_data <- .x$value[.x$cat == pair[2]]
        
        test_result <- t.test(group1_data, group2_data)
        
        tibble(
          group1 = pair[1],
          group2 = pair[2],
          p_value = test_result$p.value,
          statistic = test_result$statistic
        )
      })
    })
  ) %>%
  unnest(ttest) %>%
  select(-data)

# Apply Bonferroni correction
ttest_results <- ttest_results %>%
  mutate(
    p_adj = p.adjust(p_value, method = "bonferroni"),
    significance = case_when(
      p_adj < 0.001 ~ "***",
      p_adj < 0.01 ~ "**",
      p_adj < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )


# Plot together
bar_error_plot(plotres)
ggsave(paste0(out,"/bar_enrichment_promoter_normalised.pdf"),
         width =8, height = 5,device = cairo_pdf)

# plot within group
groups = unique(plotres$annot_group)
for(g in groups){
    print(g)
    temp = plotres %>% 
        filter(annot_group == !!g)
    
    bar_error_plot(temp)
    ggsave(paste0(out,"/bar_enrichment_promoter_normalised-", g, ".pdf"),
            width =8*length(unique(temp$annot))/3, height = 5,device = cairo_pdf)
}

# Do ENCODE enhancers combined alone
print("ENCODE grouped")
temp = normres %>% filter(annot == "Enhancer-like")
bar_error_plot(temp)
ggsave(paste0(out,"/bar_enrichment_promoter_normalised-ENCODE_grouped.pdf"),
            width =8*length(unique(temp$annot))/3, height = 5,device = cairo_pdf)

temp = normres %>% filter(annot %in% c("Enhancer", "Enhancer-like"))
bar_error_plot(temp)
ggsave(paste0(out,"/bar_enrichment_promoter_normalised-ENCODE_grouped_FANTOM.pdf"),
            width =8*length(unique(temp$annot))/3, height = 5,device = cairo_pdf)

# Print t-test results
print("..T-test results")
print(as.data.frame(ttest_results))
