### Bradley June 2025
library(tidyverse)
library(ggplot2)
out="eqtl_plot_out"
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

# Load in the files
inputf = list.files("eqtl_props_out", pattern = ".enrichment.CI", full=TRUE)
res = do.call(rbind, lapply(inputf, function(x){
    read.delim(x, sep=",") %>% 
        mutate(
            hit_or_match = ifelse(grepl(".matched", cat), "match", "hit"),
            level = gsub(".matched", "", cat))
}))

# Calculate enrichments of hits over matched for each annotation and each level
res = res %>%
    filter(annot %in% c("ENCODE_PLS", "ENCODE_dELS", "ENCODE_pELS", "FANTOM_enhancer", "FANTOM_promoter")) %>% 
    rowwise() %>%
    mutate(
        annot_group = unlist(strsplit(annot, "\\_"))[c(T,F)],
        annot_group = factor(annot_group, levels = c("FANTOM", "ENCODE")),
        level = gsub("\\-", "\\ ", level),
        level = factor(level, levels = c("All Cells", "Major population", "Cell type")),
        annot = recode(annot,
            "FANTOM_promoter" = "Promoter",
            "FANTOM_enhancer" = "Enhancer",
            "ENCODE_PLS" = "Promoter-like",
            "ENCODE_pELS" = "Enhancer-like\n(TSS proximal)",
            "ENCODE_dELS" = "Enhancer-like\n(TSS distal)"
        ),
        annot = factor(annot, levels = c("Promoter", "Enhancer", "Promoter-like", "Enhancer-like\n(TSS proximal)", "Enhancer-like\n(TSS distal)")),
        cat = gsub("\\-", "\\ ", cat),
        cat = gsub("\\.", "\\-", cat),
        cat = factor(cat, levels = c("All Cells", "All Cells-matched", "Major population", "Major population-matched", "Cell type", "Cell type-matched"))
    )

ggplot(res, aes(x=annot, y=mean, fill=cat)) + 
    geom_col(position = position_dodge(width = 0.8), width = 0.7) + 
    theme_classic() +
    scale_fill_manual(
        values=annot.class.palette
    ) + 
    geom_hline(yintercept = 1, linetype = "dashed", color = "black") + 
    facet_grid(. ~ annot_group, scales = "free_x", space = "free_x") +
    theme(
        legend.position = "bottom",
        strip.background = element_blank()
    ) + 
    labs(
        y="Enrichment\n(Compared to random SNPs)",
        x="",
        fill = NULL
    ) +
    geom_errorbar(
        aes(ymin = lower, ymax = upper),
        position = position_dodge(width = 0.8),
        width = 0.2
    )

ggsave(paste0(out,"/bar_enrichment.pdf"),
         width =8, height = 5,device = cairo_pdf)


# Express this as enhancers relative to promoters
encode_ref = res %>% filter(annot_group == "ENCODE" & annot == "Promoter-like") %>% rename(refcat = cat)
fantom_ref = res %>% filter(annot_group == "FANTOM" & annot == "Promoter") %>% rename(refcat = cat)

normres <- res %>%
    group_by(annot_group, cat) %>%
    rowwise() %>% 
    mutate(
        reference_mean = ifelse(annot_group == "ENCODE", encode_ref[encode_ref$refcat == cat,]$mean, fantom_ref[fantom_ref$refcat == cat,]$mean),
        mean_normalised = mean / reference_mean
    )
  
ggplot(normres %>% filter(!grepl("Promoter", annot) & hit_or_match == "hit"), aes(x=annot, y=mean_normalised, fill=cat)) + 
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
        y="Relative enrichment in enhancers vs promoters \n(Compared to random SNPs)",
        x="",
        fill = NULL
    )

ggsave(paste0(out,"/bar_enrichment_promoter_normalised.pdf"),
         width =8, height = 5,device = cairo_pdf)

# Divide by annot group
for(g in unique(normres$annot_group)){
    dat = normres %>% filter(!grepl("Promoter", annot) & hit_or_match == "hit", annot_group == !!(g))
    print(g)
    if(g == "ENCODE"){
        ggplot(dat, aes(x=annot, y=mean_normalised, fill=cat)) + 
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
            y="Relative enrichment in enhancers vs promoters \n(Compared to random SNPs)",
            x="",
            fill = NULL
        )
    } else {
        ggplot(dat, aes(x=cat, y=mean_normalised, fill=cat)) + 
        geom_bar(stat="identity") + 
        theme_classic() +
        scale_fill_manual(
            values=annot.class.palette
        ) + 
        theme(
            legend.position = "none"
        ) + 
        labs(
            y="Relative enrichment in enhancers vs promoters \n(Compared to random SNPs)",
            x="",
            fill = NULL
        ) + 
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size=12))
    }
    ggsave(paste0(out,"/bar_enrichment_promoter_normalised-", g, ".pdf"),
            width =4*length(unique(dat$annot)), height = 5,device = cairo_pdf)   
}
