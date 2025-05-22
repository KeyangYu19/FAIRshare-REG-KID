# df1: genes
# df2: exp variants
# df3: control variants
# df4: cCREs

gene_in <- read_csv("chddb_chd_all_genes.csv")
CM_active_cCRE <- read_table2("LVOT_v3_all_active.bed", col_names = FALSE)
CONO_v3_all_active <- read_table("CONO_v3_all_active.bed", col_names = FALSE)
CHD_kf_all_simp <- read_table("CHD_KF_dnvs.vcf", col_names = FALSE)
CHD_gelb_all_simp <- read_table("CHD_Gelb_dnvs.vcf", col_names = FALSE)
GelbStudy_hg38_dnsnps_control <- read_table("GelbStudy_hg38_dnsnps_control", col_names = FALSE)


mart<-useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl",mirror = "www")
gene_set = toupper(gene_in$X1) # gelb
gene_set = chddb_chd_all_genes$Gene #chddb
listAttributes(mart=mart)
searchAttributes(mart = mart, pattern = "hgnc")
out_genes_all <- getBM(attributes=c( "ensembl_gene_id", "hgnc_symbol", "chromosome_name", "start_position", "end_position", "external_gene_name","description"),
                   filters="hgnc_symbol", 
                   values=gene_set,
                   mart=mart, 
                   uniqueRows=TRUE)

df1 <- out_genes_all[,c(2,3,4,5)]
colnames(df1)<-c("name","chr","start","end")

df1 <- df1 %>%
  mutate(
    start_ext = start - 10000,
    end_ext = end + 10000
  )

# LVOT cCREs: CM_active_cCRE
# CONO cCREs: CONO_v3_all_active

df4 <- CM_active_cCRE # for lvot
df4 <- CONO_v3_all_active # for conotruncal

df4$chr <- gsub("chr","",df4$chr)
df4<-df4[df4$chr!="X"&df4$chr!="Y",]
df4$chr<-gsub("chr","",df4$chr)

# CHD variants(kf): kf_case_simp_subtype_match
# CHD variants(Gelb): gelb_case_simp_subtype_match
# control variants: GelbStudy_hg38_dnsnps_control

CHD_kf_all_simp <- CHD_kf_all_simp[,c(1,2,4,5)]
CHD_gelb_all_simp <- CHD_gelb_all_simp[,1:4]
                                         
df2 <- rbind(CHD_kf_all_simp,CHD_gelb_all_simp)
colnames(df2) <- c("chr","pos","ref","alt")
df2 <- df2[df2$chr != "X" & df2$chr != "Y",]
df3 <- GelbStudy_hg38_dnsnps_control[GelbStudy_hg38_dnsnps_control$chr != "X" & GelbStudy_hg38_dnsnps_control$chr != "Y",]

df2 <- df2 %>%
  semi_join(
    df1,
    by = "chr"
  ) %>%
  rowwise() %>%
  filter(any(pos >= df1$start_ext & pos <= df1$end_ext)) %>%
  ungroup()

df3 <- df3 %>%
  semi_join(
    df1,
    by = "chr"
  ) %>%
  rowwise() %>%
  filter(any(pos >= df1$start_ext & pos <= df1$end_ext)) %>%
  ungroup()

#if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
#  BiocManager::install("GenomicRanges")
#}

calculate_distance_to_RE_fast <- function(variants, REs) {
  
  # Create GRanges for variants
  variant_gr <- GRanges(
    seqnames = variants$chr,
    ranges = IRanges(start = variants$pos, end = variants$pos)
  )
  
  # Create GRanges for regulatory elements (use midpoint)
  RE_gr <- GRanges(
    seqnames = REs$chr,
    ranges = IRanges(
      start = (REs$start + REs$end) / 2,
      end = (REs$start + REs$end) / 2
    )
  )
  
  # Find nearest RE for each variant
  nearest_RE <- GenomicRanges::distanceToNearest(variant_gr, RE_gr)
  
  # Extract distance values
  distances <- rep(NA_real_, length(variant_gr))
  distances[queryHits(nearest_RE)] <- mcols(nearest_RE)$distance
  
  # Add distance back into variants
  variants$distance <- distances
  
  variants
}
df2_scored <- calculate_distance_to_RE_fast(df2, df4)
df3_scored <- calculate_distance_to_RE_fast(df3, df4)

df1 <- as_tibble(df1)

gene_test_results <- df1 %>%
  dplyr::select(name, chr, start_ext, end_ext) %>%
  mutate(
    n_exp = NA_integer_,
    n_ctrl = NA_integer_,
    ks_pval = NA_real_,
    wilcox_pval = NA_real_
  )

# Loop through each gene
for (i in seq_len(nrow(df1))) {
  
  gene_chr <- df1$chr[i]
  gene_start <- df1$start_ext[i]
  gene_end <- df1$end_ext[i]
  
  # Filter variants for gene windows
  exp_in_gene <- df2_scored %>%
    filter(chr == gene_chr, pos >= gene_start, pos <= gene_end)
  
  ctrl_in_gene <- df3_scored %>%
    filter(chr == gene_chr, pos >= gene_start, pos <= gene_end)
  
  # Save n variants
  gene_test_results$n_exp[i] <- nrow(exp_in_gene)
  gene_test_results$n_ctrl[i] <- nrow(ctrl_in_gene)
  
  # Only do tests if both groups have at least 2 points
  if (nrow(exp_in_gene) >= 2 && nrow(ctrl_in_gene) >= 2) {
    
    ks_test <- ks.test(exp_in_gene$distance, ctrl_in_gene$distance)
    wilcox_test <- wilcox.test(exp_in_gene$distance, ctrl_in_gene$distance,exact = FALSE)
    
    gene_test_results$ks_pval[i] <- ks_test$p.value
    gene_test_results$wilcox_pval[i] <- wilcox_test$p.value
    
  }
}

gene_test_results <- gene_test_results[!is.na(gene_test_results$ks_pval) & !is.na(gene_test_results$wilcox_pval),]

# Save LVOT model results
lvot_genes_ks <- gene_test_results[gene_test_results$ks_pval < 0.05,]
lvot_genes_wilcox <- gene_test_results[gene_test_results$wilcox_pval < 0.05,]
gene_test_results_lvot <- gene_test_results
write.table(gene_test_results_lvot,"gene_base_chddb_all_lvot.txt",sep = "\t",quote = F,col.names = T,row.names = F)

# Save CONO model results
cono_genes_ks <- gene_test_results[gene_test_results$ks_pval < 0.05,]
cono_genes_wilcox <- gene_test_results[gene_test_results$wilcox_pval < 0.05,]
gene_test_results_cono <- gene_test_results
write.table(gene_test_results_cono,"gene_base_chddb_all_cono.txt",sep = "\t",quote = F,col.names = T,row.names = F)

# Gene set-wide burden test for gene-based analysis results

out_genes_clean <- df1
out_genes_clean <- out_genes_clean[out_genes_clean$name %in% lvot_genes_ks$name,] # for lvot-ks genes

CONO_v3_chd_active <- CONO_v3_all_active[sapply(1:nrow(CONO_v3_all_active), function(i) { any(CONO_v3_all_active$chr[i] == out_genes_clean$chr & CONO_v3_all_active$start[i] >= out_genes_clean$start_ext & CONO_v3_all_active$end[i] <= out_genes_clean$end_ext)}), ]
CM_active_cCRE$chr <- gsub("chr","",CM_active_cCRE$chr)
CM_active_cCRE_chd <- CM_active_cCRE[sapply(1:nrow(CM_active_cCRE), function(i) { any(CM_active_cCRE$chr[i] == out_genes_clean$chr & CM_active_cCRE$start[i] >= out_genes_clean$start_ext & CM_active_cCRE$end[i] <= out_genes_clean$end_ext)}), ]

reg_nonts_ls_els_pls = reg_ls_bed[grepl("PLS",reg_ls_bed$type) | grepl("ELS",reg_ls_bed$type),] # reg_ls_bed from burden_test.r

reg_nonts_ls_els_pls <- reg_nonts_ls_els_pls[sapply(1:nrow(reg_nonts_ls_els_pls), function(i) { any(reg_nonts_ls_els_pls$chr[i] == out_genes_clean$chr & reg_nonts_ls_els_pls$start[i] >= out_genes_clean$start_ext & reg_nonts_ls_els_pls$end[i] <= out_genes_clean$end_ext)
}), ]

exp_in = CHD_kf_all_simp[,1:4]
exp_in_filtered <- exp_in %>%
  group_by(across(everything())) %>%
  filter(n() <= 5) %>%
  ungroup()

exp_all = GelbStudy_hg38_dnsnps_case[,1:4] # GelbStudy_hg38_dnsnps_case from burden_tests.r
exp_all = rbind(exp_in,exp_all) # Combine

exp_all = exp_all[order(exp_all$chr,exp_all$pos),]
exp_all = exp_all[exp_all$chr != "X" & exp_all$chr != "Y",]

control_in = GelbStudy_hg38_dnsnps_control[!GelbStudy_hg38_dnsnps_control$chr %in% c("X","Y"),1:4]

exp_in = exp_all

exp_in$gene_d = apply(exp_in,1,find_closest,out_genes_clean[,c(2,5,6)])

exp_in_range = exp_in[!(is.na(exp_in$gene_d)) & exp_in$gene_d == 0,1:4]

control_in$gene_d = apply(control_in,1,find_closest,out_genes_clean[,c(2,5,6)])

control_in_range = control_in[!(is.na(control_in$gene_d)) & control_in$gene_d == 0,1:4]

exp_in_range$cen_els_pls_lvot = apply(exp_in_range,1,find_closest,CM_active_cCRE_chd[,1:3],"mid")
control_in_range$cen_els_pls_lvot = apply(control_in_range,1,find_closest,CM_active_cCRE_chd[,1:3],"mid")

exp_in_range$cen_els_pls_cono = apply(exp_in_range,1,find_closest,CONO_v3_chd_active[,1:3],"mid")
control_in_range$cen_els_pls_cono = apply(control_in_range,1,find_closest,CONO_v3_chd_active[,1:3],"mid")

exp_in_range$cen_els_pls_nonTS = apply(exp_in_range,1,find_closest,reg_nonts_ls_els_pls[,1:3],"mid")
control_in_range$cen_els_pls_nonTS = apply(control_in_range,1,find_closest,reg_nonts_ls_els_pls[,1:3],"mid")

exp_in_range$source <- "CHD"
control_in_range$source <- "Control"

combined_df <- bind_rows(exp_in_range, control_in_range) %>%
  pivot_longer(cols = c(cen_els_pls_lvot, cen_els_pls_cono, cen_els_pls_nonTS),
               names_to = "group",
               values_to = "distance") %>%
  mutate(d_log = log10(distance+1)) 

combined_df$group <- recode(combined_df$group,
                            cen_els_pls_lvot = "LVOT",
                            cen_els_pls_cono = "CONO",
                            cen_els_pls_nonTS = "non-TS")

jpeg(file=".jpg",width=7,height=4,units="in",res = 1200)
ggplot(combined_df, aes(x = d_log, color = group, linetype = source)) +
  stat_ecdf(size = 1) +
  scale_color_brewer(palette = "Dark2") +
  theme_classic(base_size = 12) +
  theme(legend.position = "right",
        legend.title = element_blank()) +
  xlab("distance to closest cRE center") +
  ylab("fraction") +
  scale_x_continuous(breaks = c(3, 4, 5, 6),
                     labels = c("1Kbp", "10Kbp", "100Kbp", "1Mbp"),
                     limits = c(0, 6))
dev.off()

ks.test(exp_in_range$cen_els_pls_lvot,control_in_range$cen_els_pls_lvot,alternative = "greater") 
ks.test(exp_in_range$cen_els_pls_cono,control_in_range$cen_els_pls_cono,alternative = "greater") 
ks.test(exp_in_range$cen_els_pls_nonTS,control_in_range$cen_els_pls_nonTS,alternative = "greater") 

wilcox.test(x=exp_in_range$cen_els_pls_lvot, y=control_in_range$cen_els_pls_lvot, mu=0, alt="less", conf.int=T, conf.level=0.95,paired=FALSE, exact=F, correct=T)
wilcox.test(x=exp_in_range$cen_els_pls_cono, y=control_in_range$cen_els_pls_cono, mu=0, alt="less", conf.int=T, conf.level=0.95,paired=FALSE, exact=F, correct=T)
wilcox.test(x=exp_in_range$cen_els_pls_nonTS, y=control_in_range$cen_els_pls_nonTS, mu=0, alt="less", conf.int=T, conf.level=0.95,paired=FALSE, exact=F, correct=T)




