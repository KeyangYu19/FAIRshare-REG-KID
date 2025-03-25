# libraries
library(readr)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(httr)
library(jsonlite)
library(tidyverse)
require(biomaRt)

# functions
read_f <- function(filename){ 
  data = read_delim(filename, delim="\t", escape_double = FALSE,  col_names = FALSE, comment = "#", trim_ws = TRUE)
  data = data[,1:4]
  colnames(data) = c("chr","pos","ref","alt")
  data = data |> filter(ref %in% c("A","C","T","G") & alt %in% c("A","C","T","G"))
  data$chr = gsub("chr","", data$chr)
  return(data)
}

read_re_ls <- function(re_ls_f){
  re_ls = read_csv(re_ls_f, col_names = FALSE)
  re_ls_name = re_ls$X1
  return(re_ls_name)
}

find_closest <- function(x, in_bed, mode = "end") {
  chr_v <- as.character(x[1])
  end_v <- as.numeric(x[2])
  
  if (!(chr_v %in% in_bed$chr)) {
    return(NA)
  } 
  
  subset <- in_bed[in_bed$chr == chr_v, ]
  ent_distances <- apply(subset, 1, get_distance, coord = end_v, mode = mode)
  
  return(min(ent_distances))
}


get_distance <- function(in_bed, coord, mode = "end") {
  bed_start <- as.numeric(in_bed[2])
  bed_end <- as.numeric(in_bed[3])
  
  if (mode == "end") {
    if (coord <= bed_start || coord >= bed_end) {
      return(abs(bed_start - coord))
    } else {
      return(0)
    }
  }
  
  if (mode == "mid") {
    bed_mid <- (bed_start + bed_end) / 2
    return(abs(bed_mid - coord))
  }
}

# data import

GelbStudy_hg38_dnsnps_case = read_f("GelbStudy_hg38_dnvs_case.vcf") # case dnvs from Gelb paper
GelbStudy_hg38_dnsnps_control = read_f("GelbStudy_hg38_dnvs_control.vcf") # control/healthy dnvs from Gelb paper

name = "Priority_HLHS_gene_rank_1_cat_1"
ldhf_dir = "HLHS_results_v2/Priority_HLHS_gene_rank_1_cat_1/"

# SCREEN cCRE coordinates
GRCh38_cCREs <- read_table2("GRCh38-cCREs.bed", col_names = FALSE)

colnames(GRCh38_cCREs) = c("chr","start","end","sample","id","type")
GRCh38_cCREs = GRCh38_cCREs[,-4]
GRCh38_cCREs = GRCh38_cCREs %>% separate(type, c("type","CTCF"),sep = ",",fill = "right")
GRCh38_cCREs$CTCF[is.na(GRCh38_cCREs$CTCF)] <- "."

re_ls_exp = read_re_ls(paste(ldhf_dir,name,".RE.all.uniq.clean.txt",sep=""))
re_ts_ls = read_csv(paste(ldhf_dir,name,".RE.active.TS.uniq.txt",sep=""), col_names = FALSE)

# Filter for LVOT specific cREs using detailed parsed API outputs
reg_ts_ls2 <- read_table2("HLHS_results_v2/Priority_HLHS_gene_rank_1_cat_1/Priority_HLHS_gene_rank_1_cat_1.RE.ENCODE.active.TS.txt", col_names = FALSE)
reg_ts_ls2 = reg_ts_ls2[grepl("left_ventricle",reg_ts_ls2$X3) | grepl("cardiac_muscle_cell",reg_ts_ls2$X3),]
reg_ts_ls2_els = reg_ts_ls2[grepl("ELS",reg_ts_ls2$X2),]
reg_ts_ls2_els = unique(reg_ts_ls2_els$X1)
reg_ts_ls2_pls = reg_ts_ls2[grepl("PLS",reg_ts_ls2$X2),]
reg_ts_ls2_pls = unique(reg_ts_ls2_pls$X1)

# Genes inputs
input_genes <- read_csv("HLHS_results_v2/Priority_HLHS_gene_rank_1_cat_1.csv", col_names = FALSE)

# Get cRE positions
reg_ls_bed = GRCh38_cCREs
reg_ls_bed = reg_ls_bed[reg_ls_bed$chr != "chrX" & reg_ls_bed$chr != "chrY",]
reg_ls_bed$center = (reg_ls_bed$start+reg_ls_bed$end)/2
reg_ts_ls_bed = reg_ls_bed[reg_ls_bed$id %in% re_ts_ls$X1,]
#reg_ts_ls_bed = reg_ts_ls_bed[grepl("ELS",reg_ts_ls_bed$type)|grepl("PLS",reg_ts_ls_bed$type),]
reg_ts_ls_bed_el = reg_ts_ls_bed[grepl("ELS",reg_ts_ls_bed$type),]
reg_ts_ls_bed_pl = reg_ts_ls_bed[grepl("PLS",reg_ts_ls_bed$type),]

reg_ts_ls_bed_el$chr = gsub("chr","",reg_ts_ls_bed_el$chr)
reg_ts_ls_bed_pl$chr = gsub("chr","",reg_ts_ls_bed_pl$chr)

reg_ts_ls_bed_el2 = reg_ts_ls_bed_el[reg_ts_ls_bed_el$id %in% reg_ts_ls2_els,]
reg_ts_ls_bed_pl2 = reg_ts_ls_bed_pl[reg_ts_ls_bed_pl$id %in% reg_ts_ls2_pls,]


# Get gene information
# Biomart API
mart<-useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl",mirror = "useast")
gene_set = toupper(input_genes$X1)
listAttributes(mart=mart)
searchAttributes(mart = mart, pattern = "hgnc")
out_genes <- getBM(attributes=c( "ensembl_gene_id", "hgnc_symbol", "chromosome_name", "start_position", "end_position", "external_gene_name","description"),
                   filters="hgnc_symbol", 
                   values=gene_set,
                   mart=mart, 
                   uniqueRows=TRUE)
out_genes_clean = out_genes[,2:5]
out_genes_clean$re_qstart_close = out_genes_clean$start_position - 10000
out_genes_clean$re_qend_close = out_genes_clean$end_position + 10000
#out_genes_clean$re_qstart_far = out_genes_clean$start_position - 100000
#out_genes_clean$re_qend_far = out_genes_clean$end_position + 100000
colnames(out_genes_clean)[2] = "chr"
out_genes_clean = out_genes_clean[order(out_genes_clean$chr,out_genes_clean$start_position),]


# Variant cleanup
control_in = GelbStudy_hg38_dnsnps_control[GelbStudy_hg38_dnsnps_control$chr != "X",1:4]
exp_in = exp_all[GelbStudy_hg38_dnsnps_case$chr != "X" & GelbStudy_hg38_dnsnps_case$chr != "Y",1:4]

exp_in$gene_d = apply(exp_in,1,find_closest,out_genes_clean[,c(2,5,6)])
exp_in_range = exp_in[!(is.na(exp_in$gene_d)) & exp_in$gene_d == 0,1:4]

control_in$gene_d = apply(control_in,1,find_closest,out_genes_clean[,c(2,5,6)])
control_in_range = control_in[!(is.na(control_in$gene_d)) & control_in$gene_d == 0,1:4]

### Burden test for LVOT model
reg_ts_ls_els2_pls2 = rbind(reg_ts_ls_bed_el2,reg_ts_ls_bed_pl2)

exp_in_range$cen_els_pls = apply(exp_in_range,1,find_closest,reg_ts_ls_els2_pls2[,1:3],"mid")
control_in_range$cen_els_pls = apply(control_in_range,1,find_closest,reg_ts_ls_els2_pls2[,1:3],"mid")

# Rank sum test
wilcox.test(x=exp_in_range$cen_els_pls, y=control_in_range$cen_els_pls, mu=0, alt="less", conf.int=T, conf.level=0.95,paired=FALSE, exact=F, correct=T)

# Fisher's exact
count_1000_els_pls_paper <-data.frame(
  "in_range" = c(8,4),
  "out" = c(60,122),
  row.names = c("chd","control"),
  stringsAsFactors = F
)
fx_test_1000_paper <- fisher.test(count_1000_els_pls_paper)

# plot
ks_plot_df_exp2 = exp_in_range[,5]
ks_plot_df_exp2$type = "case"

ks_plot_df_control2 = control_in_range[,5]
ks_plot_df_control2$type = "control"

ks_plot_df2 = rbind(ks_plot_df_exp2,ks_plot_df_control2)
ks_plot_df2$d_log = log10(ks_plot_df2$cen_els_pls)

jpeg(file="allGreb_screen_v3.jpg",width=7,height=4,units="in",res = 1200)
ggplot(ks_plot_df2, aes(x = d_log, group = type, color = type))+
    stat_ecdf(size=1) +
    scale_color_brewer(palette="Paired") +
    theme_classic(base_size = 12) +
    theme(legend.position ="right") +
    xlab("distance to closest cRE center") +
    ylab("fraction") +
    scale_x_continuous(breaks=c(3,4,5),labels= c("1Kbp","10Kbp","100Kbp"),limits = c(1, 5.1)) +
    ylab("fraction") +
    theme(legend.title=element_blank())
dev.off()

### Burden test for nonTS model
reg_ls_bed$chr = gsub("chr","",reg_ls_bed$chr)

reg_nonts_ls_els_pls = reg_ls_bed[grepl("PLS",reg_ls_bed$type) | grepl("ELS",reg_ls_bed$type),]

# Filter for in-range (10kbps) cREs
reg_nonts_ls_els_pls <- reg_nonts_ls_els_pls[sapply(1:nrow(reg_nonts_ls_els_pls), function(i) {
  any(reg_nonts_ls_els_pls$chr[i] == out_genes_clean$chr & reg_nonts_ls_els_pls$start[i] >= out_genes_clean$re_qstart_close & reg_nonts_ls_els_pls$end[i] <= out_genes_clean$re_qend_close)
}), ]

exp_in_range$cen_els_pls_nonTS = apply(exp_in_range,1,find_closest,reg_nonts_ls_els_pls[,1:3],"mid")
control_in_range$cen_els_pls_nonTS = apply(control_in_range,1,find_closest,reg_nonts_ls_els_pls[,1:3],"mid")

wilcox.test(x=exp_in_range$cen_els_pls_nonTS, y=control_in_range$cen_els_pls_nonTS, mu=0, alt="less", conf.int=T, conf.level=0.95,paired=FALSE, exact=T, correct=T)

enrichment_nonTS_paper = 29/68/(56/126)
fisher_nonts_paper <- data.frame(
  "in_range" = c(29,56),
  "out" = c(39,70),
  row.names = c("chd","control"),
  stringsAsFactors = F
)
fisher.test(fisher_nonts_paper) # p=0.8

ks_plot_df_exp_nonTS = exp_in_range[,6]
ks_plot_df_exp_nonTS$type = "CHD cohort"
ks_plot_df_control_nonTS = control_in_range[,6]
ks_plot_df_control_nonTS$type = "Control cohort"

ks_plot_df_nonTS = rbind(ks_plot_df_exp_nonTS,ks_plot_df_control_nonTS)
ks_plot_df_nonTS$d_log = log10(ks_plot_df_nonTS$cen_els_pls_nonTS)

jpeg(file="allGelb_screen_v3_nonTS.jpg",width=7,height=4,units="in",res = 1200)
ggplot(ks_plot_df_nonTS, aes(x = d_log, group = type, color = type))+
    stat_ecdf(size=1) +
    scale_color_brewer(palette="Paired") +
    theme_classic(base_size = 12) +
    theme(legend.position ="right") +
    xlab("distance to closest cRE center") +
    ylab("fraction") +
    scale_x_continuous(breaks=c(3,4,5),labels= c("1Kbp","10Kbp","100Kbp"),limits = c(1, 5.1)) +
    ylab("fraction") +
    theme(legend.title=element_blank())
dev.off()




