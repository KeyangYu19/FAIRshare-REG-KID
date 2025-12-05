library(readr)
library(dplyr)
library(tidyverse)
require(biomaRt)
#BiocManager::install("GenomicRanges")
library(GenomicRanges)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  stop("No file path provided!")
}

exp_in <- args[1]
cCREs <- args[2]
input_genes <- args[3]

read_f <- function(filename){ 
  data = read_delim(filename, delim="\t", escape_double = FALSE,  col_names = FALSE, comment = "#", trim_ws = TRUE)
  data = data[,c(1,2,4,5)]
  colnames(data) = c("chr","pos","ref","alt")
  data = data |> filter(ref %in% c("A","C","T","G") & alt %in% c("A","C","T","G"))
  data$chr = gsub("chr","", data$chr)
  return(data)
}

### Functions for calculation distance calculation
d_to_closest_gene_GR <- function(variants, genes) {
  
  variant_gr <- GRanges(
    seqnames = variants$chr,
    ranges = IRanges(start = variants$pos, end = variants$pos)
  )
  
  # Use start-10kbps and end+10kbps as region of focus
  gene_gr <- GRanges(
    seqnames = genes$chr,
    ranges = IRanges(start = genes$re_qstart_close, end = genes$re_qend_close)
  )
  
  nearest_gene <- GenomicRanges::distanceToNearest(variant_gr, gene_gr)
  
  d_gene <- rep(NA_real_, length(variant_gr))
  closest_gene <- rep(NA_character_, length(variant_gr))

  d_gene[queryHits(nearest_gene)] <- mcols(nearest_gene)$distance
  closest_gene[queryHits(nearest_gene)] <- genes$hgnc_symbol[subjectHits(nearest_gene)]
  
  variants$d_gene <- d_gene
  variants$gene <- closest_genes
  
  variants
}


d_to_closest_cRE_GR <- function(variants, REs) {
  
  variant_gr <- GRanges(
    seqnames = variants$chr,
    ranges = IRanges(start = variants$pos, end = variants$pos)
  )
  
  # Assume RE center as GR ref
  RE_gr <- GRanges(
    seqnames = REs$chr,
    ranges = IRanges(
      start = (REs$start + REs$end) / 2,
      end = (REs$start + REs$end) / 2
    )
  )
  
  nearest_RE <- GenomicRanges::distanceToNearest(variant_gr, RE_gr)
  
  distances <- rep(NA_real_, length(variant_gr))
  closest_cRE <- rep(NA_character_, length(variant_gr))

  distances[queryHits(nearest_RE)] <- mcols(nearest_RE)$distance
  closest_cRE[queryHits(nearest_RE)] <- RE_gr$id[subjectHits(nearest_gene)]
  
  variants$distance <- distances
  variants$cRE <- closest_cRE

  variants
}

exp_in = read_f(exp_in)
cRE_in = read_table2(cCREs, col_names = FALSE)
cRE_in$chr <- gsub("chr","",cRE_in$chr)

input_genes <- read_csv(input_genes, col_names = FALSE)


mart<-useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl",mirror = "www")
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


exp_in = exp_in[order(exp_all$chr,exp_all$pos),]
exp_in = exp_in[exp_all$chr != "X" & exp_all$chr != "Y",]
exp_in = d_to_closest_gene_GR(exp_in,out_genes_clean)
exp_in_range = exp_in[exp_in$d_gene == 0,] 


cRE_in <- cRE_in[sapply(1:nrow(cRE_in), function(i) {
  any(cRE_in$chr[i] == out_genes_clean$chr & cRE_in$start[i] >= out_genes_clean$re_qstart_close & cRE_in$end[i] <= out_genes_clean$re_qend_close)
}), ]
cRE_in <- unique(cRE_in)
exp_in_range = calculate_distance_to_RE_fast(exp_in_range, cRE_in)
exp_in_range = exp_in_range[order(exp_in_range$distance),]
exp_in_range_top = exp_in_range[exp_in_range$distance < 1000,]

write.table(exp_in_range, file = "exp.output.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(exp_in_range_top, file = "exp.topv.output.txt", sep = "\t", row.names = FALSE, col.names = TRUE)



