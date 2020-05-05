# This script is to call DE genes for rice TF loss-of-function or overexpression RNASeq 
# data. It's based on DESeq2. 
# The raw count was processed from SRA-->fastp-->STAR. Reads were mapped to MSU7.
# See `snakemake` pipelines for processing details.

# Author: Ji Huang
# Date: 2020-05-04

# Rice gene: bHLH148 (LOC_Os03g53020).
# Paper: NA


# 0. Prep ---------------------------------------------------------------------------

library(here)
library(magrittr)
library(tidyverse)
library(DESeq2)

# Specify SRP numeber and DESeq2 pvalue cutoff.
SRP_num <- "SRP052309"
p_cutoff <- 0.05
gene_name <- "BHLH148"
loc_id <- "LOC_Os03g53020"



# 1. Read experiment meta-data ------------------------------------------------------

lib_table <- read_tsv(here("data", "20200430_11SRP_RNASeq",
                           "SraRunTable_20200430_11SRP.tsv")) %>% 
  filter(`SRA Study` == SRP_num, `Assay Type` == "RNA-Seq")
  

## manually change
lib_table <- lib_table %>% select(Run, Genotype, treatment) %>% 
  mutate(genotype = case_when(Genotype == "wild type" ~ "wt",
                              Genotype == "OsbHLH148 loss-of-function" ~ "bhlh148_mut")) %>% 
  mutate(condition = paste0(genotype, "_", treatment)) %>% 
  mutate(condition = factor(condition, 
                            levels=c("wt_control", "wt_drought",
                                     "bhlh148_mut_control", "bhlh148_mut_drought")))


# 2. Read raw count -----------------------------------------------------------------


raw_data <-read.table(here("data", "20200430_11SRP_RNASeq", 
                           "star_osa_0430_11SRP_SE.featureCount.gz"),
                       sep="\t",header=T,row.names=1)

raw_data[1:5] <- NULL

colnames(raw_data) %<>% str_extract("SRR.+Ali") %<>% str_replace("Ali", "")


raw_data <- raw_data %>% select(lib_table$Run)


# 3. Call DE genes by DESeq2 --------------------------------------------------------

dds_rice <- DESeqDataSetFromMatrix(countData = raw_data,
                                   colData = lib_table,
                                   design = ~ 0+condition)

keep <- rowSums(counts(dds_rice)) >= 10
sum(keep)

dds_rice <- dds_rice[keep,]

dds <- DESeq(dds_rice, quiet = F)
resultsNames(dds)

de_mutVSwt_control <- results(dds, contrast = c(-1,0,1,0), tidy = T, 
                alpha = p_cutoff) %>% filter(padj < p_cutoff)

de_mutVSwt_drought <- results(dds, contrast = c(0,-1,0,1), tidy = T, 
                              alpha = p_cutoff) %>% filter(padj < p_cutoff)


de_mutVSwt_interaction <- results(dds, contrast = c(1,-1,-1,1), tidy = T, 
                              alpha = p_cutoff) %>% filter(padj < p_cutoff)




# 4. Save result --------------------------------------------------------------------
fn_mutVSwt_control <- paste0("rice_", gene_name, "_", loc_id, "_lof_DESeq2_", nrow(de_mutVSwt_control), 
                    "_mutVSwt_control",".tsv")
write_tsv(de_mutVSwt_control, here("result", "rnaseq","DE_result", fn_mutVSwt_control))


fn_mutVSwt_drought <- paste0("rice_", gene_name, "_", loc_id, "_lof_DESeq2_", nrow(de_mutVSwt_drought), 
                             "_mutVSwt_drought",".tsv")
write_tsv(de_mutVSwt_drought, here("result", "rnaseq","DE_result", fn_mutVSwt_drought))


fn_mutVSwt_interaction <- paste0("rice_", gene_name, "_", loc_id, "_lof_DESeq2_", nrow(de_mutVSwt_interaction), 
                             "_mutVSwt_interaction",".tsv")
write_tsv(de_mutVSwt_interaction, here("result", "rnaseq","DE_result", fn_mutVSwt_interaction))

# 5. Quick PCA plot -----------------------------------------------------------------

vsd <- vst(dds_rice, blind=FALSE)
plotPCA(vsd, intgroup=c("condition"))


# 6. Session info -------------------------------------------------------------------

sessionInfo()

# R version 3.6.1 (2019-07-05)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 18363)
# 
# Matrix products: default
# 
# locale:
#   [1] LC_COLLATE=English_United States.1252 
# [2] LC_CTYPE=English_United States.1252   
# [3] LC_MONETARY=English_United States.1252
# [4] LC_NUMERIC=C                          
# [5] LC_TIME=English_United States.1252    
# 
# attached base packages:
#   [1] parallel  stats4    stats     graphics  grDevices utils     datasets 
# [8] methods   base     
# 
# other attached packages:
#   [1] DESeq2_1.24.0               SummarizedExperiment_1.14.1
# [3] DelayedArray_0.10.0         BiocParallel_1.18.1        
# [5] matrixStats_0.55.0          Biobase_2.44.0             
# [7] GenomicRanges_1.36.1        GenomeInfoDb_1.20.0        
# [9] IRanges_2.18.2              S4Vectors_0.22.1           
# [11] BiocGenerics_0.30.0         forcats_0.4.0              
# [13] stringr_1.4.0               dplyr_0.8.3                
# [15] purrr_0.3.2                 readr_1.3.1                
# [17] tidyr_1.0.0                 tibble_2.1.3               
# [19] ggplot2_3.2.1               tidyverse_1.2.1            
# [21] magrittr_1.5                here_0.1                   
# 
# loaded via a namespace (and not attached):
#   [1] nlme_3.1-141           bitops_1.0-6           bit64_0.9-7           
# [4] lubridate_1.7.4        RColorBrewer_1.1-2     httr_1.4.1            
# [7] rprojroot_1.3-2        tools_3.6.1            backports_1.1.4       
# [10] R6_2.4.0               rpart_4.1-15           Hmisc_4.2-0           
# [13] DBI_1.0.0              lazyeval_0.2.2         colorspace_1.4-1      
# [16] nnet_7.3-12            withr_2.1.2            tidyselect_0.2.5      
# [19] gridExtra_2.3          bit_1.1-14             compiler_3.6.1        
# [22] cli_1.1.0              rvest_0.3.4            htmlTable_1.13.2      
# [25] xml2_1.2.2             labeling_0.3           scales_1.0.0          
# [28] checkmate_1.9.4        genefilter_1.66.0      digest_0.6.21         
# [31] foreign_0.8-72         XVector_0.24.0         base64enc_0.1-3       
# [34] pkgconfig_2.0.3        htmltools_0.3.6        htmlwidgets_1.3       
# [37] rlang_0.4.2            readxl_1.3.1           rstudioapi_0.10       
# [40] RSQLite_2.1.2          generics_0.0.2         jsonlite_1.6          
# [43] acepack_1.4.1          RCurl_1.95-4.12        GenomeInfoDbData_1.2.1
# [46] Formula_1.2-3          Matrix_1.2-17          Rcpp_1.0.2            
# [49] munsell_0.5.0          lifecycle_0.1.0        stringi_1.4.3         
# [52] zlibbioc_1.30.0        blob_1.2.0             grid_3.6.1            
# [55] crayon_1.3.4           lattice_0.20-38        haven_2.1.1           
# [58] splines_3.6.1          annotate_1.62.0        hms_0.5.1             
# [61] locfit_1.5-9.1         zeallot_0.1.0          knitr_1.25            
# [64] pillar_1.4.2           geneplotter_1.62.0     XML_3.98-1.20         
# [67] glue_1.3.1             latticeExtra_0.6-28    data.table_1.12.2     
# [70] modelr_0.1.5           vctrs_0.2.0            cellranger_1.1.0      
# [73] gtable_0.3.0           assertthat_0.2.1       xfun_0.9              
# [76] xtable_1.8-4           broom_0.5.2            survival_2.44-1.1     
# [79] memoise_1.1.0          AnnotationDbi_1.46.1   cluster_2.1.0      
