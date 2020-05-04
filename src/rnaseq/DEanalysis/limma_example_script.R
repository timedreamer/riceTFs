library(here)
library(limma)
library(magrittr)
library(tidyverse)
library(edgeR)


lib_table <- read_tsv(here("data", "20200430_11SRP_RNASeq",
                           "SraRunTable_20200430_11SRP.tsv")) %>% 
  filter(`SRA Study` == "SRP118965", `Assay Type` == "RNA-Seq") %>% 
  select(Run, tissue, Genotype)

lib_table <- lib_table %>% mutate(genotype = case_when(Genotype == "wild type" ~ "wt",
                                                       Genotype == "wox11 mutant" ~ "wox11_mut")) %>% 
  filter(genotype != "NA")


raw_data <-read.table(here("data", "20200430_11SRP_RNASeq", 
                           "star_osa_0430_11SRP_PE.featureCount.gz"),
                      sep="\t",header=T,row.names=1)

raw_data[1:5] <- NULL

colnames(raw_data) %<>% str_extract("SRR.+Ali") %<>% 
  str_replace("Ali", "")


raw_data <- raw_data %>% select(lib_table$Run)

design <- model.matrix(~0+genotype, data = lib_table)

y <- DGEList(counts = raw_data)
keep <- filterByExpr(y)
sum(keep)
y <- y[keep,,keep.lib.size=FALSE]

background <- y$genes$genes

v <- voom(y, design, plot=T, normalize="quantile")
plotMDS(v)

fit <- lmFit(v, design)
contrast_m <- makeContrasts(genotypewox11_mut-genotypewt, levels=design)
fit2 <- contrasts.fit(fit, contrast_m)
fit2 <- eBayes(fit2)
summary(decideTests(fit2))

de_gene <- topTable(fit2, p.value = 0.05, number = Inf, 
                    adjust.method= "BH") %>% rownames_to_column(var = "geneid")

file_name <- paste0("rice_wox11_lof_limma_", nrow(de_gene), ".tsv")
write_tsv(de_gene, here("result", "rnaseq","DE_result", file_name))
