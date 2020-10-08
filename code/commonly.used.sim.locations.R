library(tidyverse)
library(data.table)
library(valr)
setwd("~/Documents/projects/NemaScan_Performance/data")
rough.genome <- data.table::fread("genome.bed.tsv")%>%
  `colnames<-`(c("chrom","start","end"))

# Divergent Regions
divergent.regions <- data.table::fread("Common_divergent_regions_clustered.tsv") %>%
  `colnames<-`(c("chrom","start","end"))
int.divergent.regions <- data.table::fread("Intermediate_divergent_regions_clustered.tsv")
rare.divergent.regions <- data.table::fread("Rare_divergent_regions_clustered.tsv")
hard.filtered.variants <- data.table::fread("hard.filtered.snps.bim") %>%
  dplyr::filter(V1 != "MtDNA") %>%
  dplyr::mutate(end = V4) %>%
  dplyr::select(chrom = V1, start = V4, end) %>%
  dplyr::mutate(chrom = as.factor(chrom))
levels(hard.filtered.variants$chrom) <- c("X","I","II","III","IV","V")
non.divergent.regions <- valr::bed_subtract(rough.genome, divergent.regions)


write.table(non.divergent.regions, "nondivergent.regions.bed", quote = F, col.names = F, row.names = F)
write.table(divergent.regions, "common.divergent.regions.bed", quote = F, col.names = F, row.names = F)
write.table(rare.divergent.regions, "rare.divergent.regions.bed", quote = F, col.names = F, row.names = F)
write.table(int.divergent.regions, "int.divergent.regions.bed", quote = F, col.names = F, row.names = F)
write.table(hard.filtered.variants, "hard.filtered.variants.only.bed", quote = F, col.names = F, row.names = F)
