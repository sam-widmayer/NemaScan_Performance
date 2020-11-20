library(tidyverse)
library(data.table)
library(valr)
library(nationalparkcolors)
library(beepr)
setwd("~/Documents/projects/NemaScan_Performance")

rough.genome <- data.table::fread("data/genome.bed.tsv")%>%
  `colnames<-`(c("chrom","start","end")) %>%
  dplyr::mutate(chrom = as.factor(chrom))
levels(rough.genome$chrom) <- c(1,2,3,4,5,6)

# Divergent Regions
divergent.regions <- data.table::fread("data/Common_divergent_regions_clustered.tsv") %>%
  `colnames<-`(c("chrom","start","end"))%>%
  dplyr::mutate(chrom = as.factor(chrom))
levels(divergent.regions$chrom) <- c(1,2,3,4,5,6)

int.divergent.regions <- data.table::fread("data/Intermediate_divergent_regions_clustered.tsv") %>%
  `colnames<-`(c("chrom","start","end"))%>%
  dplyr::mutate(chrom = as.factor(chrom))
levels(int.divergent.regions$chrom) <- c(1,2,3,4,5,6)

rare.divergent.regions <- data.table::fread("data/Rare_divergent_regions_clustered.tsv") %>%
  `colnames<-`(c("chrom","start","end"))%>%
  dplyr::mutate(chrom = as.factor(chrom))
levels(rare.divergent.regions$chrom) <- c(1,2,3,4,5,6)






hard.filtered.variants <- data.table::fread("data/hard.filtered.snps.bim") %>%
  dplyr::filter(V1 != "MtDNA") %>%
  dplyr::mutate(end = V4) %>%
  dplyr::select(chrom = V1, start = V4, end) %>%
  dplyr::mutate(chrom = as.factor(chrom))
levels(hard.filtered.variants$chrom) <- c("X","I","II","III","IV","V")
levels(hard.filtered.variants$chrom) <- c(6,1,2,3,4,5)
hard.filtered.variants <- hard.filtered.variants %>%
  dplyr::mutate(marker = paste(chrom,start,sep = ":"))

imputed.variants <- data.table::fread("data/imputed.snps.bim") %>%
  dplyr::filter(V1 != "MtDNA") %>%
  dplyr::mutate(end = V4) %>%
  dplyr::select(chrom = V1, start = V4, end) %>%
  dplyr::mutate(chrom = as.factor(chrom))
levels(imputed.variants$chrom) <- c("X","I","II","III","IV","V")
levels(imputed.variants$chrom) <- c(6,1,2,3,4,5)
imputed.variants <- imputed.variants %>%
  dplyr::mutate(marker = paste(chrom,start,sep = ":"))

common.div.regions.hard.filtered <- valr::bed_intersect(hard.filtered.variants, divergent.regions)
int.div.regions.hard.filtered <- valr::bed_intersect(hard.filtered.variants, int.divergent.regions)
rare.div.regions.hard.filtered <- valr::bed_intersect(hard.filtered.variants, rare.divergent.regions)

hard.filtered.variants.annot <- hard.filtered.variants %>%
  dplyr::mutate(in.imputed = hard.filtered.variants$marker %in% imputed.variants$marker,
                in.imputed = if_else(condition = in.imputed == TRUE, 
                                     true = "Imputed", 
                                     false = "Hard-Filtered Only"),
                common.div = hard.filtered.variants$marker %in% common.div.regions.hard.filtered$marker.x,
                int.div = hard.filtered.variants$marker %in% int.div.regions.hard.filtered$marker.x,
                rare.div = hard.filtered.variants$marker %in% rare.div.regions.hard.filtered$marker.x) %>%
  tidyr::pivot_longer(cols = c(common.div, int.div, rare.div), 
                      names_to = "frequency",
                      values_to = "divergent",
                      names_repair = "unique") %>%
  dplyr::mutate(divergent = if_else(divergent == TRUE, true = "Divergent", false = "Non-Divergent"))
hard.filtered.variants.annot$frequency <- case_when(
  hard.filtered.variants.annot$frequency == "common.div" ~ "Common",
  hard.filtered.variants.annot$frequency == "int.div" ~ "Intermediate",
  hard.filtered.variants.annot$frequency == "rare.div" ~ "Rare")
levels(hard.filtered.variants.annot$chrom) <- c("X","I","II","III","IV","V")
hard.filtered.variants.annot$chrom <- factor(hard.filtered.variants.annot$chrom, 
                                             levels = levels(hard.filtered.variants.annot$chrom)[c(2:6,1)])
hard.filtered.variants.annot$marker <- paste(hard.filtered.variants.annot$chrom, 
                                             hard.filtered.variants.annot$start, 
                                             sep = ":")
variant.set.totals <- hard.filtered.variants.annot %>%
  dplyr::group_by(in.imputed, chrom) %>%
  dplyr::summarise(n()) %>%
  `colnames<-`(c("in.imputed","chrom","total"))

hard.filtered.variants.annot %>%
  dplyr::group_by(in.imputed, chrom, frequency, divergent) %>%
  dplyr::summarise(n()) %>%
  dplyr::full_join(.,variant.set.totals) %>%
  dplyr::mutate(cat.perc = `n()`/total) %>%
  dplyr::filter(divergent == "Divergent") %>%
  ggplot(., mapping = aes(x = in.imputed, y = cat.perc*100, fill = in.imputed)) + 
  theme_bw() + 
  geom_col() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.title.x = element_blank()) + 
  scale_fill_manual(values = nationalparkcolors::park_palettes$RockyMountains[c(2,4)], name = "Variant Set") + 
  facet_grid(frequency~chrom) +
  labs(y = "Divergent Variants (%)") + 
  ggsave("output/perc.divergent.variants.in.variant.set.png", width = 8, height = 5.5)

hard.filtered.variants.annot %>%
  dplyr::group_by(in.imputed, chrom, frequency, divergent) %>%
  dplyr::summarise(n()) %>%
  dplyr::full_join(.,variant.set.totals) %>%
  dplyr::mutate(cat.perc = `n()`/total) %>%
  dplyr::filter(divergent == "Divergent") %>%
  ggplot(., mapping = aes(x = in.imputed, y = `n()`, fill = in.imputed)) + 
  theme_bw() + 
  geom_col() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.title.x = element_blank()) + 
  scale_fill_manual(values = nationalparkcolors::park_palettes$RockyMountains[c(2,4)], name = "Variant Set") + 
  facet_grid(frequency~chrom) +
  labs(y = "Divergent Variants (%)") + 
  ggsave("output/perc.divergent.variants.in.variant.set.png", width = 8, height = 5.5)

hf.div.nest <- hard.filtered.variants.annot %>%
  dplyr::arrange(chrom) %>%
  dplyr::group_by(frequency, chrom, in.imputed) %>%
  tidyr::nest()

chrom.breaks <- seq(0, 21000000, by = 100000)
binned.divergence <- function(chrom, in.imputed, freq,  div){
  div.binned <- div %>%
    dplyr::mutate(div.bins = cut(start, breaks = chrom.breaks, include.lowest = T)) %>%
    dplyr::group_by(div.bins) %>%
    tidyr::nest()
  
  bin.counts <- list()
  for(i in 1:length(div.binned$data)){
    bin.counts[[i]] <- div.binned$data[[i]] %>%
      dplyr::group_by(divergent) %>%
      dplyr::summarise(n()) %>%
      `colnames<-`(c("divergent","frequency")) %>%
      dplyr::mutate(bin = div.binned$div.bins[[i]]) %>%
      suppressMessages()
  }
  bin.counts.df <- Reduce(rbind, bin.counts) %>%
    dplyr::mutate(chrom = chrom,
                  in.imputed = in.imputed,
                  freq = freq)
  
}
hf.div.binned.df <- pmap(.l = list(hf.div.nest$chrom,
                                   hf.div.nest$in.imputed,
                                   hf.div.nest$frequency,
                                   hf.div.nest$data), 
                         .f = binned.divergence) %>%
  Reduce(rbind,.)
beep(sound = "coin")


bin.sums <- hf.div.binned.df %>%
  dplyr::group_by(chrom, in.imputed, bin, freq) %>%
  dplyr::summarise(sum(frequency))

hf.div.binned.df %>%
  dplyr::full_join(., bin.sums) %>%
  dplyr::mutate(pct = frequency/`sum(frequency)`) %>%
  tidyr::pivot_wider(names_from = in.imputed, values_from = pct) %>%
  View()
  ggplot(., mapping = aes(x = bin, y = pct, fill = divergent)) + 
  theme_bw() + 
  geom_point(shape = 21) + 
  facet_grid(freq + in.imputed ~ chrom) + 
  theme(axis.text.x = element_blank())



bed_intersect(int.divergent.regions, divergent.regions)

  
non.divergent.regions <- valr::bed_subtract(rough.genome, rare.divergent.regions) %>%
  valr::bed_subtract(., int.divergent.regions) %>%
  valr::bed_subtract(., divergent.regions)
  



write.table(non.divergent.regions, "data/nondivergent.regions.bed", quote = F, col.names = F, row.names = F)
write.table(divergent.regions, "data/common.divergent.regions.bed", quote = F, col.names = F, row.names = F)
write.table(rare.divergent.regions, "data/rare.divergent.regions.bed", quote = F, col.names = F, row.names = F)
write.table(int.divergent.regions, "data/int.divergent.regions.bed", quote = F, col.names = F, row.names = F)
write.table(hard.filtered.variants, "data/hard.filtered.variants.only.bed", quote = F, col.names = F, row.names = F)
