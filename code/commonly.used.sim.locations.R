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
  `colnames<-`(c("chrom","start","end")) %>%
  dplyr::mutate(chrom = as.factor(chrom))
levels(rare.divergent.regions$chrom) <- c(1,2,3,4,5,6)

arms.centers <- data.table::fread("data/ARMS_CENTERS.tsv") %>%
  `colnames<-`(c("chrom","start","end","region")) %>%
  dplyr::mutate(chrom = as.factor(chrom))
levels(arms.centers$chrom) <- c(1,2,3,4,5,6)

non.divergent.regions <- valr::bed_subtract(rough.genome, rare.divergent.regions) %>%
  valr::bed_subtract(., int.divergent.regions) %>%
  valr::bed_subtract(., divergent.regions)

all.divergent.regions <- divergent.regions %>%
  dplyr::full_join(., int.divergent.regions) %>%
  dplyr::full_join(., rare.divergent.regions)

write.table(non.divergent.regions, "data/nondivergent.regions.bed", quote = F, col.names = F, row.names = F)
write.table(divergent.regions, "data/common.divergent.regions.bed", quote = F, col.names = F, row.names = F)
write.table(rare.divergent.regions, "data/rare.divergent.regions.bed", quote = F, col.names = F, row.names = F)
write.table(int.divergent.regions, "data/int.divergent.regions.bed", quote = F, col.names = F, row.names = F)
write.table(hard.filtered.variants, "data/hard.filtered.variants.only.bed", quote = F, col.names = F, row.names = F)


divergent.regions.chrom.regions <- valr::bed_intersect(all.divergent.regions, arms.centers) %>%
  dplyr::select(chrom, start.x, end.x, region.y) %>%
  dplyr::rename(start = start.x, end = end.x, region = region.y) %>%
  dplyr::group_by(chrom, region) %>%
  tidyr::nest()

nondivergent.regions.chrom.regions <- valr::bed_intersect(non.divergent.regions, arms.centers) %>%
  dplyr::select(chrom, start.x, end.x, region.y) %>%
  dplyr::rename(start = start.x, end = end.x, region = region.y) %>%
  dplyr::group_by(chrom, region) %>%
  tidyr::nest()

write.beds.divergent <- function(data, chrom, region){
  data %>%
    dplyr::mutate(chrom = chrom) %>%
    dplyr::select(chrom, start, end) %>%
    dplyr::arrange(start) %>%
    write.table(., paste("data/Chr",chrom, region, "divergent.regions.bed", sep = "."), quote = F, col.names = F, row.names = F)
    
}
pmap(.l = list(divergent.regions.chrom.regions$data,
               divergent.regions.chrom.regions$chrom,
               divergent.regions.chrom.regions$region), write.beds.divergent)

write.beds.nondivergent <- function(data, chrom, region){
  data %>%
    dplyr::mutate(chrom = chrom) %>%
    dplyr::select(chrom, start, end) %>%
    dplyr::arrange(start) %>%
    write.table(., paste("data/Chr",chrom, region, "nondivergent.regions.bed", sep = "."), quote = F, col.names = F, row.names = F)
  
}
pmap(.l = list(nondivergent.regions.chrom.regions$data,
               nondivergent.regions.chrom.regions$chrom,
               nondivergent.regions.chrom.regions$region), write.beds.nondivergent)

int.divergent.chrom.regions <- valr::bed_intersect(int.divergent.regions, arms.centers) %>%
  dplyr::select(chrom, start.x, end.x, region.y) %>%
  dplyr::rename(start = start.x, end = end.x, region = region.y)

rare.chrom.regions <- valr::bed_intersect(rare.divergent.regions, arms.centers) %>%
  dplyr::select(chrom, start.x, end.x, region.y) %>%
  dplyr::rename(start = start.x, end = end.x, region = region.y)






hard.filtered.variants <- data.table::fread("data/hard.filtered.snps.bim") %>%
  dplyr::filter(V1 != "MtDNA") %>%
  dplyr::mutate(end = V4) %>%
  dplyr::select(chrom = V1, start = V4, end) %>%
  dplyr::mutate(chrom = as.factor(chrom))
levels(hard.filtered.variants$chrom) <- c("X","I","II","III","IV","V")
levels(hard.filtered.variants$chrom) <- c(6,1,2,3,4,5)
hard.filtered.variants <- hard.filtered.variants %>%
  dplyr::mutate(marker = paste(chrom,start,sep = ":"))


common.div.regions.hard.filtered <- valr::bed_intersect(hard.filtered.variants, divergent.regions)
int.div.regions.hard.filtered <- valr::bed_intersect(hard.filtered.variants, int.divergent.regions)
rare.div.regions.hard.filtered <- valr::bed_intersect(hard.filtered.variants, rare.divergent.regions)


hard.filtered.chrom.armscenters.nested <- valr::bed_intersect(hard.filtered.variants, arms.centers) %>%
  dplyr::select(chrom, start.x, end.x, region.y) %>%
  dplyr::rename(start = start.x, end = end.x, region = region.y) %>%
  dplyr::group_by(chrom, region) %>%
  tidyr::nest()

markers <- hard.filtered.chrom.armscenters.nested$data[[1]]
div.regions <- divergent.regions.chrom.regions$data[[1]]
non.div.regions <- nondivergent.regions.chrom.regions$data[[1]]
chrom <- hard.filtered.chrom.armscenters.nested$chrom[[1]]
region <- hard.filtered.chrom.armscenters.nested$region[[1]]

pct.div.markers <- function(markers, div.regions, non.div.regions, chrom, region){

  markers.C <- markers %>%
    dplyr::mutate(chrom = chrom) %>%
    dplyr::select(chrom, start, end)
  divergent.C <- div.regions %>%
    dplyr::mutate(chrom = chrom) %>%
    dplyr::select(chrom, start, end)
  non.divergent.C <- non.div.regions %>%
    dplyr::mutate(chrom = chrom) %>%
    dplyr::select(chrom, start, end)
  
  div.markers <- valr::bed_intersect(markers.C, divergent.C) %>%
    dplyr::select(chrom, start.x) %>%
    tidyr::unite("marker", chrom:start.x, sep = ":") %>%
    nrow()
  non.div.markers <- valr::bed_intersect(markers.C, non.divergent.C) %>%
    dplyr::select(chrom, start.x) %>%
    tidyr::unite("marker", chrom:start.x, sep = ":") %>%
    nrow()
  
  data.frame(div.markers, non.div.markers) %>%
    dplyr::mutate(chrom = chrom, region = region, pct.div = round(div.markers/(div.markers + non.div.markers), 3),
                  reps.non.divergent = 100,
                  reps.divergent = if_else(condition = (reps.non.divergent*pct.div)*3 > 100, 
                                           true = 100,
                                           false = (reps.non.divergent*pct.div)*3),
                  reps.divergent = if_else(condition = reps.divergent < 30, 
                                           true = 30,
                                           false = round(reps.divergent,0)))
  
}
pmap(.l = list(hard.filtered.chrom.armscenters.nested$data,
               divergent.regions.chrom.regions$data,
               nondivergent.regions.chrom.regions$data,
               hard.filtered.chrom.armscenters.nested$chrom,
               hard.filtered.chrom.armscenters.nested$region), pct.div.markers) %>%
  Reduce(rbind,.) %>%
  tidyr::pivot_longer(cols = contains("reps"), names_to = "type", values_to = "reps") %>%
  dplyr::mutate(type = gsub(type, pattern = "reps.", replacement = "")) %>%
  data.frame()


# 
# imputed.variants <- data.table::fread("data/imputed.snps.bim") %>%
#   dplyr::filter(V1 != "MtDNA") %>%
#   dplyr::mutate(end = V4) %>%
#   dplyr::select(chrom = V1, start = V4, end) %>%
#   dplyr::mutate(chrom = as.factor(chrom))
# levels(imputed.variants$chrom) <- c("X","I","II","III","IV","V")
# levels(imputed.variants$chrom) <- c(6,1,2,3,4,5)
# imputed.variants <- imputed.variants %>%
#   dplyr::mutate(marker = paste(chrom,start,sep = ":"))
# 

# hard.filtered.variants.annot <- hard.filtered.variants %>%
#   dplyr::mutate(in.imputed = hard.filtered.variants$marker %in% imputed.variants$marker,
#                 in.imputed = if_else(condition = in.imputed == TRUE, 
#                                      true = "Imputed", 
#                                      false = "Hard-Filtered Only"),
#                 common.div = hard.filtered.variants$marker %in% common.div.regions.hard.filtered$marker.x,
#                 int.div = hard.filtered.variants$marker %in% int.div.regions.hard.filtered$marker.x,
#                 rare.div = hard.filtered.variants$marker %in% rare.div.regions.hard.filtered$marker.x) %>%
#   tidyr::pivot_longer(cols = c(common.div, int.div, rare.div), 
#                       names_to = "frequency",
#                       values_to = "divergent",
#                       names_repair = "unique") %>%
#   dplyr::mutate(divergent = if_else(divergent == TRUE, true = "Divergent", false = "Non-Divergent"))
# hard.filtered.variants.annot$frequency <- case_when(
#   hard.filtered.variants.annot$frequency == "common.div" ~ "Common",
#   hard.filtered.variants.annot$frequency == "int.div" ~ "Intermediate",
#   hard.filtered.variants.annot$frequency == "rare.div" ~ "Rare")
# levels(hard.filtered.variants.annot$chrom) <- c("X","I","II","III","IV","V")
# hard.filtered.variants.annot$chrom <- factor(hard.filtered.variants.annot$chrom, 
#                                              levels = levels(hard.filtered.variants.annot$chrom)[c(2:6,1)])
# hard.filtered.variants.annot$marker <- paste(hard.filtered.variants.annot$chrom, 
#                                              hard.filtered.variants.annot$start, 
#                                              sep = ":")
# variant.set.totals <- hard.filtered.variants.annot %>%
#   dplyr::group_by(in.imputed, chrom) %>%
#   dplyr::summarise(n()) %>%
#   `colnames<-`(c("in.imputed","chrom","total"))
# 
# hard.filtered.variants.annot %>%
#   dplyr::group_by(in.imputed, chrom, frequency, divergent) %>%
#   dplyr::summarise(n()) %>%
#   dplyr::full_join(.,variant.set.totals) %>%
#   dplyr::mutate(cat.perc = `n()`/total) %>%
#   dplyr::filter(divergent == "Divergent") %>%
#   ggplot(., mapping = aes(x = in.imputed, y = cat.perc*100, fill = in.imputed)) + 
#   theme_bw() + 
#   geom_col() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
#         axis.title.x = element_blank()) + 
#   scale_fill_manual(values = nationalparkcolors::park_palettes$RockyMountains[c(2,4)], name = "Variant Set") + 
#   facet_grid(frequency~chrom) +
#   labs(y = "Divergent Variants (%)") + 
#   ggsave("output/perc.divergent.variants.in.variant.set.png", width = 8, height = 5.5)
# 
# hard.filtered.variants.annot %>%
#   dplyr::group_by(in.imputed, chrom, frequency, divergent) %>%
#   dplyr::summarise(n()) %>%
#   dplyr::full_join(.,variant.set.totals) %>%
#   dplyr::mutate(cat.perc = `n()`/total) %>%
#   dplyr::filter(divergent == "Divergent") %>%
#   ggplot(., mapping = aes(x = in.imputed, y = `n()`, fill = in.imputed)) + 
#   theme_bw() + 
#   geom_col() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
#         axis.title.x = element_blank()) + 
#   scale_fill_manual(values = nationalparkcolors::park_palettes$RockyMountains[c(2,4)], name = "Variant Set") + 
#   facet_grid(frequency~chrom) +
#   labs(y = "Divergent Variants (%)") + 
#   ggsave("output/perc.divergent.variants.in.variant.set.png", width = 8, height = 5.5)
# 
# hf.div.nest <- hard.filtered.variants.annot %>%
#   dplyr::arrange(chrom) %>%
#   dplyr::group_by(frequency, chrom, in.imputed) %>%
#   tidyr::nest()
# 
# chrom.breaks <- seq(0, 21000000, by = 100000)
# binned.divergence <- function(chrom, in.imputed, freq,  div){
#   div.binned <- div %>%
#     dplyr::mutate(div.bins = cut(start, breaks = chrom.breaks, include.lowest = T)) %>%
#     dplyr::group_by(div.bins) %>%
#     tidyr::nest()
#   
#   bin.counts <- list()
#   for(i in 1:length(div.binned$data)){
#     bin.counts[[i]] <- div.binned$data[[i]] %>%
#       dplyr::group_by(divergent) %>%
#       dplyr::summarise(n()) %>%
#       `colnames<-`(c("divergent","frequency")) %>%
#       dplyr::mutate(bin = div.binned$div.bins[[i]]) %>%
#       suppressMessages()
#   }
#   bin.counts.df <- Reduce(rbind, bin.counts) %>%
#     dplyr::mutate(chrom = chrom,
#                   in.imputed = in.imputed,
#                   freq = freq)
#   
# }
# hf.div.binned.df <- pmap(.l = list(hf.div.nest$chrom,
#                                    hf.div.nest$in.imputed,
#                                    hf.div.nest$frequency,
#                                    hf.div.nest$data), 
#                          .f = binned.divergence) %>%
#   Reduce(rbind,.)
# beep(sound = "coin")
# 
# 
# bin.sums <- hf.div.binned.df %>%
#   dplyr::group_by(chrom, in.imputed, bin, freq) %>%
#   dplyr::summarise(sum(frequency))
# 
# hf.div.binned.df %>%
#   dplyr::full_join(., bin.sums) %>%
#   dplyr::mutate(pct = frequency/`sum(frequency)`) %>%
#   tidyr::pivot_wider(names_from = in.imputed, values_from = pct) %>%
#   View()
#   ggplot(., mapping = aes(x = bin, y = pct, fill = divergent)) + 
#   theme_bw() + 
#   geom_point(shape = 21) + 
#   facet_grid(freq + in.imputed ~ chrom) + 
#   theme(axis.text.x = element_blank())
# 
# 
# bed_intersect(int.divergent.regions, divergent.regions)
# 
#   

