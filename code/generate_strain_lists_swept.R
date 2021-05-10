library(tidyverse)
require(ids)
require(data.table)
require(RColorBrewer)
setwd("~/Documents/projects/NemaScan_Performance/")
chrom.swept.bin <- function(data, chrom){
  lowest <- data %>%
    dplyr::filter(pct.swept < 0.20) %>%
    dplyr::mutate(swept.bin = "0-0.20") 
  
  data %>%
    dplyr::filter(pct.swept > 0.20) %>%
    mutate(swept.bin = cut(pct.swept, c(0.2, 0.5, 1), include.lowest = F)) %>%
    dplyr::mutate(swept.bin = gsub(swept.bin, pattern = "\\(", replacement = "") %>%
                    gsub(., pattern = "\\]", replacement = "") %>%
                    gsub(., pattern = ",", replacement = "-")) %>%
    dplyr::full_join(., lowest) %>%
    dplyr::mutate(chrom = chrom, 
                  swept.bin = as.factor(swept.bin))
}
chrom.swept.liberal.bin <- function(data, chrom){
  lowest <- data %>%
    dplyr::filter(pct.swept < 0.30) %>%
    dplyr::mutate(swept.bin = "0-0.30") 
  
  data %>%
    dplyr::filter(pct.swept > 0.30) %>%
    mutate(swept.bin = cut(pct.swept, c(0.3, 1), include.lowest = F)) %>%
    dplyr::mutate(swept.bin = gsub(swept.bin, pattern = "\\(", replacement = "") %>%
                    gsub(., pattern = "\\]", replacement = "") %>%
                    gsub(., pattern = ",", replacement = "-")) %>%
    dplyr::full_join(., lowest) %>%
    dplyr::mutate(chrom = chrom, 
                  swept.bin = as.factor(swept.bin))
}
# setwd("~/Documents/AndersenLab/NemaScan_Performance/")
strain.data <- data.table::fread("data/CelegansStrainData.tsv")
isotypes.20210121 <- data.table::fread("data/complete.isotypes.20210121.txt", header = F) %>%
  dplyr::rename(isotype = V1)
sweeps <- data.table::fread("data/sweep_summary_20210121.tsv")
hahnel_204 <- data.table::fread("data/hahnel_209.tsv") %>%
  dplyr::select(strain)
today <- format(Sys.time(), '%Y%m%d')

sweptness.nested <- sweeps %>%
  tidyr::pivot_longer(cols = -isotype, names_to = "chrom", values_to = "pct.swept") %>%
  dplyr::group_by(chrom) %>%
  tidyr::nest()



swept.coarse.bins <- purrr::map2(sweptness.nested$data, sweptness.nested$chrom, chrom.swept.bin) %>%
  Reduce(rbind, .) %>%
  dplyr::filter(!chrom %in% c("II", "III")) %>%
  dplyr::group_by(isotype, swept.bin) %>%
  dplyr::summarise(n()) %>%
  tidyr::pivot_wider(names_from = swept.bin, values_from = `n()`) 

conservative.unswept.isotypes <- swept.coarse.bins %>% 
  dplyr::filter(`0-0.20` >= 3) %>%
  dplyr::select(isotype) %>%
  dplyr::distinct(.)
genome.unswept.list <- paste(conservative.unswept.isotypes$isotype,sep = "",collapse = ",")

conservative.swept.isotypes <- swept.coarse.bins %>% 
  dplyr::filter(`0.5-1` >= 3) %>%
  dplyr::select(isotype) %>%
  dplyr::distinct(.)
genome.swept.list <- paste(conservative.swept.isotypes$isotype,sep = "",collapse = ",")




lib.swept.coarse.bins <- purrr::map2(sweptness.nested$data, sweptness.nested$chrom, chrom.swept.liberal.bin) %>%
  Reduce(rbind, .) %>%
  dplyr::filter(!chrom %in% c("II", "III")) %>%
  dplyr::group_by(isotype, swept.bin) %>%
  dplyr::summarise(n()) %>%
  tidyr::pivot_wider(names_from = swept.bin, values_from = `n()`) 

liberal.swept.isotypes <- lib.swept.coarse.bins %>% 
  dplyr::filter(`0.3-1` >= 1) %>%
  dplyr::select(isotype) %>%
  dplyr::distinct(.)
liberal.swept.isotypes.list <- paste(liberal.swept.isotypes$isotype,sep = "",collapse = ",")
write.table(data.frame("liberal.swept", liberal.swept.isotypes.list),
            file = "output/liberal.swept.isotypes.txt",
            quote = F, row.names = F, col.names = F)

liberal.unswept.isotypes <- lib.swept.coarse.bins %>% 
  dplyr::select(isotype) %>%
  dplyr::filter(!isotype %in% liberal.swept.isotypes$isotype) %>%
  dplyr::distinct(.)
liberal.unswept.isotypes.list <- paste(liberal.unswept.isotypes$isotype,sep = "",collapse = ",")
write.table(data.frame("liberal.unswept", liberal.unswept.isotypes.list),
            file = "output/liberal.unswept.isotypes.txt",
            quote = F, row.names = F, col.names = F)



  
# sweeps %>%
#   dplyr::mutate(sweep.group = case_when(isotype %in% conservative.swept.isotypes$isotype ~ "swept",
#                                         isotype %in% conservative.unswept.isotypes$isotype ~ "unswept",
#                                         TRUE ~ "NA")) %>%
#   tidyr::pivot_longer(cols = c("I","II","III","IV","V","X"), names_to = "CHROM", values_to = "pct.swept") %>%
#   dplyr::mutate(dummy = "dummy") %>%
#   dplyr::arrange(isotype) %>%
#   ggplot(., mapping = aes(x = dummy, y = isotype, fill = pct.swept)) +
#   theme_bw(base_size = 9) + 
#   geom_tile() + 
#   facet_grid(sweep.group ~ CHROM, scales = "free", space = "free", drop = TRUE) + 
#   theme(axis.text.x = element_blank(),
#         axis.ticks.x = element_blank(),
#         axis.title.x = element_blank()) + 
#   scale_fill_distiller(palette = "PRGn") + 
#   ggsave(filename = "output/sweep.pops.png", height = 40, width = 8)


# Joining Strain Info
metadata <- sweeps %>%
  dplyr::full_join(., strain.data) %>%
  dplyr::mutate(I = dplyr::if_else(condition = I == TRUE, 
                                   true = "Swept", 
                                   false = "Unswept"),
                II = dplyr::if_else(condition = II == TRUE, 
                                    true = "Swept", 
                                    false = "Unswept"),
                III = dplyr::if_else(condition = III == TRUE, 
                                     true = "Swept", 
                                     false = "Unswept"),
                IV = dplyr::if_else(condition = IV == TRUE, 
                                    true = "Swept", 
                                    false = "Unswept"),
                V = dplyr::if_else(condition = V == TRUE, 
                                   true = "Swept", 
                                   false = "Unswept"),
                X = dplyr::if_else(condition = X == TRUE, 
                                   true = "Swept", 
                                   false = "Unswept"))

complete <- metadata %>%
  dplyr::select(isotype) %>%
  dplyr::filter(!duplicated(isotype))
complete.list <- paste(complete$isotype,sep = "",collapse = ",")

complete.20210121 <- isotypes.20210121 %>%
  dplyr::select(isotype) %>%
  dplyr::filter(!duplicated(isotype))
complete.20210121.list <- paste(complete.20210121$isotype,sep = "",collapse = ",")



# # Chromosome-specific sweep-based strain lists
# chrom.swept.list <- purrr::map2(sweptness.nested$data, 
#                                 sweptness.nested$chrom, 
#                                 chrom.swept.bin)
# sample.sizes <- list()
# for(i in 1:length(chrom.swept.list)){
#   conservative.swept <- chrom.swept.list[[i]] %>%
#     dplyr::filter(swept.bin == "0.5-1") %>%
#     dplyr::select(isotype, pct.swept, chrom) %>%
#     as.data.frame()
#   
#   conservative.unswept <- chrom.swept.list[[i]] %>%
#     dplyr::filter(swept.bin == "0-0.20") %>%
#     dplyr::select(isotype, pct.swept, chrom) %>%
#     as.data.frame()
#   
#   sample.sizes[[i]] <- chrom.swept.list[[i]] %>%
#     dplyr::filter(swept.bin == "0.5-1") %>%
#     dplyr::full_join(., chrom.swept.list[[i]] %>%
#                        dplyr::filter(swept.bin == "0-0.20")) %>%
#     dplyr::group_by(swept.bin, chrom) %>%
#     dplyr::summarise(n = n())
#   
#   assign(paste("chr", unique(conservative.swept$chrom), "swept.list", sep = "."), 
#          paste(conservative.swept$isotype,sep = "",collapse = ","))
#   
#   assign(paste("chr", unique(conservative.unswept$chrom), "unswept.list", sep = "."), 
#          paste(conservative.unswept$isotype,sep = "",collapse = ","))
# }
# sample.sizes %>%
#   Reduce(rbind,.) %>%
#   dplyr::ungroup() %>%
#   tidyr::pivot_wider(names_from = chrom, values_from = n)


slow.strains <- strain.data %>%
  dplyr::select(strain, isotype) %>%
  dplyr::filter(strain %in% c("ECA246","QG2075","WN2002","ECA191","CX11254","LSJ1","QG536","NIC513","ECA363","ECA348")) %>%
  dplyr::select(isotype)

summer.prerequisite.strains <- isotypes.20210121 %>%
  dplyr::filter(isotype %in% c("N2", "MY16", "CB4856", # Janneke, Sam, and Tim DRC strains
                               "JU775", "CX11314", "DL238", # Janneke DRC strains
                               "ECA36","CB4855","ECA396","RC301","XZ1516")) # Sam and Tim DRC strains; FYI - ECA248 == CB4855

# # Subsampling of CeNDR Set at Various Sweptness and Population Size
# ####
# From Complete Set
pop.sizes <- c(96,144,192)
n.populations <- rep(50, length(pop.sizes))
generate.strain.set <- function(pop.size, n.populations){
  replicate.pops <- list()
  population.key <- list()
  for(i in 1:n.populations){
    
    mapping.strains <- isotypes.20210121 %>%
      dplyr::filter(!isotype %in% summer.prerequisite.strains$isotype,
                    !isotype %in% slow.strains$isotype)
      
     strains <- c(sample(unique(mapping.strains$isotype), pop.size-nrow(summer.prerequisite.strains), replace = F), summer.prerequisite.strains$isotype)
     strain.list <- paste(strains,
                          sep = "", collapse = ",")
     population.id <- ids::adjective_animal(n = 1, max_len = c(16,16)) %>%
       gsub(., pattern = "_", replacement = ".")
     replicate.pops[[i]] <- data.frame(population.id, strain.list)

     population.key[[i]] <- data.frame(strains) %>%
       dplyr::mutate(population.id = population.id) %>%
       `colnames<-`(c("isotype","population.id")) %>%
       dplyr::left_join(.,sweeps)

  }
  strain.lists <- Reduce(rbind,replicate.pops)
  population.keys <- Reduce(rbind,population.key)
  filename <- paste0("subsampled.strains_",today,"_n",pop.size,"_reps",n.populations)
  write.table(strain.lists,
              file = paste0("output/",filename,".txt"),
              quote = F, row.names = F, col.names = F)
}
purrr::map2(pop.sizes, n.populations, generate.strain.set)

# From Conservative Swept Set
pop.sizes <- c(144)
n.populations <- rep(50, length(pop.sizes))
generate.strain.set.swept <- function(pop.size, n.populations){
  replicate.pops <- list()
  population.key <- list()
  for(i in 1:n.populations){
    
    mapping.strains <- conservative.swept.isotypes %>%
      dplyr::filter(!isotype %in% summer.prerequisite.strains$isotype,
                    !isotype %in% slow.strains$isotype)
    
    strains <- c(sample(unique(mapping.strains$isotype), pop.size-nrow(summer.prerequisite.strains), replace = F), summer.prerequisite.strains$isotype)
    strain.list <- paste(strains,
                         sep = "", collapse = ",")
    population.id <- ids::adjective_animal(n = 1, max_len = c(16,16)) %>%
      gsub(., pattern = "_", replacement = ".")
    replicate.pops[[i]] <- data.frame(population.id, strain.list)

    population.key[[i]] <- data.frame(strains) %>%
      dplyr::mutate(population.id = population.id) %>%
      `colnames<-`(c("isotype","population.id")) %>%
      dplyr::left_join(.,sweeps)

  }
  strain.lists <- Reduce(rbind,replicate.pops)
  population.keys <- Reduce(rbind,population.key)
  filename <- paste0("subsampled.swept.strains_",today,"_n",pop.size,"_reps",n.populations)
  write.table(strain.lists,
              file = paste0("output/",filename,".txt"),
              quote = F, row.names = F, col.names = F)
  # write_tsv(population.keys,
  #           file = paste0("output/",filename,"population.sweep.key.tsv"),
  #           quote = F,col_names = F)

}
purrr::map2(pop.sizes, n.populations, generate.strain.set.swept)

# From Conservative Unswept Set
pop.sizes <- c(144)
n.populations <- rep(50, length(pop.sizes))
generate.strain.set.unswept <- function(pop.size, n.populations){
  replicate.pops <- list()
  population.key <- list()
  for(i in 1:n.populations){
    
    mapping.strains <- conservative.unswept.isotypes %>%
      dplyr::filter(!isotype %in% summer.prerequisite.strains$isotype,
                    !isotype %in% slow.strains$isotype)
    
    strains <- c(sample(unique(mapping.strains$isotype), pop.size-nrow(summer.prerequisite.strains), replace = F), summer.prerequisite.strains$isotype)
    strain.list <- paste(strains,
                         sep = "", collapse = ",")
    population.id <- ids::adjective_animal(n = 1, max_len = c(16,16)) %>%
      gsub(., pattern = "_", replacement = ".")
    replicate.pops[[i]] <- data.frame(population.id, strain.list)

    population.key[[i]] <- data.frame(strains) %>%
      dplyr::mutate(population.id = population.id) %>%
      `colnames<-`(c("isotype","population.id")) %>%
      dplyr::left_join(.,sweeps)

  }
  strain.lists <- Reduce(rbind,replicate.pops)
  population.keys <- Reduce(rbind,population.key)
  filename <- paste0("subsampled.unswept.strains_",today,"_n",pop.size,"_reps",n.populations)
  write.table(strain.lists,
              file = paste0("output/",filename,".txt"),
              quote = F, row.names = F, col.names = F)
  # write_tsv(population.keys,
  #           file = paste0("output/",filename,"population.sweep.key.tsv"),
  #           quote = F,col_names = F)

}
purrr::map2(pop.sizes, n.populations, generate.strain.set.unswept)

# From Liberal Swept Set
pop.sizes <- c(192)
n.populations <- rep(50, length(pop.sizes))
generate.strain.set.liberal.swept <- function(pop.size, n.populations){
  replicate.pops <- list()
  population.key <- list()
  for(i in 1:n.populations){
    
    mapping.strains <- liberal.swept.isotypes %>%
      dplyr::filter(!isotype %in% summer.prerequisite.strains$isotype,
                    !isotype %in% slow.strains$isotype)
    
    strains <- c(sample(unique(mapping.strains$isotype), pop.size-nrow(summer.prerequisite.strains), replace = F), summer.prerequisite.strains$isotype)
    strain.list <- paste(strains,
                         sep = "", collapse = ",")
    population.id <- ids::adjective_animal(n = 1, max_len = c(16,16)) %>%
      gsub(., pattern = "_", replacement = ".")
    replicate.pops[[i]] <- data.frame(population.id, strain.list)
    
    population.key[[i]] <- data.frame(strains) %>%
      dplyr::mutate(population.id = population.id) %>%
      `colnames<-`(c("isotype","population.id")) %>%
      dplyr::left_join(.,sweeps)
    
  }
  strain.lists <- Reduce(rbind,replicate.pops)
  population.keys <- Reduce(rbind,population.key)
  filename <- paste0("subsampled.liberal.swept.strains_",today,"_n",pop.size,"_reps",n.populations)
  write.table(strain.lists,
              file = paste0("output/",filename,".txt"),
              quote = F, row.names = F, col.names = F)
  # write_tsv(population.keys,
  #           file = paste0("output/",filename,"population.sweep.key.tsv"),
  #           quote = F,col_names = F)
  
}
purrr::map2(pop.sizes, n.populations, generate.strain.set.liberal.swept)

# From Liberal Unswept Set
pop.sizes <- c(144)
n.populations <- rep(50, length(pop.sizes))
generate.strain.set.liberal.unswept <- function(pop.size, n.populations){
  replicate.pops <- list()
  population.key <- list()
  for(i in 1:n.populations){
    
    mapping.strains <- liberal.unswept.isotypes %>%
      dplyr::filter(!isotype %in% summer.prerequisite.strains$isotype,
                    !isotype %in% slow.strains$isotype)
    
    strains <- c(sample(unique(mapping.strains$isotype), pop.size-nrow(summer.prerequisite.strains), replace = F), summer.prerequisite.strains$isotype)
    strain.list <- paste(strains,
                         sep = "", collapse = ",")
    population.id <- ids::adjective_animal(n = 1, max_len = c(16,16)) %>%
      gsub(., pattern = "_", replacement = ".")
    replicate.pops[[i]] <- data.frame(population.id, strain.list)
    
    population.key[[i]] <- data.frame(strains) %>%
      dplyr::mutate(population.id = population.id) %>%
      `colnames<-`(c("isotype","population.id")) %>%
      dplyr::left_join(.,sweeps)
    
  }
  strain.lists <- Reduce(rbind,replicate.pops)
  population.keys <- Reduce(rbind,population.key)
  filename <- paste0("subsampled.liberal.unswept.strains_",today,"_n",pop.size,"_reps",n.populations)
  write.table(strain.lists,
              file = paste0("output/",filename,".txt"),
              quote = F, row.names = F, col.names = F)
  # write_tsv(population.keys,
  #           file = paste0("output/",filename,"population.sweep.key.tsv"),
  #           quote = F,col_names = F)
  
}
purrr::map2(pop.sizes, n.populations, generate.strain.set.liberal.unswept)
####

# CeNDR Strain "Sets"
colnames(metadata)
cendr.strain.sets.nested <- strain.data %>%
  dplyr::select(isotype, strain_set) %>%
  dplyr::filter(strain_set %in% c(1:8)) %>%
  dplyr::arrange(strain_set) %>%
  dplyr::group_by(strain_set) %>%
  tidyr::nest()

divergent.set <-  strain.data %>%
  dplyr::select(isotype, strain_set) %>%
  dplyr::filter(strain_set == "D") %>%
  select(isotype)

strain.sets <- c(c(1:8),c(1:8))
strain.sets <- unique(strain.sets)
gather.strain.sets <- t(combn(strain.sets, 2)) %>%
  as.data.frame() %>%
  `colnames<-`(c("set1","set2"))

cendr.strain.set.list <- list()
for(i in 1:nrow(gather.strain.sets)){
  combo <- gather.strain.sets[i,]
  
  set1 <- as.numeric(combo$set1)
  set2 <- as.numeric(combo$set2)


  combined.strain.sets <- data.frame(cendr.strain.sets.nested$data[set1]) %>%
    dplyr::full_join(.,data.frame(cendr.strain.sets.nested$data[set2]))
  
  combo.name <- paste("cendr.set",set1,"set",set2,sep = ".")
  
  strain.list <- paste(combined.strain.sets$isotype,sep = "", collapse = ",")
  cendr.strain.set.list[[i]] <- data.frame(combo.name, strain.list)
}
cendr.strain.sets.df <- do.call(rbind, cendr.strain.set.list) %>%
  `colnames<-`(c("strain.set","strains"))


# Hahnel Benzimidazoles Strains
hahnel.isotypes <- list()
for(j in 1:nrow(hahnel_204)){
  hahnel.strain <- hahnel_204$strain[j]
  hahnel.isotypes[[j]] <- strain.data %>%
    dplyr::select(strain, isotype) %>%
    dplyr::filter(strain == hahnel.strain) %>%
    dplyr::select(isotype)
}
hahnel_isotypes <- Reduce(rbind, hahnel.isotypes) %>%
  dplyr::filter(!duplicated(isotype))
hahnel_isotypes <- paste(hahnel_isotypes$isotype, sep = "", collapse = ",")


# Combine strain sets
strain.lists <- data.frame(c("complete","complete.20210121","chr.I.swept","chr.I.unswept",
                             "chr.II.swept","chr.II.unswept",
                             "chr.III.swept","chr.III.unswept",
                             "chr.IV.swept","chr.IV.unswept",
                             "chr.V.swept","chr.V.unswept",
                             "chr.X.swept","chr.X.unswept","hahnel.isotypes"), 
                           c(complete.list, complete.20210121.list, chr.I.swept.list, chr.I.unswept.list,
                             chr.II.swept.list, chr.II.unswept.list,
                             chr.III.swept.list, chr.III.unswept.list,
                             chr.IV.swept.list, chr.IV.unswept.list,
                             chr.V.swept.list, chr.V.unswept.list,
                             chr.X.swept.list, chr.X.unswept.list, hahnel_isotypes)) %>%
  `colnames<-`(c("strain.set","strains")) %>%
  dplyr::full_join(cendr.strain.sets.df)
write.table(strain.lists, "output/swept.by.chr.strain.lists.txt", 
            quote = F, row.names = F, col.names = F)

