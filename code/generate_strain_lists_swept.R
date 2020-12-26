library(tidyverse)
require(ids)
# setwd("~/Documents/projects/NemaScan_Performance/")
setwd("~/Documents/AndersenLab/NemaScan_Performance/")
strain.data <- data.table::fread("data/CelegansStrainData.tsv") 
sweeps <- data.table::fread("data/sweep_summary.tsv")
hahnel_204 <- data.table::fread("data/hahnel_209.tsv") %>%
  dplyr::select(strain)
today <- format(Sys.time(), '%Y%m%d')
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

# Sweep-based strain lists
chr.I.swept <- metadata %>%
  dplyr::filter(I == "Swept") %>%
  dplyr::select(isotype) %>%
  dplyr::filter(!duplicated(isotype))
chr.I.unswept <- metadata %>%
  dplyr::filter(I == "Unswept") %>%
  dplyr::select(isotype)%>%
  dplyr::filter(!duplicated(isotype))
chr.I.swept.list <- paste(chr.I.swept$isotype,sep = "",collapse = ",")
chr.I.unswept.list <- paste(chr.I.unswept$isotype,sep = "",collapse = ",")

chr.II.swept <- metadata %>%
  dplyr::filter(II == "Swept") %>%
  dplyr::select(isotype)%>%
  dplyr::filter(!duplicated(isotype))
chr.II.unswept <- metadata %>%
  dplyr::filter(II == "Unswept") %>%
  dplyr::select(isotype)%>%
  dplyr::filter(!duplicated(isotype))
chr.II.swept.list <- paste(chr.II.swept$isotype,sep = "",collapse = ",")
chr.II.unswept.list <- paste(chr.II.unswept$isotype,sep = "",collapse = ",")

chr.III.swept <- metadata %>%
  dplyr::filter(III == "Swept") %>%
  dplyr::select(isotype)%>%
  dplyr::filter(!duplicated(isotype))
chr.III.unswept <- metadata %>%
  dplyr::filter(III == "Unswept") %>%
  dplyr::select(isotype)%>%
  dplyr::filter(!duplicated(isotype))
chr.III.swept.list <- paste(chr.III.swept$isotype,sep = "",collapse = ",")
chr.III.unswept.list <- paste(chr.III.unswept$isotype,sep = "",collapse = ",")

chr.IV.swept <- metadata %>%
  dplyr::filter(IV == "Swept") %>%
  dplyr::select(isotype)%>%
  dplyr::filter(!duplicated(isotype))
chr.IV.unswept <- metadata %>%
  dplyr::filter(IV == "Unswept") %>%
  dplyr::select(isotype)%>%
  dplyr::filter(!duplicated(isotype))
chr.IV.swept.list <- paste(chr.IV.swept$isotype,sep = "",collapse = ",")
chr.IV.unswept.list <- paste(chr.IV.unswept$isotype,sep = "",collapse = ",")

chr.V.swept <- metadata %>%
  dplyr::filter(V == "Swept") %>%
  dplyr::select(isotype)%>%
  dplyr::filter(!duplicated(isotype))
chr.V.unswept <- metadata %>%
  dplyr::filter(V == "Unswept") %>%
  dplyr::select(isotype)%>%
  dplyr::filter(!duplicated(isotype))
chr.V.swept.list <- paste(chr.V.swept$isotype,sep = "",collapse = ",")
chr.V.unswept.list <- paste(chr.V.unswept$isotype,sep = "",collapse = ",")

chr.X.swept <- metadata %>%
  dplyr::filter(X == "Swept") %>%
  dplyr::select(isotype)%>%
  dplyr::filter(!duplicated(isotype))
chr.X.unswept <- metadata %>%
  dplyr::filter(X == "Unswept") %>%
  dplyr::select(isotype)%>%
  dplyr::filter(!duplicated(isotype))
chr.X.swept.list <- paste(chr.X.swept$isotype,sep = "",collapse = ",")
chr.X.unswept.list <- paste(chr.X.unswept$isotype,sep = "",collapse = ",")

genome.swept <- metadata %>%
  dplyr::filter(swept_chroms == 4) %>%
  dplyr::select(isotype)%>%
  dplyr::filter(!duplicated(isotype))
genome.unswept <- metadata %>%
  dplyr::filter(swept_chroms < 4) %>%
  dplyr::select(isotype)%>%
  dplyr::filter(!duplicated(isotype))
genome.swept.list <- paste(genome.swept$isotype,sep = "",collapse = ",")
genome.unswept.list <- paste(genome.unswept$isotype,sep = "",collapse = ",")

# Subsampling of Full CeNDR Set
full.subsample.96.1 <- paste(c(sample(genome.swept$isotype, 41, replace = F),
                                sample(genome.unswept$isotype, 55, replace = F)),
                              sep = "", collapse = ",")
full.subsample.96.2 <- paste(c(sample(genome.swept$isotype, 41, replace = F),
                               sample(genome.unswept$isotype, 55, replace = F)),
                             sep = "", collapse = ",")
full.subsample.96.3 <- paste(c(sample(genome.swept$isotype, 41, replace = F),
                               sample(genome.unswept$isotype, 55, replace = F)),
                             sep = "", collapse = ",")

full.subsample.200.1 <- paste(c(sample(genome.swept$isotype, 86, replace = F),
                                sample(genome.unswept$isotype, 114, replace = F)),
                              sep = "", collapse = ",")
full.subsample.300.1 <- paste(c(sample(genome.swept$isotype, 129, replace = F),
                                sample(genome.unswept$isotype, 171, replace = F)),
                              sep = "", collapse = ",")
full.subsample.200.2 <- paste(c(sample(genome.swept$isotype, 86, replace = F),
                                sample(genome.unswept$isotype, 114, replace = F)),
                              sep = "", collapse = ",")
full.subsample.300.2 <- paste(c(sample(genome.swept$isotype, 129, replace = F),
                                sample(genome.unswept$isotype, 171, replace = F)),
                              sep = "", collapse = ",")
full.subsample.200.3 <- paste(c(sample(genome.swept$isotype, 86, replace = F),
                                sample(genome.unswept$isotype, 114, replace = F)),
                              sep = "", collapse = ",")
full.subsample.300.3 <- paste(c(sample(genome.swept$isotype, 129, replace = F),
                                sample(genome.unswept$isotype, 171, replace = F)),
                              sep = "", collapse = ",")
full.subsample.200.4 <- paste(c(sample(genome.swept$isotype, 86, replace = F),
                                sample(genome.unswept$isotype, 114, replace = F)),
                              sep = "", collapse = ",")
full.subsample.300.4 <- paste(c(sample(genome.swept$isotype, 129, replace = F),
                                sample(genome.unswept$isotype, 171, replace = F)),
                              sep = "", collapse = ",")
full.subsample.200.5 <- paste(c(sample(genome.swept$isotype, 86, replace = F),
                                sample(genome.unswept$isotype, 114, replace = F)),
                              sep = "", collapse = ",")
full.subsample.300.5 <- paste(c(sample(genome.swept$isotype, 129, replace = F),
                                sample(genome.unswept$isotype, 171, replace = F)),
                              sep = "", collapse = ",")


# Subsampling of CeNDR Set at Various Sweptness and Population Size
# From Complete Set
pop.sizes <- c(seq(1:8)*48)
n.populations <- rep(200, length(pop.sizes))
generate.strain.set <- function(pop.size, n.populations){
  replicate.pops <- list()
  population.key <- list()
  for(i in 1:n.populations){
     strains <- sample(unique(metadata$isotype), pop.size, replace = F)
     strain.list <- paste(strains, 
                          sep = "", collapse = ",")
     population.id <- ids::adjective_animal(n = 1, max_len = c(10,10)) %>%
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

# From Swept Set
pop.sizes <- c(seq(1:7)*24)
n.populations <- rep(200, length(pop.sizes))
generate.strain.set.swept <- function(pop.size, n.populations){
  replicate.pops <- list()
  population.key <- list()
  for(i in 1:n.populations){
    strains <- sample(genome.swept$isotype, pop.size, replace = F)
    strain.list <- paste(strains, 
                         sep = "", collapse = ",")
    population.id <- ids::adjective_animal(n = 1, max_len = c(10,10)) %>%
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

# From Unwept Set
pop.sizes <- c(seq(1:9)*24)
n.populations <- rep(200, length(pop.sizes))
generate.strain.set.unswept <- function(pop.size, n.populations){
  replicate.pops <- list()
  population.key <- list()
  for(i in 1:n.populations){
    strains <- sample(genome.unswept$isotype, pop.size, replace = F)
    strain.list <- paste(strains, 
                         sep = "", collapse = ",")
    population.id <- ids::adjective_animal(n = 1, max_len = c(10,10)) %>%
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
strain.lists <- data.frame(c("complete","chr.I.swept","chr.I.unswept",
                             "chr.II.swept","chr.II.unswept",
                             "chr.III.swept","chr.III.unswept",
                             "chr.IV.swept","chr.IV.unswept",
                             "chr.V.swept","chr.V.unswept",
                             "chr.X.swept","chr.X.unswept",
                             "genome.swept","genome.unswept",
                             "full.subsample.200.1", "full.subsample.200.2", "full.subsample.200.3", "full.subsample.200.4", "full.subsample.200.5",
                             "full.subsample.300.1", "full.subsample.300.2", "full.subsample.300.3", "full.subsample.300.4", "full.subsample.300.5",
                             "full.subsample.96.1", "full.subsample.96.2", "full.subsample.96.3","hahnel.isotypes"), 
                           c(complete.list, chr.I.swept.list, chr.I.unswept.list,
                             chr.II.swept.list, chr.II.unswept.list,
                             chr.III.swept.list, chr.III.unswept.list,
                             chr.IV.swept.list, chr.IV.unswept.list,
                             chr.V.swept.list, chr.V.unswept.list,
                             chr.X.swept.list, chr.X.unswept.list,
                             genome.swept.list, genome.unswept.list,
                             full.subsample.200.1, full.subsample.200.2, full.subsample.200.3, full.subsample.200.4, full.subsample.200.5,
                             full.subsample.300.1, full.subsample.300.2, full.subsample.300.3, full.subsample.300.4, full.subsample.300.5,
                             full.subsample.96.1, full.subsample.96.2, full.subsample.96.3, hahnel_isotypes)) %>%
  `colnames<-`(c("strain.set","strains")) %>%
  dplyr::full_join(cendr.strain.sets.df)
write.table(strain.lists, "output/swept.by.chr.strain.lists.txt", 
            quote = F, row.names = F, col.names = F)

