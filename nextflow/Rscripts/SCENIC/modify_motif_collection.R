args <- commandArgs(trailingOnly = TRUE)

library(tidyverse)

OF_GFF <- read.table(args[1], sep = '\t')
sp <- args[2]

motif_data <- read.table("v10nr_clust_public/snapshots/motifs-v10-nr.chicken-m0.00001-o0.0.tbl", sep = '\t', 
                   header = T, comment.char="")

unique(motif_data$gene_name)  %>% length()

OF_GFF_chicken <- filter(OF_GFF, species == "Chicken" & 
                           (gff_genes %in% motif_data$gene_name)) %>% 
  filter(!is.na(Bird_cluster)) %>% 
    dplyr::select(gff_genes, Bird_cluster, GeneID, species) %>% 
    unique()
print(head(OF_GFF_chicken))
unique(OF_GFF_chicken$gff_genes)  %>% length()

species_motif_data <- list()



OF_GFF_sp <- filter(OF_GFF, species == sp,
                      Bird_cluster %in% OF_GFF_chicken$Bird_cluster) %>% 
    dplyr::select(gff_genes, Bird_cluster, GeneID, species) %>% 
    unique()
print('here')
tmp1 <- merge(OF_GFF_sp, OF_GFF_chicken, by = 'Bird_cluster', suffixes = c("_tmp", "_Chicken"), all = F)
dup1 <- tmp1$GeneID_Chicken[duplicated(tmp1$GeneID_Chicken)]
dup2 <- tmp1$GeneID_tmp[duplicated(tmp1$GeneID_tmp)]
tmp2 <- filter(tmp1, GeneID_Chicken %in% dup1 | 
                  GeneID_tmp %in% dup2) %>% 
    filter(gff_genes_tmp == gff_genes_Chicken)

mgd <- rbind(tmp2, filter(tmp1, !(GeneID_Chicken %in% dup1) & 
                             !(GeneID_tmp %in% dup2)))
  
genes <- mgd[,"gff_genes_tmp"]
names(genes) <- mgd$gff_genes_Chicken
  
species_motif_data[[sp]] <- motif_data %>% 
    mutate(gene_name = genes[gene_name]) %>% 
    filter(!is.na(gene_name))
print(length(unique(species_motif_data[[sp]]$gene_name)))
colnames(species_motif_data[[sp]])[1] <- '#motif_id'
write.table(species_motif_data[[sp]], file = paste0(sp, "_motifs.tbl"), sep = '\t', col.names = T, quote = F, row.names = F)
  

