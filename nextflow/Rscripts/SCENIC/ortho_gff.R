args <- commandArgs(trailingOnly = TRUE)


if (file.exists(paste0(args[1], "/rtracklayer")) == F){
  dir.create(file.path(args[1]), showWarnings = FALSE)
  install.packages("BiocManager", lib = args[1])
  BiocManager::install("rtracklayer", lib = args[1])
}

library(tidyverse)
library(rtracklayer, lib = args[1])

ggid <- function(x){
  out <- x[[1]][grep("GeneID", x[[1]])] %>% 
    gsub("GeneID:", "", .)
  if (length(out) == 0){
    return("")
  } else {
    return(out)
  }
}

ortho_file <- read.table("all_genes.txt", quote = "",  stringsAsFactors = FALSE, sep = '\t')
colnames(ortho_file) <- c("orthodb_geneid", "orthodb_orgid", "protein_id", "synonyms", "uniprot_id",
                          "ensembl_ids", "NCBI_gid", "description", "coords", "DNA_id", "chromosome")



species_df <- read.csv(args[2], header = F)[,c(2,7)]  %>% 
	unique() %>% 
	rbind(c("Chicken", "9031_0"))


print(species_df)

colnames(species_df) <- c("species", "orthodb_orgid")

species_df <- filter(species_df,	species	%in% c(args[5],	'Chicken'))

ortho_file <- ortho_file %>% 
  merge(species_df)

bird_OGs <- read.table(args[3])

colnames(bird_OGs) <- c("Bird_cluster", "orthodb_geneid")
bird_mam_OGs <- read.table(args[4])
colnames(bird_mam_OGs) <- c("Vert_cluster", "orthodb_geneid")


ortho_file <- ortho_file %>% 
  merge(bird_OGs, all = T) %>% 
  merge(bird_mam_OGs, all = T)


gffs <- list.files(".", full.names = T, pattern = ".gtf")

print(gffs)
tmp_ofgff <- list()
for (i in 1:nrow(species_df)){
  print(i)
  gff <- gffs[grep(species_df[i,1], gffs)]
  print(gff)
  z <- rtracklayer::readGFF(gff) %>% as.data.frame()
  print(8)
  z2 <- z %>%  filter(type %in% c("CDS")) %>% 
    dplyr::select(gene, type, db_xref, protein_id) %>% 
    unique()
  print(8)
  print(head(z2))
  species_of <- ortho_file[ortho_file$species == species_df[i,1],]
  print(head(species_of))
  gene_ID_df = data.frame(gff_genes  = z2$gene,
                          protein_id = z2$protein_id, 
                          species = species_df[i,1], 
                          GeneID = lapply(z2$db_xref, ggid) %>% unlist()) %>% 
    merge(species_of, all.y = T, by = c('species', 'protein_id'))
  print(head(gene_ID_df))
  print(dim(species_of))
  print(dim(gene_ID_df))
  
  tmp_ofgff[[i]] <- gene_ID_df
}

OF_GFF <- do.call(rbind, tmp_ofgff) %>% 
  bind_rows(ortho_file[ortho_file$species %in% c("Human", "Mouse"),])
print(head(OF_GFF))
write.table(OF_GFF, file = "OrthoDB_GFF_orthologs.tsv", sep = '\t')
