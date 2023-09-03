library(dplyr)
library(data.table)
library(tidyr)

#read in lookup table 
lookup <- fread("GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt.gz", header = T)

#read in hapmap3 snps 
hapmap <- fread("../../eQTLsummary/GTEx_v8_genotype_EUR_HM3_exclude_dups.allchr.reorder.bim", header = F)

#filter lookup table 
filtered_lookup <- lookup %>%
	filter(rs_id_dbSNP151_GRCh38p7 %in% hapmap$V2)

# format b37 variant id to chromosome:position_allele1_allele2 such that alleles are in lexographical or indel length order 
filtered_lookup <- filtered_lookup %>%
  mutate(chromosome_position = sub("^(\\d+)_(\\d+)_.*", "\\1:\\2", variant_id_b37),
         alleles = sub("^\\d+_\\d+_(.*)", "\\1", variant_id_b37))

filtered_lookup <- filtered_lookup %>%
  separate(alleles, into = c("allele1", "allele2", "b"), sep = "_") %>%
  rowwise() %>%
  mutate(sorted_alleles = paste0(min(allele1, allele2), "_", max(allele1, allele2))) %>%
  ungroup() %>%
  mutate(variant_id_format = paste0(chromosome_position, "_", sorted_alleles))

filtered_lookup <- as.data.frame(filtered_lookup)

print(head(filtered_lookup))

#write lookup table 
write.table(filtered_lookup, file = "lookup_table_reorder_alleles_b38.txt", quote = F, row.names = F, col.names = T)


