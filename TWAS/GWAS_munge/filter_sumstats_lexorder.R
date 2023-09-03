library(dplyr)
library(data.table)

#read in lookup table
lookup <- fread("lookup_table_reorder_alleles.txt", header = T)

args <- commandArgs(trailingOnly = TRUE)
sumstat_file <- args[1]
out <- args[2]

#read in sumstat
sumstat <- fread(sumstat_file, header = T)

if (ncol(sumstat) > 17){
	sumstat <- sumstat[, c(1:17)]
}

print(head(sumstat))

lookup_table <- lookup[, c(7,14)]
colnames(lookup_table) <- c("rsid", "rs_number")

result <- sumstat %>%
  left_join(lookup_table, by = "rs_number") %>%
  filter(!is.na(rsid))

result[,1] <- result$rsid
result <- result[, -18]

write.table(result, file = out, quote = F, row.names = F, col.names = T)

