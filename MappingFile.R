library(BaseSet)
library(mgsa)
library(tidyr)
library(dplyr)



##Taking the GAF file from
#https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/188/095/GCF_000188095.3_BIMP_2.2/

setwd("GSEA/")


addtofile <- getGAF("GCF_000188095.3_BIMP_2.2_gene_ontology.gaf.gz")
addtofile_df <- as.data.frame(addtofile)
colnames(addtofile_df)[2] <- c("GO")

term_def <- readGAF("GCF_000188095.3_BIMP_2.2_gene_ontology.gaf.gz")
term_def_df <- as.data.frame(term_def@setAnnotations)
term_def_df2 <- cbind(GO=rownames(term_def_df), term_def_df)
term_def_df2  <- term_def_df2[-1,]

combined_data <- merge(addtofile_df, term_def_df2, by = "GO", all.x = TRUE)
write.csv(combined_data, "combineGO.csv")
write.table(combined_data, "combineGO.txt", row.names=TRUE)


##Get File TopGO ready
gene_go_nounique <- addtofile_df %>%
  group_by(elements) %>%
  summarize(values = paste(GO, collapse = ","))

gene_go_uniquegos <- addtofile_df %>%
  group_by(elements) %>%
  summarize(values = paste(unique(GO), collapse = ","))


write.table(gene_go_nounique, "Gene2GoMap_someRepgos.txt", row.names=FALSE, quote=FALSE, col.names = FALSE)
write.table(gene_go_uniquegos, "Gene2GoMap_unique.txt", row.names=FALSE, quote=FALSE, col.names = FALSE)


#gene_go <- cbind(gene=combined_data$elements, GO=combined_data$GO)

#playing
#test_gg <- as.data.frame(gene_go[1:50,])


#collapsed_df2 <- test_gg %>%
#  group_by(gene) %>%
#  summarize(values = paste(unique(GO, collapse = ", ")))









