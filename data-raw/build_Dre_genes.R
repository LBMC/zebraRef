requireNamespace('biomaRt', quietly = T)

# get gene ids and length info from the biomart
mart <- biomaRt::useMart("ensembl", dataset = "drerio_gene_ensembl")
Dre_genes <- biomaRt::getBM(attributes = c("ensembl_gene_id", 
                                           "ensembl_transcript_id",
                                           "external_gene_name",
                                           "transcript_length"),
                            mart = mart)
colnames(Dre_genes)[1:3] <- c("ens_id", "transcript_id", "public_name")

# save object to data
save('Dre_genes', file = "data/Dre_genes.RData", compress = "xz")
rm(mart, Dre_genes)