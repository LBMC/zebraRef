sapply(c("biomaRt", "Biobase", "GEOquery", "limma", "utils", "RAPToR"), 
       requireNamespace, quietly = T)

# gene ids
load("data/Dre_genes.RData")

# fetch data
geo_id <- "GSE24616"
geo_obj <- GEOquery::getGEO(geo_id)[[1]]

probe_ids <- Biobase::pData(Biobase::featureData(geo_obj))

pdat <- Biobase::pData(geo_obj)
pdat <- pdat[, c("title", "geo_accession", "source_name_ch1",
                 "developmental stage:ch1", "developmental timing:ch1", 
                 "number of individuals per sample:ch1")]
colnames(pdat)[4:6] <- c("devstage", "age_char", "nb_indiv")

time_values <- data.frame(age_char = unique(pdat$age_char), 
                          age_min = c(0, 15, 45, 75, 105, 135, 165, 200, 240,
                                      280, 320, 360, 420, 480, 540, 600, 620, 660,
                                      700, 720, 780, 840, 900, 960, 1020, 1080, 1140, 
                                      1200, 1260, 1320, 1380, 1500, 1620, 1800, 2040, 2280,
                                      2520, 2880, 3600, 4320, 5760, 8640, 11520, 14400, 20160,
                                      25920, 34560, 43200, 57600, 64800, 79200, 93600, 115200, 129600,
                                      151200, 172800, 302400, 388800, 604800, 777600, 907200),
                          stringsAsFactors = F)
pdat$age_min <- time_values[match(pdat$age_char, time_values$age_char), "age_min"]
pdat$age_hr <- pdat$age_min/60
pdat$nb_indiv <- as.numeric(pdat$nb_indiv)
pdat$devstage <- as.factor(pdat$devstage)

pdat$title <- as.character(pdat$title)
pdat$sex <- as.factor(gsub("[^_]*_[^_]*_([^_]*)_rep(\\d)", "\\1", pdat$title))
pdat$rep <- as.factor(gsub("[^_]*_[^_]*_([^_]*)_rep(\\d)", "\\2", pdat$title))


gdat <- Biobase::exprs(geo_obj)
rownames(gdat) <- probe_ids$ID


gdat <- RAPToR::format_ids(gdat, probe_ids, from = "ID", to = "ENSEMBL_ID")
gdat <- RAPToR::format_ids(gdat, Dre_genes, from = "transcript_id", to = "ens_id")

colnames(gdat) <- pdat$title

str(pdat)
# keep non-adult, mixed samples for reference
p <- pdat[pdat$sex == "mixed" & pdat$age_hr < 1e4, ]

w_s <- match(unique(p$devstage), p$devstage)
w_s <- w_s[-length(w_s)]
w_e <-c(w_s[-1], nrow(p))
# plot(p$age_hr+1, 1:103, log = "x")
# abline(v=1+p$age_hr[w_s])

# write file with devstages
utils::write.csv(data.frame(
  name = as.character(p$devstage[w_s]), 
  tstart = round(p$age_hr[w_s], 2), tend = round(p$age_hr[w_e],2),
  tunit = rep("hours post-fertilization", length(w_s)),
  stringsAsFactors = F),
  file = "data-raw/devstages.csv", row.names = F, quote = F)

p <- p[, c("title", "geo_accession", "age_hr")]
colnames(p) <- c("sname", "accession", "age")
p$age_ini <- p$age

gdat <- log1p(limma::normalizeBetweenArrays(gdat[,p$sname], method = "quantile"))

pca <- summary(stats::prcomp(t(gdat), rank = 40, center = TRUE, scale = FALSE))
nc <- sum(pca$importance[3,] < .95) + 1



Dre_emb_larv <- list(g = gdat, p = p, 
                     geim_params = list(
                       formula = "X~s(age, bs='cr')",
                       method = "gam",
                       dim_red = "pca",
                       nc = nc),
                     t.unit = "h post fertilization",
                     cov.levels = NULL,
                     metadata = list("organism" = "D. rerio",
                                     "profiling" = "whole-organism, bulk",
                                     "technology" = "Microarray")
)

save("Dre_emb_larv", file = "data/Dre_emb_larv.RData", compress = "xz")
rm(pdat, gdat, p, 
   geo_id, geo_obj, probe_ids, time_values, 
   Dre_genes, Dre_emb_larv,
   w_s, w_e, pca, nc)
