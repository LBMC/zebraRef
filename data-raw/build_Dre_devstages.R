datasets_s <- "name, tstart, tend, tunit
Dre_emb_larval, 0, 1080,hours post-fertilization
"


devst <- read.table("data-raw/devstages.csv", h = T, sep = ',', stringsAsFactors = F)
devst$tunit <- factor(devst$tunit)

dss <- read.table(text = datasets_s, h = T, sep = ',', stringsAsFactors = F)
dss$tunit <- factor(dss$tunit)

Dre_devstages <- list(
  devstages = devst, 
  datasets = dss)

save("Dre_devstages", file = "data/Dre_devstages.RData")

rm(devst,
   datasets_s, dss,
   Dre_devstages)