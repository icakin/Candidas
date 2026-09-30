# Extract isolate-level posterior draws for the translation analysis.
# Reads the brms fits and writes flat CSVs of the columns 31_translation.py needs.
# Run:  Rscript scripts/30_extract_draws.R
rds <- "results/rds"
extract <- function(file){
  f <- readRDS(file.path(rds,file)); sim <- attr(f$fit,"sim")
  ch <- sim$samples; warm <- sim$warmup; pn <- names(ch[[1]])
  M <- sapply(pn, function(p) unlist(lapply(ch, function(c) c[[p]][(warm+1):sim$iter])))
  as.data.frame(M, check.names=FALSE)
}
G <- extract("bayes_growth_ss.rds"); R <- extract("bayes_resp_arr.rds")
gcols <- grep("^b_E_Group|^b_Eh_Group|^b_Th_Group|^r_Isolate__E\\[|^r_Isolate__Th\\[", names(G), value=TRUE)
rcols <- grep("^b_E_Group|^r_Isolate__E\\[", names(R), value=TRUE)
dir.create("results/tables", showWarnings=FALSE, recursive=TRUE)
write.csv(G[,gcols], "results/tables/growth_draws_isolate.csv", row.names=FALSE)
write.csv(R[,rcols], "results/tables/resp_draws_isolate.csv", row.names=FALSE)
cat("wrote", nrow(G), "draws\n")
