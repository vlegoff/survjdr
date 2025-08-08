library(survival)
library(data.table)

library(mlr3)
library(mlr3proba)
library(paradox)
load("BLCA.RData")
x <- fread("results_glbl_mirna_cnv_rna_mut_BLCA.csv",header=TRUE)
omiques = c("mutation", "rna", "cnv", "mirna_geom_mean")
long <- melt(x,
             id.vars = c("bcr_patient_barcode", "ensembl_gene_id"),
             measure.vars = omiques,
             variable.name = "omique",
             value.name = "value")

long[, gene_omique := paste0(ensembl_gene_id, "_", omique)]
wide <- dcast(long, bcr_patient_barcode ~ gene_omique, value.var = "value")
mat <- as.matrix(wide[, -1, with = FALSE])
rownames(mat) <- wide$bcr_patient_barcode
clinical <- dat[match(rownames(mat), dat$bcr_patient_barcode), c("time", "status")]
mat_with_surv <- cbind(mat, clinical)

mirna_colonnes <- gsub("_mirna_geom_mean$","_mirna",colnames(mat_with_surv))
setnames(mat_with_surv, old = colnames(mat_with_surv), new = mirna_colonnes)

task = as_task_surv(as.data.frame(mat_with_surv), event="status", time="time", id="BLCA")
learner <- LearnertFOBI$new()
learner$param_set$values <- list(
  blocks = c("rna", "mirna", "mutation","cnv"),
  clinical_fav = FALSE,
  CV_measure = "V&VH",
  seed = 123,
  rankPCA = 5
)

learner$train(task)


