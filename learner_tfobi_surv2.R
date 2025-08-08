#' @export
LearnertFOBI = R6::R6Class("LearnertFOBI",
                           inherit = LearnerSeqMod,
                           public = list(
                             #' @description
                             #' Creates a new instance of this [R6][R6::R6Class] class.
                             initialize = function() {
                               param_set = ps(
                                 blocks = p_uty(tags = c("train", "predict")),
                                 clinical_fav = p_lgl(default = TRUE, tags = c("train", "predict")),
                                 nfolds = p_int(lower = 1L, default = 10L, tags = c("train")),
                                 rankPCA = p_int(lower = 1L, default = 5L, tags = c("train", "predict")),
                                 nlambdas = p_int(lower = 10L, upper = 1000L, default = 100L, tags = c("train")),
                                 CV_measure = p_fct(
                                   levels = c("cindex", "ibs", "ibsRR", "auc", "C", "deviance", "basic", "V&VH", "linpred"),
                                   default = "cindex",
                                   tags = c("train")
                                 ),
                                 cv_save_path = p_uty(default = NULL, tags = c("predict")),
                                 seed = p_int(lower = 0L, special_vals = list(NULL), default = NULL, tags = c("train"))
                               )
                               
                               super$initialize(
                                 id = "surv.tfobi",
                                 packages = c("tensorBSS", "mlr3misc"),
                                 feature_types = c("integer", "numeric", "factor"),
                                 predict_types = c("crank", "lp", "distr"),
                                 param_set = param_set,
                                 properties = c(),
                                 man = "",
                                 label = ""
                               )
                             }
                           ),
                           
                           private = list(
                             
                             tfobi_args = NULL,
                             
                             train_jdr = function(x, y, pars) {
                               library(abind)
                               rna_t     <- t(x$rna)     
                               mirna_t   <- t(x$mirna)   
                               mutation_t<- t(x$mutation)
                               cnv_t     <- t(x$cnv)     
                               
                               tensor_array <- abind(
                                 rna_t,
                                 mirna_t,
                                 mutation_t,
                                 cnv_t,
                                 along = 2
                               )
                               
                               
                               x_centered <- tensorBSS::tensorCentering(tensor_array)
                               tensor_mean <- attr(x_centered, "location")
                               gene_vars <- apply(x_centered, 1, function(x) var(as.vector(x)))
                               top_genes <- order(gene_vars, decreasing = TRUE)[1:1200]
                               x_centered <- x_centered[top_genes,,]
                               
                               tpca_fit = tensorBSS::tPCA(
                                 x = x_centered,
            
                                 d = pars$rankPCA
                               )
                               
                               x_tpca = tpca_fit$S
                               
                               if (is.null(private$tfobi_args)) {
                                 private$tfobi_args = list(
                                   center = TRUE,
                                   whiten = TRUE,
                                   n.comp = pars$rankJ,
                                   norm = TRUE
                                 )
                               }
                               tfobi_fit = mlr3misc::invoke(
                                 tensorBSS::tFOBI,
                                 x = x_tpca,
                                 .args = private$tfobi_args
                               )
                               
                               dims <- dim(tfobi_fit$S)
                               permuted <- aperm(tfobi_fit$S, c(3, 1, 2))
                               matrix(permuted, nrow = dims[3], ncol = dims[1] * dims[2])
                               
                               return(list(
                                 mean = tensor_mean,
                                 tpca = tpca_fit,
                                 tfobi = tfobi_fit
                               ))
                             },
                             
                             predict_jdr = function(new_x, jdr_model, pars) {
                               library(abind)
                               rna_t     <- t(new_x$rna)    
                               mirna_t   <- t(new_x$mirna)  
                               mutation_t<- t(new_x$mutation)
                               cnv_t     <- t(new_x$cnv)     
                               
                             
                               tensor_array <- abind(
                                 rna_t,
                                 mirna_t,
                                 mutation_t,
                                 cnv_t,
                                 along = 2
                               )
                               
                               x_centered <- sweep(tensor_array, MARGIN = c(1, 2), STATS = jdr_model$mean, FUN = "-")
                               gene_vars <- apply(x_centered, 1, function(x) var(as.vector(x)))
                               top_genes <- order(gene_vars, decreasing = TRUE)[1:1200]
                               x_centered <- x_centered[top_genes,,]
                               for (m in seq_along(jdr_model$tpca$W)) {
                                 x_centered <- tensorBSS::tensorTransform(x_centered, t(jdr_model$tpca$U[[m]]), m)
                               }
                               
                               for (m in seq_along(jdr_model$tfobi$W)) {
                                 x_centered <- tensorBSS::tensorTransform(x_centered, jdr_model$tfobi$W[[m]], m)
                               }
                               
                               dims <- dim(x_centered)
                               permuted <- aperm(x_centered, c(3, 1, 2))
                               matrix(permuted, nrow = dims[3], ncol = dims[1] * dims[2])
                             }
                           ))
