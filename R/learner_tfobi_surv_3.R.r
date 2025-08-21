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
                               
                               n_patients <- nrow(x[[1]])
                               n_genes    <- ncol(x[[1]])
                               n_blocks   <- length(x)
                               
                               # Construction du tenseur (gènes × blocs × patients)
                               tensor_array <- array(NA, dim = c(n_genes, n_blocks, n_patients))
                               
                               for (i in seq_along(x)) {
                                 tensor_array[, i, ] <- t(x[[i]]) 
                               }
                               # Centrage global par tensorBSS
                               x_centered <- tensorBSS::tensorCentering(tensor_array)
                               tensor_mean <- attr(x_centered, "location")
                               
                               # Centrage et normalisation gène par gène
                               gene_sds   <- apply(x_centered, 1, sd)
                              
                               x_centered <- sweep(x_centered, 1, gene_means, "-")
                               x_centered <- sweep(x_centered, 1, gene_sds, "/")
                               
                               # PCA tensorielle
                               tpca_fit = tensorBSS::tPCA(
                                 x = x_centered,

                                 d = c(pars$rankPCA, n_blocks, n_patients) # dimensions après réduction
                               )
                               # tenseur transformé
                               x_tpca = tpca_fit$S

                               # FOBI tensoriel
                               tfobi_fit = mlr3misc::invoke(
                                 tensorBSS::tFOBI,
                                 x = x_tpca
                               )
                               # Mode 3 unfold en forme patients × (gènes*blocs)
                               dims <- dim(tfobi_fit$S)
                               permuted <- aperm(tfobi_fit$S, c(3, 1, 2))
                               latent <- matrix(permuted, nrow = dims[3], ncol = dims[1] * dims[2])
                               
                               return(list(
                                 x= latent,
                                 jdr = list(
                                   mean = tensor_mean,
                                   sd = gene_sds,
                                   tpca = tpca_fit,
                                   tfobi = tfobi_fit
                                 )
                                 
                               ))
                             },
                             
                             predict_jdr = function(new_x, jdr_model, pars) {
                               
                               n_patients <- nrow(new_x[[1]])
                               n_genes    <- ncol(new_x[[1]])
                               n_blocks   <- length(new_x)
                               
                               tensor_array <- array(NA, dim = c(n_genes, n_blocks, n_patients))
                               
                               for (i in seq_along(new_x)) {
                                 tensor_array[, i, ] <- t(new_x[[i]]) 
                               }
                               
                               # Appliquer les mêmes normalisations que dans train
                               x_centered <- sweep(tensor_array, MARGIN = c(1, 2), STATS = jdr_model$mean, FUN = "-")
                               
                               x_centered <- sweep(x_centered, 1, jdr_model$sd , "/")
                              
                               # Reprojection avec tPCA
                               for (m in seq_along(jdr_model$tpca$U)) {
                                 x_centered <- tensorBSS::tensorTransform(x_centered, t(jdr_model$tpca$U[[m]]), m)
                               }
                               
                               # Reprojection avec tFOBI
                               for (m in seq_along(jdr_model$tfobi$W)) {
                                 x_centered <- tensorBSS::tensorTransform(x_centered, jdr_model$tfobi$W[[m]], m)
                               }
                               
                               dims <- dim(x_centered)
                               permuted <- aperm(x_centered, c(3, 1, 2))
                               latent <- matrix(permuted, nrow = dims[3], ncol = dims[1] * dims[2])
                               
                               return (list(x=latent))
                             }
                           ))
