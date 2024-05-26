#' @export
LearnerSurvPCAI = R6::R6Class("LearnerSurvPCAI",
  inherit = LearnerSeqMod,
  public = list(
    #' @description
    #' Creates a new instance of this [R6][R6::R6Class] class.
    initialize = function() {
      param_set = ps(
        blocks=p_uty(tags=c("train", "predict")),
        clinical_fav=p_lgl(default=TRUE, tags=c("train", "predict")),
        rankA=p_int(1L, default=1L, tags=c("train", "predict")),
        nfolds=p_int(1L, default=10L, tags=c("train")),
        nlambdas=p_int(10L, 1000L, default=100L, tags=c("train")),
        CV_measure=p_fct(c("cindex", "ibs", "ibsRR",  "auc", "C", "deviance",
                           "basic", "V&VH", "linpred"), default="cindex",
                          tags=c("train")),
        cv_save_path=p_uty(default=NULL, tags=c("predict")),
        seed=p_int(0L, special_vals=list(NULL), default=NULL, tags=c("train"))
      )
      param_set$values = param_set$default

      super$initialize(
        id = "surv.pca_joint",
        packages = c("ade4", "mlr3misc"),
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

    pca_args = NULL,

    train_jdr = function(x, y, pars) {
        
        if(is.null(private$pca_args)) {
            private$pca_args = list(
              center=TRUE,
              scale=TRUE,
              scannf=FALSE,
              nf=pars$rankA
            )
        }

      pca = lapply(x,
        function(block) {
          mlr3misc::invoke(ade4::dudi.pca,
            .args=c(list(df=block), private$pca_args))
        }
      )

      latent_space = Reduce(cbind, lapply(pca, `[[`, "li"))
      latent_space = as.matrix(latent_space)
      colnames(latent_space) = paste(
          rep(pars$blocks, each=private$pca_args$nf),
          rep(1:private$pca_args$nf, private$pca_args$nf),
          sep="."
      )

      return(list(x=as.matrix(latent_space), jdr=pca))
    },

    predict_jdr = function(newx, jdr, pars) {

      pca = lapply(seq_along(newx),
        function(i) ade4::suprow(jdr[[i]], newx[[i]]))
      latent_space_new = Reduce(cbind, lapply(pca, `[[`, "lisup"))
      latent_space_new = as.matrix(latent_space_new)
      colnames(latent_space_new) = paste(
          rep(pars$blocks, each=private$pca_args$nf),
          rep(1:private$pca_args$nf, private$pca_args$nf),
          sep="."
      )

      return(list(x=as.matrix(latent_space_new)))
    }

  )
)

.extralrns_dict$add("surv.pca_indiv", LearnerSurvPCAI)
