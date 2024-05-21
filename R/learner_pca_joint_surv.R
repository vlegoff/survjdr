#' @export
LearnerSurvPCAJ = R6::R6Class("LearnerSurvPCAJ",
  inherit = LearnerSeqMod,
  public = list(
    #' @description
    #' Creates a new instance of this [R6][R6::R6Class] class.
    initialize = function() {
      param_set = ps(
        blocks=p_uty(tags=c("train", "predict")),
        clinical_fav=p_lgl(default=TRUE, tags=c("train", "predict")),
        rankJ=p_int(1L, default=1L, tags=c("train", "predict")),
        nfolds=p_int(1L, default=10L, tags=c("train")),
        nlambdas=p_int(10L, 1000L, default=100L, tags=c("train")),
        CV_measure=p_fct(c("cindex", "ibs", "ibsRR",  "auc", "C", "deviance",
                           "basic", "V&VH", "linpred"), default="cindex",
                          tags=c("train")),
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
              nf=pars$rankJ
            )
        }

        x_joint = Reduce(cbind, x)
        pca = mlr3misc::invoke(ade4::dudi.pca,
            .args=c(list(df=x_joint), private$pca_args))
        return(list(x=as.matrix(pca$li), jdr=pca))
    },

    predict_jdr = function(newx, jdr, pars) {
        newx_joint = Reduce(cbind, newx)
        newpca = ade4::suprow(jdr, newx_joint)
        return(list(x=as.matrix(newpca$li)))
    }

  )
)

.extralrns_dict$add("surv.pca_joint", LearnerSurvPCAJ)
