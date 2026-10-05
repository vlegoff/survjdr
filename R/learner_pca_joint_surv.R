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
        center=p_lgl(default=TRUE, tags=c("train", "predict")),
        scale=p_lgl(default=FALSE, tags=c("train", "predict")),
        norm=p_lgl(default=TRUE, tags=c("train", "predict")),
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
              center=FALSE,
              scale=FALSE,
              scannf=FALSE,
              nf=pars$rankJ
            )
        }

        n = names(x)
        cols = lapply(x, function(xi) apply(xi, 2, stats::var)!=0)
        x = lapply(names(x), function(n) x[[n]][,cols[[n]]])

        if (pars$center) {
            means = lapply(x, \(block) apply(block, 2, mean))
            x = lapply(seq_along(x),
                function(i) scale(x[[i]], center=means[[i]], scale=FALSE))
        } else {
            means = NULL
        }

        if (pars$scale) {
            standevs = lapply(x, \(block) apply(block, 2, sd))
            x = lapply(seq_along(x),
                function(i) scale(x[[i]], center=FALSE, scale=standevs[[i]]))
        } else {
            standevs = NULL
        }

        if (pars$norm) {
            norms = lapply(x, \(block) norm(block, "F"))
            x = lapply(seq_along(x), \(i) x[[i]]/norms[[i]])
        } else {
            norms = NULL
        }
        #print(sapply(x, apply, 2, mean))
        #print(sapply(x, norm, "F"))

        names(x) = n

        x_joint = Reduce(cbind, x)
        pca = mlr3misc::invoke(ade4::dudi.pca,
            .args=c(list(df=x_joint), private$pca_args))
        return(list(x=as.matrix(pca$li), jdr=list(pca=pca, means=means,
            norms=norms, cols=cols, standevs=standevs)))
    },

    predict_jdr = function(newx, jdr, pars) {

        n = names(newx)
        newx = lapply(names(newx),
                      function(n) newx[[n]][,jdr$cols[[n]]])
        if (pars$center){
            newx = lapply(seq_along(newx),
                function(i) scale(newx[[i]], center=jdr$means[[i]], scale=FALSE))
        }
        if (pars$scale){
            newx = lapply(seq_along(newx),
                function(i) scale(newx[[i]], center=FALSE, scale=jdr$standevs[[i]]))
        }
        if (pars$norm) {
            newx = lapply(seq_along(newx), \(i) newx[[i]]/jdr$norms[[i]])
        }
        names(newx) = n

        newx_joint = Reduce(cbind, newx)
        newpca = ade4::suprow(jdr$pca, newx_joint)
        return(list(x=as.matrix(newpca$li)))
    }

  )
)

.extralrns_dict$add("surv.pca_joint", LearnerSurvPCAJ)
