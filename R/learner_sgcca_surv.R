#' @export
LearnerSurvSGCCA = R6::R6Class("LearnerSurvSGCCA",
  inherit = LearnerSeqMod,
  public = list(
    #' @description
    #' Creates a new instance of this [R6][R6::R6Class] class.
    initialize = function() {
      param_set = ps(
        blocks=p_uty(tags=c("train", "predict")),
        sparsity=p_dbl(0, 1, default=.5, tags=c("train")),
        clinical_fav=p_lgl(default=TRUE, tags=c("train", "predict")),
        supervised=p_lgl(default=FALSE, tags=c("train", "predict")),
        ncomp=p_int(1L, default=1L, tags=c("train")),
        scheme=p_fct(c("horst", "factorial", "centroid"), default="factorial",
                      tags=c("train")),
        nfolds=p_int(1L, default=10L, tags=c("train")),
        nlambdas=p_int(10L, 1000L, default=100L, tags=c("train")),
        CV_measure=p_fct(c("cindex", "ibs", "ibsRR",  "auc", "C", "deviance",
                           "basic", "V&VH", "linpred"), default="cindex",
                          tags=c("train")),
        seed=p_int(0L, special_vals=list(NULL), default=NULL, tags=c("train"))
      )
      param_set$values = param_set$default

      super$initialize(
        id = "surv.rgcca",
        packages = c("RGCCA", "mlr3misc"),
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

    rgcca_args = NULL,

    train_jdr = function(x, y, pars) {
        
        if(is.null(private$rgcca_args)) {
            if(!pars$supervised) {
                complete_matrix = matrix(1, length(pars$blocks),
                                       length(pars$blocks))
                diag(complete_matrix) = 0
                rownames(complete_matrix) = pars$blocks
                colnames(complete_matrix) = pars$blocks
            }

            # sparsity percentage per block
            print(sapply(x, ncol))
            sparsity = pars$sparsity/sapply(x, ncol) + (1-pars$sparsity)
            if (pars$supervised) sparsity = c(sparsity, 1)
            print(sparsity)

            private$rgcca_args = list(
                response=if(pars$supervised) length(x) + 1,
                connection=if(!pars$supervised) complete_matrix,
                sparsity=sparsity,
                ncomp=pars$ncomp,
                scheme=pars$scheme,
                method="rgcca",
                scale=TRUE,
                scale_block="inertia",
                verbose=F
            )
        }

        if(pars$supervised) {
            null_mod = survival::coxph(y~1)
            x[["residuals"]] = residuals(null_mod, type="deviance")
        }

        rgcca_fit = mlr3misc::invoke(RGCCA::rgcca,
            .args=c(list(blocks=x), private$rgcca_args))
        comps = Reduce(cbind,
            rgcca_fit$Y[names(rgcca_fit$Y)!="residuals"])
        colnames(comps) = paste(
            rep(pars$blocks, each=private$rgcca_args$ncomp),
            rep(1:private$rgcca_args$ncomp, private$rgcca_args$ncomp),
            sep="."
        )
        #for (b in names(rgcca_fit$blocks)) { # trick to avoid storing blocks, not working though
            #rgcca_fit$blocks[[b]] = 1
            #rgcca_fit$call$blocks[[b]] = 1
        #}
        return(list(x=comps, jdr=rgcca_fit))
    },

    predict_jdr = function(newx, jdr, pars) {
        if(pars$supervised) {
            newx[["residuals"]] = as.matrix(rep(1, nrow(newx[[1]])))
            colnames(newx[["residuals"]]) = "residuals"
        }
        rgcca_pred = RGCCA::rgcca_transform(jdr, newx)
        pred_space = Reduce(cbind, rgcca_pred[names(rgcca_pred)!="residuals"])
        colnames(pred_space) = paste(
            rep(pars$blocks, each=private$rgcca_args$ncomp),
            rep(1:private$rgcca_args$ncomp, private$rgcca_args$ncomp),
            sep="."
        )
        return(list(x=pred_space))
    }

  )
)

.extralrns_dict$add("surv.sgcca", LearnerSurvSGCCA)
