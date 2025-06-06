#' @export
LearnerSurvNMF = R6::R6Class("LearnerSurvNMF",
  inherit = LearnerSeqMod,
  public = list(
    #' @description
    #' Creates a new instance of this [R6][R6::R6Class] class.
    initialize = function() {
      param_set = ps(
        blocks=p_uty(tags=c("train", "predict")),
        clinical_fav=p_lgl(default=TRUE, tags=c("train", "predict")),
        k=p_int(2L, 100L, default=2L, tags=c("train", "predict")),
        weights=p_uty(default=NULL, tags=c("train")),
        nfolds=p_int(1L, default=10L, tags=c("train")),
        nlambdas=p_int(10L, 1000L, default=100L, tags=c("train")),
        CV_measure=p_fct(c("cindex", "ibs", "ibsRR",  "auc", "C", "deviance",
                           "basic", "V&VH", "linpred"), default="cindex",
                          tags=c("train")),
        cv_save_path=p_uty(default=NULL, tags=c("predict")),
        n.ini=p_int(1L, default=30L, tags=c("train")),
        maxiter=p_int(10L, default=200L, tags=c("train")),
        seed=p_int(0L, special_vals=list(NULL), default=NULL, tags=c("train"))
      )
      param_set$values = param_set$default

      super$initialize(
        id = "surv.intnmf",
        packages = c("IntNMF", "mlr3misc"),
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

    nmf_args = NULL,
    train_jdr = function(x, y, pars) {

        if(is.null(private$nmf_args)) {
            private$nmf_args = list(
              k=pars$k,
              maxiter=pars$maxiter,
              st.count=20,
              n.ini=pars$n.ini,
              ini.nndsvd=FALSE,
              seed=TRUE
            )
        }

        #Min-Max normalization, cantini style
        n = names(x)

        mins = sapply(x, min)
        mins = abs(mins)*(mins<0)
        x = lapply(seq_along(x), function(i) x[[i]] + mins[i])

        maxs = sapply(x, max)
        x = lapply(seq_along(x), function(i) x[[i]] / maxs[i])

        names(x) = n

        #Remove bad columns
        cols = lapply(x, function(xi) apply(xi, 2, stats::var)!=0)
        x = lapply(names(x), function(n) x[[n]][,cols[[n]]])

        # computing weights
        if(is.null(pars$weights)) weights=rep(1, length(pars$blocks))
        else if(pars$weights=="frob") { 
            norms=sapply(x, norm, type="F")^2
            weights=max(norms)/norms
        } else if(pars$weights=="frob_p") {
            p=sapply(x, ncol)
            norms=sapply(x, norm, type="F")^2
            norms=norms/p
            weights=max(norms)/norms
        }
        else if(length(blocks)==length(pars$weights)) weights=pars$weights
        else stop("weights and blocks should have the same length")

        nmf_fit = mlr3misc::invoke(IntNMF::nmf.mnnals,
          .args=c(list(dat=x, wt=weights), private$nmf_args))

        return(list(x=nmf_fit$W, jdr=list(fit=nmf_fit$H, cols=cols, wt=weights,
            mins=mins, maxs=maxs)))
    },

    predict_jdr = function(newx, jdr, pars) {
        newx = lapply(names(newx),
                      function(n) newx[[n]][,jdr$cols[[n]]])
        n = names(newx)
        newx = lapply(seq_along(newx), function(i) newx[[i]] + jdr$mins[i])
        newx = lapply(seq_along(newx), function(i) newx[[i]] / jdr$maxs[i])
        names(newx) = n

        latent_space = IntNMF:::W.fcnnls(x=jdr$fit, y=newx, weight=jdr$wt)
        latent_space = t(latent_space$coef)
        return(list(x=latent_space))
    }

  )
)

.extralrns_dict$add("surv.intnmf", LearnerSurvNMF)
