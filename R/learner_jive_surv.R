
#' @title Survival Rgcca For Dimension Reduction Followed By Glmnet For Surv Prediction Learner
#' @author Unknown
#' @name mlr_learners_surv.rgcca
#'
#' @description
#' FIXME: BRIEF DESCRIPTION OF THE LEARNER.
#' Calls [RGCCA::RGCCA()] from FIXME: (CRAN VS NO CRAN): \CRANpkg{RGCCA} | 'RGCCA'.
#'
#' @section Initial parameter values:
#' FIXME: DEVIATIONS FROM UPSTREAM PARAMETERS. DELETE IF NOT APPLICABLE.
#'
#' @section Custom mlr3 parameters:
#' FIXME: DEVIATIONS FROM UPSTREAM DEFAULTS. DELETE IF NOT APPLICABLE.
#'
#' @templateVar id surv.rgcca
#' @template learner
#'
#' @references
#' `r format_bib(FIXME: ONE OR MORE REFERENCES FROM bibentries.R)`
#'
#' @template seealso_learner
#' @template example
#' @export
LearnerSurvJIVE = R6::R6Class("LearnerSurvJIVE",
  inherit = LearnerSeqMod,
  public = list(
    #' @description
    #' Creates a new instance of this [R6][R6::R6Class] class.
    initialize = function() {
      param_set = ps(
        blocks=p_uty(tags=c("train", "predict")),
        clinical_fav=p_lgl(default=TRUE, tags=c("train", "predict")),
        rankJ=p_int(1L, default=1L, tags=c("train", "predict")),
        rankA=p_int(1L, default=1L, tags=c("train", "predict")),
        nfolds=p_int(1L, default=10L, tags=c("train")),
        nlambdas=p_int(10L, 1000L, default=100L, tags=c("train")),
        CV_measure=p_fct(c("cindex", "ibs", "ibsRR",  "auc", "C", "deviance",
                           "basic", "V&VH", "linpred"), default="cindex",
                          tags=c("train")),
        seed=p_int(0L, special_vals=list(NULL), default=NULL, tags=c("train"))
      )
      param_set$values = param_set$default

      super$initialize(
        id = "surv.jive",
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

    jive_args = NULL,

    train_jdr = function(x, y, pars) {
        
        if(is.null(private$jive_args)) {
            private$jive_args = list(
              rankJ=pars$rankJ,
              rankA=rep(pars$rankA, length(x)),
              method="given",
              orthIndiv=FALSE,
              scale=TRUE,
              center=TRUE
            )
        }

        x = lapply(x, t) # why oh god why

        jive_fit = mlr3misc::invoke(r.jive::jive,
          .args=c(list(data=x), private$jive_args))
        predict = r.jive::jive.predict(data.new=x, jive_fit)

        latent_space = matrix(0,
          nrow=ncol(x[[1]]),
          ncol=private$jive_args$rankJ+length(x)*private$jive_args$rankA)

        latent_space[,1:private$jive_args$rankJ] = t(predict$joint.scores)
        
        latent_space[,(private$jive_args$rankJ+1):ncol(latent_space)] = 
          t(Reduce(rbind, predict$indiv.scores)) # cbind was not working

        colnames(latent_space) = c(
            paste0("common.", 1:pars$rankJ),
            paste0(
                rep(pars$blocks, each=pars$rankA),
                ".",
                rep(1:pars$rankA, pars$rankA)
            )
        )

        return(list(x=latent_space, jdr=jive_fit))
    },

    predict_jdr = function(newx, jdr, pars) {
        newx = lapply(newx, t)

        predict = r.jive::jive.predict(data.new=newx, jdr)

        latent_space = matrix(0,
          nrow=ncol(newx[[1]]),
          ncol=private$jive_args$rankJ+length(newx)*private$jive_args$rankA)

        latent_space[,1:private$jive_args$rankJ] = t(predict$joint.scores)
        
        latent_space[,(private$jive_args$rankJ+1):ncol(latent_space)] = 
          t(Reduce(rbind, predict$indiv.scores)) # cbind was not working

        colnames(latent_space) = c(
            paste0("common.", 1:pars$rankJ),
            paste0(
                rep(pars$blocks, each=pars$rankA),
                ".",
                rep(1:pars$rankA, pars$rankA)
            )
        )

        return(list(x=latent_space))
    }

  )
)

mlr3::mlr_learners$add("surv.jive", LearnerSurvJIVE)
