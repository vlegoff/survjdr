
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
LearnerSurvSJIVE = R6::R6Class("LearnerSurvSJIVE",
  inherit = LearnerSeqMod,
  public = list(
    #' @description
    #' Creates a new instance of this [R6][R6::R6Class] class.
    initialize = function() {
      param_set = ps(
        blocks=p_uty(tags=c("train", "predict")),
        clinical_fav=p_lgl(default=TRUE, tags=c("train", "predict")),
        rankJ=p_int(0L, default=1L, tags=c("train", "predict")),
        rankA=p_int(0L, default=1L, tags=c("train", "predict")),
        eta=p_dbl(0, 1, default=0, tags=c("train")),
        nfolds=p_int(1L, default=10L, tags=c("train")),
        nlambdas=p_int(10L, 1000L, default=100L, tags=c("train")),
        CV_measure=p_fct(c("cindex", "ibs", "ibsRR",  "auc", "C", "deviance",
                           "basic", "V&VH", "linpred"), default="cindex",
                          tags=c("train")),
        seed=p_int(0L, special_vals=list(NULL), default=NULL, tags=c("train")),
        max.iter=p_int(1L, default=1000L, tags=c("train"))
      )
      param_set$values = param_set$default

      super$initialize(
        id = "surv.sjive",
        packages = c("sup.r.jive", "mlr3misc"),
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
    cols = NULL,

    train_jdr = function(x, y, pars) {
        
        if(is.null(private$jive_args)) {
            private$jive_args = list(
              rankJ=pars$rankJ,
              rankA=rep(pars$rankA, length(x)),
              eta=pars$eta,
              center.scale=TRUE,
              reduce.dim=TRUE,
              max.iter=pars$max.iter
            )
        }

        x = lapply(x, t) # why oh god why
        # removing all 0 variables to avoid sjive bug
        #private$cols = lapply(x, function(xi) as.logical(rowSums(xi!=0)))
        #private$cols = lapply(x, function(xi) apply(xi, 1, var)!=0)
        cols = lapply(x, function(xi) apply(xi, 1, var)!=0)
        x = lapply(names(x), function(n) x[[n]][cols[[n]],])

        null_mod = survival::coxph(y~1)
        resi = residuals(null_mod, type="deviance")

        jive_fit = mlr3misc::invoke(sup.r.jive::sJIVE,
          .args=c(list(X=x, Y=resi), private$jive_args))

        latent_space = matrix(0,
          nrow=ncol(x[[1]]),
          ncol=private$jive_args$rankJ+length(x)*private$jive_args$rankA[1])

        if(private$jive_args$rankA[1]==0) {
          latent_space[,1:private$jive_args$rankJ] = t(jive_fit$S_J)
          colnames(latent_space) = paste0("common.", 1:pars$rankJ)
        } else if (private$jive_args$rankJ==0) {
          latent_space[,1:ncol(latent_space)] = 
            t(Reduce(rbind, jive_fit$S_I)) # cbind was not working

          colnames(latent_space) =
              paste0(
                  rep(pars$blocks, each=pars$rankA),
                  ".",
                  rep(1:pars$rankA, pars$rankA)
              )
        } else {
          latent_space[,1:private$jive_args$rankJ] = t(jive_fit$S_J)
          
          latent_space[,(private$jive_args$rankJ+1):ncol(latent_space)] = 
            t(Reduce(rbind, jive_fit$S_I)) # cbind was not working
          colnames(latent_space) = c(
              paste0("common.", 1:pars$rankJ),
              paste0(
                  rep(pars$blocks, each=pars$rankA),
                  ".",
                  rep(1:pars$rankA, pars$rankA)
              )
          )
        }

        return(list(x=latent_space, jdr=list(fit=jive_fit, cols=cols)))
    },

    predict_jdr = function(newx, jdr, pars) {
        newx = lapply(newx, t)
        #if(!is.null(private$cols)) {
        newx = lapply(names(newx),
                      function(n) newx[[n]][jdr$cols[[n]],])
        #}

        predict = predict(jdr$fit, newdata=newx)

        latent_space = matrix(0,
          nrow=ncol(newx[[1]]),
          ncol=private$jive_args$rankJ+length(newx)*private$jive_args$rankA[1])

        #latent_space[,1:private$jive_args$rankJ] = t(predict$Sj)
        
        #latent_space[,(private$jive_args$rankJ+1):ncol(latent_space)] = 
          #t(Reduce(rbind, predict$Si)) # cbind was not working

        if(private$jive_args$rankA[1]==0) {
          latent_space[,1:private$jive_args$rankJ] = t(predict$Sj)
          colnames(latent_space) = paste0("common.", 1:pars$rankJ)
        } else if (private$jive_args$rankJ==0) {
          latent_space[,1:ncol(latent_space)] = 
            t(Reduce(rbind, predict$Si)) # cbind was not working

          colnames(latent_space) =
              paste0(
                  rep(pars$blocks, each=pars$rankA),
                  ".",
                  rep(1:pars$rankA, pars$rankA)
              )
        } else {
          print(predict$Sj)
          latent_space[,1:private$jive_args$rankJ] = t(predict$Sj)
          
          latent_space[,(private$jive_args$rankJ+1):ncol(latent_space)] = 
            t(Reduce(rbind, predict$Si)) # cbind was not working
          colnames(latent_space) = c(
              paste0("common.", 1:pars$rankJ),
              paste0(
                  rep(pars$blocks, each=pars$rankA),
                  ".",
                  rep(1:pars$rankA, pars$rankA)
              )
          )
        }

        return(list(x=latent_space))
    }

  )
)

mlr3::mlr_learners$add("surv.sjive", LearnerSurvSJIVE)
