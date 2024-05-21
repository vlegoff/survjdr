#' @export
LearnerSurvIPF = R6::R6Class("LearnerSurvIPF",
  inherit = mlr3proba::LearnerSurv,
  public = list(
    #' @description
    #' Creates a new instance of this [R6][R6::R6Class] class.
    initialize = function() {
      param_set = ps(
        blocks=p_uty(tags=c("train", "predict")),
        standardize=p_lgl(default=TRUE, tags=c("train")),
        alpha=p_dbl(0, 1, default=1, tags=c("train")),
        nfolds=p_int(3L, default=5L, tags=c("train")),
        ncv=p_int(1L, default=10L, tags=c("train")),
        seed=p_int(0L, special_vals=list(NULL), default=NULL, tags=c("train"))
      )
      param_set$values$standardize = TRUE
      param_set$values$alpha = 1
      param_set$values$nfolds = 5L
      param_set$values$ncv = 10L

      super$initialize(
        id = "surv.ipflasso",
        packages = c("ipflasso", "data.table", "pec"),
        feature_types = c("integer", "numeric", "factor"),
        predict_types = c("crank", "lp", "distr"),
        param_set = param_set,
        properties = c(),
        man = "",
        label = "Regularized Generalized Canonical Correlation Analysis for Survival"
      )
    }
  ),
  private = list(
    .train = function(task) {
      # get parameters for training
      pars = self$param_set$get_values(tags="train")

      # Load data
      x = task$data(cols=task$feature_names)
      ysurv = task$truth()

      set.seed(pars$seed)

      block_list = lapply(paste0("_", pars$blocks, "$"), grep, colnames(x))

      #remove unused predictors from x
      preds = Reduce(c, block_list)
      x_pred = as.matrix(x[,..preds])
      block_list = lapply(paste0("_", pars$blocks, "$"), grep, colnames(x_pred))
      names(block_list)

      ipf_args = list(
        X=x_pred,
        Y=ysurv,
        family="cox",
        type.measure="deviance", # make it a parameter 
        standardize=pars$standardize,
        alpha=0, # ridge regression in first step
        type.step1="step",
        blocks=block_list,
        nfolds=5,
        ncv=10
      )

      #ipf_fit = do.call(ipflasso::cvr.adaptive.ipflasso, ipf_args)

      ipf_fit = ipflasso::cvr.adaptive.ipflasso(X=x_pred, Y=ysurv,
        family="cox", type.measure="deviance", standardize=pars$standardize,
        alpha=0, blocks=block_list, type.step1="sep",
        nfolds=pars$nfolds, ncv=pars$ncv)

      print(ipf_fit)

      # get non zero coefs and create a coxph model for output purpose
      coeffs = ipf_fit$coeff[-1, ipf_fit$ind.bestlambda]
      coeffs[is.na(coeffs)] <- 0 # from hermann
      nz = coeffs!=0 
      if (!any(nz)) stop("Null model fitted")
      coxph = survival::coxph(
        formula=ysurv ~ .,
        data=x[,..nz],
        init=coeffs[nz],
        iter.max=0
      )

      list(ipf=ipf_fit,
           coxph=coxph,
           y=ysurv,
           nz=nz,
           x=x[,..nz])
    },

    .predict = function(task) {
      # get parameters with tag "predict"
      pars = self$param_set$get_values(tags="predict")

      # get newdata and ensure same ordering in train and predict
      newx = mlr3extralearners:::ordered_features(task, self)
      #newx = task$data()[,.SD, .SDcols=task$feature_names] # ensure target is not in data
      fit = survival::survfit(self$model$coxph,
                              newdata=newx[,self$model$nz, with=FALSE],
                              se.fit=FALSE)
      lp = predict(self$model$coxph,
                   newdata=newx[,self$model$nz, with=FALSE])

      return(mlr3proba::.surv_return(times=fit$time,
                                     surv=t(fit$surv),
                                     lp=lp))

    }
  )
)

.extralrns_dict$add("surv.ipflasso", LearnerSurvIPF)
