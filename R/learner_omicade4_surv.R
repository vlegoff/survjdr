#' @export
LearnerSurvMCIA = R6::R6Class("LearnerSurvMCIA",
  inherit = LearnerSeqMod,
  public = list(
    #' @description
    #' Creates a new instance of this [R6][R6::R6Class] class.
    initialize = function() {
      param_set = ps(
        blocks=p_uty(tags=c("train", "predict")),
        clinical_fav=p_lgl(default=TRUE, tags=c("train", "predict")),
        cia.nf=p_int(2L, default=2L, tags=c("train")),
        pred_mode=p_fct(c("Ntrain", "Ntest"), default="Ntest",
                        tags=c("train", "predict")),
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
        id = "surv.omicade4",
        packages = c("omicade4", "mlr3misc"),
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

    mcia_args = NULL,

    train_jdr = function(x, y, pars) {
        
        if(is.null(private$mcia_args)) {
            private$mcia_args = list(
              cia.nf=pars$cia.nf,
              cia.scan=FALSE,
              nsc=TRUE,
              svd=TRUE
            )
        }

        x = lapply(x, t) # why oh god why
        # columns to remove
        ind = lapply(x, apply, 1, function(d) all(d==min(d)))
        to_remove = lapply(ind, function(d) which(d))
        x = mapply(function(X, Y) if (length(Y)>0) X[-Y,] else X,
            X=x, Y=to_remove)
        mins = lapply(x, min) # it's a surprise tool that will help us later

        mcia_fit = mlr3misc::invoke(omicade4::mcia,
          .args=c(list(df.list=x), private$mcia_args))

        SynVar = as.matrix(mcia_fit$mcoa$SynVar)
        Tl1 = Reduce(cbind, lapply(pars$blocks, \(x) {
            ret = mcia_fit$mcoa$Tl1[grepl(paste0("\\.", x, "$"),
                rownames(mcia_fit$mcoa$Tl1)),]
            colnames(ret) = sub("Axis", x, colnames(ret))
            ret
          }))
        latent_space = as.matrix(cbind(SynVar, Tl1))

        return(list(x=latent_space,
                    jdr=list(fit=mcia_fit, to_remove=to_remove, mins=mins)))
    },

    predict_jdr = function(newx, jdr, pars) {
        newx = lapply(newx, t)
        newx = mapply(function(X, Y) if (length(Y)>0) X[-Y,] else X,
            X=newx, Y=jdr$to_remove)
        newx = lapply(newx, as.data.frame)

        mode = ifelse(pars$pred_mode=="Ntrain", "1", "2")
        predict = predict_omicade4(newx, jdr$fit, jdr$mins, mode=mode)
        SynVar = as.matrix(predict$SynVar)
        Tl1 = Reduce(cbind, lapply(pars$blocks, \(x) {
            ret = predict$Tl1[grepl(paste0("\\.", x, "$"),
                rownames(predict$Tl1)),]
            colnames(ret) = sub("Axis", x, colnames(ret))
            rownames(ret) = NULL
            ret
          }))
        latent_space = cbind(SynVar, Tl1)

        return(list(x=latent_space))
    }

  )
)

.extralrns_dict$add("surv.omicade4", LearnerSurvMCIA)
