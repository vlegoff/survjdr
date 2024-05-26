#' @export
LearnerSurvMOFA = R6::R6Class("LearnerSurvMOFA",
  inherit = LearnerSeqMod,
  public = list(
    #' @description
    #' Creates a new instance of this [R6][R6::R6Class] class.
    initialize = function() {
      param_set = ps(
        blocks=p_uty(tags=c("train", "predict")),
        clinical_fav=p_lgl(default=TRUE, tags=c("train", "predict")),
        scale_views=p_lgl(default=TRUE, tags=c("train")),
        likelihoods=p_uty(default=NULL, tags=c("train")),
        num_factors=p_int(1L, 100L, default=15L, tags=c("train")),
        spikeslab_factors=p_lgl(default=FALSE, tags=c("train")),
        spikeslab_weights=p_lgl(default=TRUE, tags=c("train")),
        ard_factors=p_lgl(default=FALSE, tags=c("train")),
        ard_weights=p_lgl(default=TRUE, tags=c("train")),
        maxiter=p_int(1L, default=1000L, tags=c("train")),
        convergence_mode=p_fct(c("fast", "medium", "slow"),
                               default="fast", tags=c("train")),
        startELBO=p_int(default=1L, tags=c("train")),
        freqELBO=p_int(default=1L, tags=c("train")),
        gpu_mode=p_lgl(default=FALSE, tags=c("train")),
        stochastic=p_lgl(default=FALSE, tags=c("train")),
        verbose=p_lgl(default=FALSE, tags=c("train")),
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
        id = "surv.mofa",
        packages = c("MOFA2", "mlr3misc"),
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

    data_opts = NULL,
    model_opts = NULL,
    train_opts = NULL,

    train_jdr = function(x, y, pars) {

        x = lapply(x, t) # why oh god why

        MOFAobject = MOFA2::create_mofa(x)
        print(MOFAobject)

        if(is.null(private$data_opts)) {
            private$data_opts = list(
              scale_views=pars$scale_views,
              scale_groups=FALSE,
              center_groups=TRUE,
              use_float32=FALSE,
              views=pars$blocks
            )
            # Making it more robust
            default_options = MOFA2::get_default_data_options(MOFAobject)
            private$data_opts = c(
                private$data_opts,
                default_options[
                    !names(default_options) %in% names(private$data_opts)
                ]
            )
            private$data_opts = private$data_opts[names(default_options)]
        }
        if(is.null(private$model_opts)) {
            private$model_opts = list(
                likelihoods=pars$likelihoods,
                num_factors=pars$num_factors,
                spikeslab_factors=pars$spikeslab_factors,
                spikeslab_weights=pars$spikeslab_weights,
                ard_factors=pars$ard_factors,
                ard_weights=pars$ard_weights
            )
            default_options = MOFA2::get_default_model_options(MOFAobject)
            private$model_opts = c(
                private$model_opts,
                default_options[
                    !names(default_options) %in% names(private$model_opts)
                ]
            )
            private$model_opts = private$model_opts[names(default_options)]
        }
        if(is.null(private$train_opts)) {
            private$train_opts = list(
                maxiter=pars$maxiter,
                convergence_mode=pars$convergence_mode,
                verbose=pars$verbose,
                startELBO=pars$startELBO,
                freqELBO=pars$freqELBO,
                gpu_mode=pars$gpu_mode,
                stochastic=pars$stochastic
            )
            default_options = MOFA2::get_default_training_options(MOFAobject)
            private$train_opts = c(
                private$train_opts,
                default_options[
                    !names(default_options) %in% names(private$train_opts)
                ]
            )
            private$train_opts = private$train_opts[names(default_options)]
        }

        #print(MOFA2::get_default_data_options(MOFAobject))
        #print("***********")
        #print(private$data_opts)

        MOFAobject = MOFA2::prepare_mofa(
            object=MOFAobject,
            data_options=private$data_opts,
            model_options=private$model_opts,
            training_options=private$train_opts
        )

        # Manage temp file used for training...

        outfile = file.path(tempdir(), "mofa_model.hdf5")
        MOFAobject_trained = MOFA2::run_mofa(MOFAobject, outfile,
            use_basilisk=TRUE)
        #MOFAobject_trained = MOFA2::run_mofa(MOFAobject, use_basilisk=TRUE)

        latent_space = Reduce(rbind, MOFAobject_trained@expectations$Z)
        #print(latent_space)

        return(list(x=latent_space, jdr=MOFAobject_trained))
    },

    predict_jdr = function(newx, jdr, pars) {
        #newx = lapply(newx, t)
        newx = Reduce(cbind, newx)

        #print(lapply(jdr@expectations$W, dim))
        W = Reduce(cbind, lapply(jdr@expectations$W, t))
        #W = Reduce(rbind, jdr@expectations$W)
        #print(dim(newx))
        #print(dim(W))
        #lm1 = lm(t(newx) ~ t(W) - 1)
        #print(coef(lm1))
        #Z_new = t(as.matrix(coef(lm1)))
        #Z2_new = newx %*% t(W) %*% MASS::ginv(W %*% t(W)) 
        Z2_new = newx %*% MASS::ginv(W) # fastest way
        #print(dim(Z2_new))
        #print(Z_new - Z2_new)

        return(list(x=Z2_new))
    }

  )
)

.extralrns_dict$add("surv.mofa", LearnerSurvMOFA)
