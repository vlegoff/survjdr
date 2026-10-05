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
        hvg=p_dbl(0, 1, default=1, tags=c("train", "predict")),
        center=p_lgl(default=TRUE, tags=c("train", "predict")),
        scale_views=p_lgl(default=TRUE, tags=c("train", "predict")),
        likelihoods=p_uty(default=NULL, tags=c("train")),
        num_factors=p_int(1L, 100L, default=15L, tags=c("train", "predict")),
        spikeslab_factors=p_lgl(default=FALSE, tags=c("train")),
        spikeslab_weights=p_lgl(default=TRUE, tags=c("train")),
        ard_factors=p_lgl(default=FALSE, tags=c("train")),
        ard_weights=p_lgl(default=TRUE, tags=c("train")),
        iter=p_int(1L, default=1000L, tags=c("train")),
        convergence_mode=p_fct(c("fast", "medium", "slow"),
                               default="fast", tags=c("train")),
        startELBO=p_int(default=1L, tags=c("train", "predict")),
        freqELBO=p_int(default=1L, tags=c("train", "predict")),
        gpu_mode=p_lgl(default=FALSE, tags=c("train")),
        stochastic=p_lgl(default=FALSE, tags=c("train")),
        MaxIterations=p_int(1L, Inf, default=10000L, tags=c("train", "predict")),
        MinIterations=p_int(1L, Inf, default=2L, tags=c("train", "predict")), 
        ConvergenceIts=p_int(1L, Inf, default=2L, tags=c("train", "predict")),
        ConvergenceTH=p_dbl(0, 1, default = 0.0005, tags=c("train", "predict")),
        CenterTrg=p_lgl(default=FALSE, tags=c("train", "predict")),
        verbose=p_lgl(default=FALSE, tags=c("train")),
        quiet=p_lgl(default=FALSE, tags=c("train")),
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
        packages = c("MOFA2", "mlr3misc", "reticulate"),
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
    mofapy = NULL,

    # copied from MOFA2 R package
    .infer_likelihoods = function(x) {
      
      # Gaussian by default
      likelihood <- rep(x="gaussian", times=length(x))
      names(likelihood) <- names(x)
      
      for (m in names(x)) {
        data <- x[[m]]
        
        # bernoulli
        if (length(unique(data[!is.na(data)]))==2) {
          likelihood[m] <- "bernoulli"
        # poisson
        } else if (all(data[!is.na(data)]%%1==0)) {
          likelihood[m] <- "poisson"
        }
      }
      
      return(likelihood)
    },

    # Bernoulli intercept for one block
    .bernoulli_intercept_block = function(x, Z, W) {
        means = colMeans(x, na.rm=TRUE)
        intercept_naive = log(means/(1-means))
        ZW = Z %*% t(W)
        intercept = Reduce(rbind, lapply(
            1:ncol(ZW),
            \(d) private$.bernoulli_intercept_variable(x[,d], ZW[,d], intercept_naive[d])
        ))
        rownames(intercept) = colnames(x)
        return(intercept)
    },

    .bernoulli_intercept_variable = function(x, zw, start) {
        loglik = function(beta0) {
            ll = dbinom(x, size=1, prob=plogis(zw + beta0))
            return(-sum(log(ll[ll!=0])))
        }
        intercept = tryCatch({
            fit = stats4::mle(loglik, start=list(beta0=start))@coef[1]
            if (!is.finite(fit)) stop()
            return(data.frame(intercept=fit, method="MLE"))
        }, error=\(e) {
            return(data.frame(intercept=start, method="naive"))
        })
        return(intercept)
    },

    # not working for now, due to the behaviour of stats4::mle "fixed" argument
    .bernoulli_loglik = function(beta0, zw, x) {
        prob = plogis(zw + beta0)
        ll = dbinom(x, size=1, prob=prob)
        -sum(log(dens[dens!=0]))
    },

    .compute_intercepts = function(x, likelihoods, Z, W, pars) {
        intercept = list()

        for (b in names(x)) {
            if(likelihoods[[b]]=="gaussian") {
                if (pars$center) {
                    intercept[[b]] = rep(0, ncol(x[[b]]))
                    names(intercept[[b]]) = colnames(x[[b]])
                } else {
                    intercept[[b]] = colMeans(x[[b]])
                }
            } else if (likelihoods[[b]]=="bernoulli") {
                intercept[[b]] = private$.bernoulli_intercept_block(x[[b]], Z, W[[b]])
            } else {
                stop(paste0("intercept for likelihood ", likelihoods[[b]],
                    " has not been implemented"))
            }
        }

        return(intercept)

    },

    train_jdr = function(x, y, pars) {

        n = names(x)
        vars = lapply(x, function(xi) apply(xi, 2, var))
        cols_nonzero = lapply(vars, \(v) names(v[v!=0]))
        if(pars$hvg<1) {
            vars = lapply(vars, sort, decreasing=TRUE)
            cols = lapply(vars, function(vari)
                names(vari[1:ceiling(pars$hvg*length(vari))]))
            cols = mapply(function(x, y) intersect(x, y),
                x=cols, y=cols_nonzero)
            #x = lapply(x, as.matrix)
        } else {
            #cols = lapply(x, function(xi) apply(xi, 2, stats::var)!=0)
            cols = cols_nonzero
        }
        x = lapply(names(x), function(n) x[[n]][,cols[[n]]])
        names(x) = n

        likelihoods = pars$likelihoods
        if(is.null(likelihoods)) {
            likelihoods = private$.infer_likelihoods(x)
        }

        means = vector("list", length(x))
        names(means) = names(x)
        norms = vector("numeric", length(x))
        names(norms) = names(x)
        for(n in names(x)) {
            means[[n]] = colMeans(x[[n]])
            print(length(means[[n]]))
            if(likelihoods[n]=="gaussian" & pars$center) {
                x[[n]] = scale(x[[n]], center=means[[n]], scale=FALSE)
            }
            # sqrt(n*p) is needed because MOFA divides by the standard dev.
            #norms[n] = norm(x[[n]], type="F") / sqrt(ncol(x[[n]])*nrow(x[[n]]))
            # alternatively, sd works even if the data is not centered
            norms[n] = sd(x[[n]])
            if(likelihoods[n]=="gaussian" & pars$scale_views) {
                x[[n]] =  x[[n]] / norms[n]
            }
        }

        if(is.null(private$data_opts)) {
            private$data_opts = list(
              scale_views=FALSE, # managed by hand
              scale_groups=FALSE,
              center_groups=FALSE,
              use_float32=FALSE,
              views=pars$blocks
            )
        }
        if(is.null(private$model_opts)) {
            private$model_opts = list(
                likelihoods=likelihoods,
                num_factors=pars$num_factors,
                spikeslab_factors=pars$spikeslab_factors,
                spikeslab_weights=pars$spikeslab_weights,
                ard_factors=pars$ard_factors,
                ard_weights=pars$ard_weights
            )
            # Just for the likelihoods trick to work later
            private$model_opts = private$model_opts[
                !sapply(private$model_opts, is.null)
            ]
        }
        if(is.null(private$train_opts)) {
            private$train_opts = list(
                iter=pars$iter,
                convergence_mode=pars$convergence_mode,
                verbose=pars$verbose,
                quiet=pars$quiet,
                startELBO=pars$startELBO,
                freqELBO=pars$freqELBO,
                gpu_mode=pars$gpu_mode,
                stochastic=pars$stochastic
            )
        }

        if (is.null(private$mofapy)) {
            CONDA = Sys.getenv("CONDA_PREFIX")
            if (CONDA != "") {
                reticulate::use_python(file.path(CONDA, "bin/python"))
            } # else i dont know, hope for the best
            private$mofapy = reticulate::import("mofapy2.run.entry_point")
            message("reticulate loaded")
        }

        ent = private$mofapy$entry_point()

        ent$set_data_options(
            scale_views=private$data_opts$scale_views,
            scale_groups=private$data_opts$scale_groups,
            center_groups=private$data_opts$center_groups,
            use_float32=private$data_opts$use_float32
        )


        message("setting data matrix")
        samples = rownames(x[[1]])
        if (is.null(samples)) {
            samples = paste0("sample_", 1:nrow(x[[1]]))
            for (n in names(x)) rownames(x[[n]]) = samples
        }
        ent$set_data_matrix(
            data=lapply(x, \(bl) list(bl)),
            likelihoods=private$model_opts$likelihoods,
            views_names=names(x),
            samples_names=list(samples),
            features_names=unname(lapply(x, colnames))
        )
        message("data matrix set")

        ent$set_model_options(
            factors=private$model_opts$num_factors,
            spikeslab_factors=private$model_opts$spikeslab_factors,
            spikeslab_weights=private$model_opts$spikeslab_weights,
            ard_factors=private$model_opts$ard_factors,
            ard_weights=private$model_opts$ard_weights
        )

        ent$set_train_options(
            iter=private$train_opts$iter,
            startELBO=private$train_opts$startELBO,
            freqELBO=private$train_opts$freqELBO,
            convergence_mode=private$train_opts$convergence_mode,
            verbose=private$train_opts$verbose,
            quiet=private$train_opts$quiet,
            gpu_mode=private$train_opts$gpu_mode,
            dropR2=-1 # just to be sure not factor is dropped
        )

        ent$build()
        ent$run()

        # Manage temp file used for training...
        #outfile = tempfile(patter="mofa_model_", fileext=".hdf5")
        #ent$save(outfile=outfile, save_data=FALSE, expectations="all")

        #fctrzn = MOFA2::load_model(file=outfile)

        Z = ent$model$nodes$Z$getExpectation()
        print(Z)
        colnames(Z) = paste0("Factor", 1:ncol(Z))
        rownames(Z) = samples

        W = ent$model$nodes$W$getExpectation()
        names(W) = pars$blocks
        for (n in names(W)) {
            rownames(W[[n]]) = colnames(x[[n]])
            colnames(W[[n]]) = colnames(Z)
        }

        Tau = ent$model$nodes$Tau$getExpectation()
        names(Tau) = pars$blocks
        for (n in names(Tau)) {
            colnames(Tau[[n]]) = colnames(x[[n]])
            rownames(Tau[[n]]) = samples
        }


        intercepts = private$.compute_intercepts(x, likelihoods, Z, W, pars)

        return(list(x=Z, jdr=list(cols=cols, means=means,
            norms=norms, intercepts=intercepts, Z=Z, W=W, Tau=Tau,
            likelihoods=likelihoods)))
    },

    predict_jdr = function(newx, jdr, pars) {

        print("predict")

        #mofa = load_model(file=jdr$mofa)
        
        if(!is.null(jdr$cols)) {
            n = names(newx)
            newx = lapply(names(newx),
                          function(n) newx[[n]][,jdr$cols[[n]]])
            names(newx) = n
        }

        likelihoods = jdr$likelihoods
        for(n in names(newx)) {
            if(likelihoods[n]=="gaussian" & pars$center) {
                newx[[n]] = scale(newx[[n]], center=jdr$means[[n]], scale=FALSE)
            }
            if(likelihoods[n]=="gaussian" & pars$scale_view) {
                newx[[n]] =  newx[[n]] / jdr$norms[n]
            }
        }

        # Create list of parameters for MOTL
        TL_param = list()
        TL_param$YTrg = newx
        TL_param$Fctrzn_Lrn_W0 = lapply(jdr$intercepts, \(x) {
            if (class(x) == "data.frame"){
                cur = x$intercept
                names(cur) = rownames(x)
                return(cur)
            } else return(x) 
            })
        TL_param$Tau = jdr$Tau

        for (view in names(newx)){
            TL_param$Fctrzn_Lrn_W[[view]] = jdr$W[[view]]
            TL_param$Fctrzn_Lrn_WSq[[view]] = TL_param$Fctrzn_Lrn_W[[view]]^2
            print(any(!is.finite(TL_param$Tau[[view]])))
            TL_param$Tau[[view]] = colMeans(TL_param$Tau[[view]], na.rm = T)
            TL_param$Tau[[view]] = matrix(TL_param$Tau[[view]],
                nrow = dim(newx[[view]])[1], ncol = dim(newx[[view]])[2],
                byrow = T)
            colnames(TL_param$Tau[[view]]) = colnames(newx[[view]])
            TL_param$TauLn[[view]] = numeric()
            if(likelihoods[[view]] == "gaussian"){
                TL_param$TauLn[[view]] = log(TL_param$Tau[[view]])
            }
        }
        print("TL param ready")

        TL_data=transferLearning_function(TL_param = TL_param, 
            MaxIterations=pars$MaxIterations, 
            MinIterations=pars$MinIterations, 
            minFactors=pars$num_factors,       # To be sure it is always verified
            StartDropFactor=1e10,    # To be almost sure it never happens
            FreqDropFactor=1e10,    # If it happens it happens really not often
            StartELBO=pars$startELBO, 
            FreqELBO=pars$freqELBO, 
            DropFactorTH=0,       # No factor can have an explained variance below zero so no factor can be dropped
            ConvergenceIts=pars$ConvergenceIts, 
            ConvergenceTH=pars$ConvergenceTH, 
            CenterTrg=pars$CenterTrg,
            likelihoods=likelihoods,
            Z=jdr$Z)

        return(list(x=TL_data$ZMu))

        #newx = lapply(newx, t)

        #newx = Reduce(cbind, newx)
        #W = Reduce(cbind, lapply(jdr$mofa@expectations$W, t))
        #Z_new = newx %*% MASS::ginv(W) # fastest way

        #return(list(x=Z_new))
    }

  )
)

.extralrns_dict$add("surv.mofa", LearnerSurvMOFA)
