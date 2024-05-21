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
        nfolds=p_int(1L, default=10L, tags=c("train")),
        nlambdas=p_int(10L, 1000L, default=100L, tags=c("train")),
        CV_measure=p_fct(c("cindex", "ibs", "ibsRR",  "auc", "C", "deviance",
                           "basic", "V&VH", "linpred"), default="cindex",
                          tags=c("train")),
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
    mins = NULL,
    cols = NULL,

    train_jdr = function(x, y, pars) {

        if(is.null(private$nmf_args)) {
            private$nmf_args = list(
              k=pars$k,
              maxiter=200,
              st.count=20,
              n.ini=30,
              ini.nndsvd=FALSE,
              seed=TRUE
            )
        }

        if(is.null(private$mins)) {
            private$mins = sapply(x, min)
            private$mins = abs(private$mins)*(private$mins<0)
        }

        #Make the matrices non negative
        n = names(x)
        x = lapply(seq_along(x), function(i) x[[i]] + private$mins[i])
        names(x) = n
        #Remove bad columns
        #private$cols = lapply(x, function(xi) apply(xi, 2, var)!=0)
        cols = lapply(x, function(xi) apply(xi, 2, var)!=0)
        x = lapply(names(x), function(n) x[[n]][,cols[[n]]])

        nmf_fit = mlr3misc::invoke(IntNMF::nmf.mnnals,
          .args=c(list(dat=x), private$nmf_args))

        #print(nmf_fit$W)
        print("trained")

        return(list(x=nmf_fit$W, jdr=list(fit=nmf_fit$H, cols=cols)))
    },

    predict_jdr = function(newx, jdr, pars) {
        #if(!is.null(private$cols)) {
          #newx = lapply(names(newx),
                        #function(n) newx[[n]][,private$cols[[n]]])
        #}
        newx = lapply(names(newx),
                      function(n) newx[[n]][,jdr$cols[[n]]])
        n = names(newx)
        newx = lapply(seq_along(newx), function(i) newx[[i]] + private$mins[i])
        names(newx) = n

        XHt = 0
        HHt = 0
        for(i in seq_along(newx)) {
            XHt = XHt + newx[[i]] %*% t(jdr$fit[[i]]) 
            HHt = HHt + jdr$fit[[i]] %*% t(jdr$fit[[i]])
        }
        latent_space2 = XHt %*% MASS::ginv(HHt)
        #latent_space2 = latent_space2 + abs(min(latent_space2))

        latent_space4 = IntNMF:::W.fcnnls(x=jdr$fit, y=newx,
                                          weight=rep(1, length(jdr$fit)))
        latent_space4 = t(latent_space4$coef)

        newx = Reduce(cbind, newx)
        H = Reduce(cbind, jdr$fit)
        print(dim(H))
        print(dim(newx))
        # On the right track but need to use non negative least squares
        #lm1 = lm(t(newx) ~ t(H) - 1)
        # casting matrix nnls as a serie of nnls:
        #nnlm = vector(mode="list", length=private$nmf_args$k)
        latent_space = matrix(0, ncol=private$nmf_args$k, nrow=nrow(newx))
        for(i in 1:nrow(latent_space)) {
            nn = nnls::nnls(t(H), newx[i,])
            #latent_space[,i] = coef(nnls::nnls(t(H), t(newx)[,i]))
            latent_space[i,] = coef(nn)
        }
        latent_space3 = newx %*% t(H) %*% MASS::ginv(H %*% t(H))
        #print(abs(min(latent_space3)))
        #latent_space3 = latent_space3 + abs(min(latent_space3))
        print(latent_space2)
        print(latent_space3)
        print(latent_space)
        print(latent_space4)
        print(latent_space3 - latent_space)
        print(latent_space2 - latent_space4)

        # latent_space4 seems closer to latent_space_train than the others
        return(list(x=latent_space4))
    }

  )
)

.extralrns_dict$add("surv.intnmf", LearnerSurvNMF)
