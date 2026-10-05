#' @export
LearneriCluster2 = R6::R6Class("LearneriCluster2",
  inherit = LearnerSeqMod,
  public = list(
    #' @description
    #' Creates a new instance of this [R6][R6::R6Class] class.
    initialize = function() {
      param_set = ps(
        blocks=p_uty(tags=c("train", "predict")),
        clinical_fav=p_lgl(default=TRUE, tags=c("train", "predict")),
        K=p_int(2L, 100L, default=2L, tags=c("train", "predict")),
        s=p_dbl(lower=0, upper=1, special_vals=list(NULL), default=NULL,
            tags=c("train")),
        lambda=p_dbl(lower=0, default=0, tags=c("train")), # for iCluster
        nlambdas=p_int(10L, 1000L, default=100L, tags=c("train")), # for glmnet
        nfolds=p_int(1L, default=10L, tags=c("train")),
        CV_measure=p_fct(c("cindex", "ibs", "ibsRR",  "auc", "C", "deviance",
                           "basic", "V&VH", "linpred"), default="cindex",
                          tags=c("train")),
        cv_save_path=p_uty(default=NULL, tags=c("train", "predict")),
        maxiter=p_int(10L, default=200L, tags=c("train")),
        eps=p_dbl(lower=0, upper=Inf, default=1e-04, tags=c("train")),
        eps2=p_dbl(lower=0, upper=Inf, default=1e-08, tags=c("train")),
        center=p_lgl(default=TRUE, tags=c("train")),
        seed=p_int(0L, special_vals=list(NULL), default=NULL, tags=c("train"))
      )
      param_set$values = param_set$default

      super$initialize(
        id = "surv.icluster2",
        packages = c("iClusterPlus", "mlr3misc"),
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

    iclust2_args = NULL,
    train_jdr = function(x, y, pars) {

        message("train")
        if(is.null(private$iclust2_args)) {
            cancers = c("LAML", "ESCA", "PAAD", "SARC", "LIHC", "COAD", "KIRP",
                "OV", "KIRC", "SKCM", "STAD", "BLCA", "UCEC", "LGG", "HNSC",
                "LUSC", "LUAD", "BRCA")
            if(!is.null(pars$s) & (pars$task_id %in% cancers) & (pars$K %in% c(2, 6, 11))) {
                pars$s = 1-pars$s # to have as input: 0 least sparse, 1 most sparse
                if (pars$s==1) { # no sparsity at all
                    message("no sparsity asked")
                    lambda = 0
                } else {
                    # linearly interpolating lambda from known values
                    message("linear interpolation")
                    lambda_refs = read.csv(
                        system.file("icluster_sparsity_results.csv", package="survjdr")
                    )
                    lambda_refs = subset(lambda_refs, K==pars$K & problem==pars$task_id)
                    lambda_refs = lambda_refs[order(lambda_refs$lambda, decreasing=TRUE),]
                    below_idx = max(which(lambda_refs$nzero_prop<pars$s))
                    above_idx = min(which(lambda_refs$nzero_prop>pars$s))
                    below_s = lambda_refs[below_idx,]$nzero_prop
                    above_s = lambda_refs[above_idx,]$nzero_prop
                    below_lambda = lambda_refs[below_idx,]$lambda
                    above_lambda = lambda_refs[above_idx,]$lambda
                    lambda = below_lambda + (pars$s - below_s) *
                        ((above_lambda - below_lambda) / (above_s - below_s))
                }
            } else {
                lambda = pars$lambda
            }
            if (is.null(lambda)) lambda = rep(0, length(x))
            private$iclust2_args = list(
              K=pars$K,
              lambda=rep(lambda, length(x)),
              chr=NULL,
              method=rep("lasso", length(x)),
              maxiter=pars$maxiter,
              eps=pars$eps,
              eps2=pars$eps2
            )
        }

        #Remove bad columns
        n = names(x)
        cols = lapply(x, function(xi) apply(xi, 2, stats::var)!=0)
        x = lapply(names(x), function(n) x[[n]][,cols[[n]]])

        means = lapply(x, \(block) apply(block, 2, mean))
        x = lapply(seq_along(x),
            function(i) scale(x[[i]], center=means[[i]], scale=FALSE))

        names(x) = n

        iclust2_fit = mlr3misc::invoke(iClusterPlus::iCluster2,
          .args=c(list(x=x), private$iclust2_args))
        message(paste0("iter: ", iclust2_fit$iter, " dif: ", iclust2_fit$dif))

        return(list(x=t(iclust2_fit$meanZ), jdr=list(beta=iclust2_fit$beta,
            phi=iclust2_fit$Phivec, cols=cols, means=means,
            x=t(iclust2_fit$meanZ))))
    },

    predict_jdr = function(newx, jdr, pars) {

        message("predict")
        n = names(newx)
        newx = lapply(names(newx),
                      function(n) newx[[n]][,jdr$cols[[n]]])
        newx = lapply(seq_along(newx),
            function(i) scale(newx[[i]], center=jdr$means[[i]], scale=FALSE))
        names(newx) = n

        newx = Reduce(cbind, newx)
        B = as.matrix(Reduce(rbind, jdr$beta))
        btp = B/jdr$phi
        btpb = t(btp) %*% B
        btpb = btpb + diag(nrow(btpb))
        tempm1 = btp %*% solve(btpb)
        #latent_space = t(tempm1) %*% t(newx)
        latent_space = newx %*% tempm1
        if(nrow(jdr$x)==nrow(latent_space)) {
            print(cor(latent_space, jdr$x))
            print(apply(abs(latent_space - jdr$x), 2, mean))
            #plot(latent_space[,1], jdr$x[,1])
            abline(a=0, b=1)
        }

        return(list(x=latent_space))
    }

  )
)

.extralrns_dict$add("surv.icluster2", LearneriCluster2)
