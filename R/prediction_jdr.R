#' @export
PredictionJDR = R6::R6Class("PredictionJDR",
  inherit = mlr3proba::PredictionSurv,
  public = list(
    initialize = function(task = NULL, row_ids = task$row_ids, 
      truth = task$truth(), crank = NULL, distr = NULL, lp = NULL,
      response = NULL, lambdas = NULL, best_lambda = NULL, 
      cv_grid = NULL, cv_results = NULL, check = TRUE) {

      private$.best_lambda = best_lambda
      private$.cv_grid = cv_grid
      private$.cv_results = cv_grid
      super$initialize(
        task = task,
        row_ids = row_ids,
        truth = truth,
        crank = crank,
        distr = distr,
        lp = lp,
        response = response,
        check = TRUE
      )
    }
  ),

  active = list(
    lambdas = function() {
      #self$.lambdas %??% rep(NA_real_, length(self$data$row_ids))
      self$.lambdas
    },
    best_lambda = function() {
      #self$.lambdas %??% rep(NA_real_, length(self$data$row_ids))
      self$.best_lambdas
    },
    cv_grid = function() {
      self$.cv_grid
    },
    cv_results = function() {
      self$.cv_results
    }
  ),

  private = list(
    .lambdas = NULL,
    .best_lambda = NULL,
    .cv_grid = NULL,
    .cv_results = NULL 
  )
)

#' @export
as_prediction_surv = function(x, ...) {
  UseMethod("as_prediction_surv")
}

#' @export
as_prediction_surv.PredictionJDR = function(x, ...) {
  x
}

#' @export
as_prediction.PredictionJDR = function(x, ...) {
  x
}
