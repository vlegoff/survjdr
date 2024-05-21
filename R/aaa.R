.extralrns_dict = R6::R6Class(".extralrns_dict",
  public = list(
    lrns = list(),
    add = function(key, learn) {
      checkmate::assert_character(key, len = 1)
      lst = list(key = learn)
      names(lst) = key
      self$lrns = mlr3misc::insert_named(self$lrns, lst)
    }
  )
)$new()
