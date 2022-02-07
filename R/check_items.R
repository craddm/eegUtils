#' Check that all classes in a list match
#'
#' @param data list of objects to check
#' @keywords internal

check_classes <- function(data) {

  stopifnot(is.list(data))

  dat_classes <- lapply(data,
                        class)
  check_class <- sapply(dat_classes,
                        identical,
                        dat_classes[[1]])
  all(check_class)
}

#' Check that all objects have the same number of trials
#'
#' @param data list of objects to check
#' @keywords internal

check_epochs <- function(data) {

  stopifnot(is.list(data))

  dat_epochs <- lapply(data,
                        function(x) nrow(x$signals))
  check_epochs <- sapply(dat_epochs,
                        identical,
                        dat_epochs[[1]])
  all(check_epochs)
}


#' Check that all conditions in a list match
#' @noRd

check_conds <- function(data_list) {

  get_names <- lapply(data_list,
                      function(x) names(x$signals))
  check_names <- sapply(get_names,
                        identical,
                        get_names[[1]])
  all(check_names)
}

#' Check consistency of labels
#'
#' Internal function for checking 1) whether the labels submitted are a mixture
#' of hierarchical and non-hierarchical types 2) whether the labels submitted
#' are present in the data
#'
#' @author Matt Craddock \email{matt@@mattcraddock.com}
#' @param cond_labs labels submitted by the user
#' @param data_labs labels from the actual data
#' @keywords internal

label_check <- function(cond_labs,
                        data_labs) {

  if (all(grepl("/", cond_labs))) {
    lab_check <- cond_labs %in% data_labs
  } else if (any(grepl("/",
                       cond_labs))) {
    stop("Do not mix hierarchical and non-hierarchical event labels.")
  } else {
    # Check if there is a hierarchical separator "/". If so,
    # split the labels
    if (any(grepl("/",
                  data_labs))) {
      split_labels <- strsplit(data_labs,
                               "/")

      lab_check <- lapply(cond_labs,
                          function(x) vapply(split_labels,
                                             function(i) x %in% i,
                                             logical(1)))
      #condense to a single TRUE or FALSE for each label
      lab_check <- vapply(lab_check,
                          any,
                          logical(1))
    } else {
      lab_check <- cond_labs %in% data_labs
    }
  }
}
