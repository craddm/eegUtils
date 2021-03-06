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

#' Convert to 3d matrix
#'
#' @param data data to be converted
#' @param ... additional parameters
#' @keywords internal
conv_to_mat <- function(data,...) {
  UseMethod("conv_to_mat", data)
}

conv_to_mat.default <- function(data, ...) {
  stop("Not implemented for objects of class", class(data))
}

#' @describeIn conv_to_mat Convert eeg_epochs to 3D matrix
conv_to_mat.eeg_epochs <- function(data, ...) {
  n_epochs <- length(unique(data$timings$epoch))
  n_channels <- ncol(data$signals)
  n_times <- length(unique(data$timings$time))
  data <- array(as.matrix(data$signals),
                dim = c(n_times, n_epochs, n_channels))
  data
}


is.true <- function(x) {
  !is.na(x) & x
}

#' Generate a zero centred vector
#'
#' @param vec_length how many points you want the vector to have
#' @importFrom stats median
#' @keywords internal
zero_vec <- function(vec_length) {
  out_seq <- seq(1, vec_length, by = 1)
  out_seq - stats::median(out_seq)
}

#' Pad a vector with zeros
#'
#' Pads a vector with zeros at the beginning and end.
#'
#' @param x vector to pad
#' @param n number of zeros to pad
#' @param startval value to add at the start of the vector - defaults to zero
#' @param endval value to add at the end of the vector - defaults to zero
#' @keywords internal
pad <- function(x,
                n,
                startval = 0,
                endval = 0) {c(rep(startval, n),
                         x,
                         rep(endval, n))}

#' Unpad a vector
#'
#' remove a specified number of points from the beginning and end of a vector
#' @param x Vector to unpad/trim
#' @param n Number of points to remove
#' @keywords internal
unpad <- function(x, n) {
  start <- n + 1
  end <- length(x) - n
  x <- x[start:end]
  x
}

#' Fix group delay
#'
#' Corrects a signal for the group delay of an FIR filter.
#'
#' @keywords internal
fix_grpdelay <- function(x, n, grp_delay) {
  start <- n + 1 + grp_delay
  end <- length(x) - n + grp_delay
  x <- x[start:end]
  x
}

#' Get number of samples
#' @param .data Object to get total number of sampling points from
#' @keywords internal
samples <- function(.data) {
  nrow(.data$signals)
}

#' Generate circle as radians
#' @keywords internal
circ_rad_fun <- function() {
  seq(0,
      2 * pi,
      length.out = 101)
}

#' Update radius
#' @keywords internal
update_r <-
  function(r = 95,
           data,
           interp_limit) {

    max_elec <- calc_max_elec(data)
    r <- switch(interp_limit,
                "head" = min(max_elec * 1.10, max_elec + 15),
                "skirt" = r) # mm are expected for coords, 95 is good approx for Fpz - Oz radius
    r
  }

#' Calculate maximum electrode distance from origin.
#' @keywords internal
calc_max_elec <- function(data) max(sqrt(data$x^2 + data$y^2), na.rm = TRUE)
