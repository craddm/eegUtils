#' Modified ggplot2 code for contours
#'
#' This file contains modified versions of source code from ggplot2. As the
#' `geom_contour` code relies on several internal ggplot2 functions, I have
#' reproduced them here to be able to modify the `geom_contour` code
#' appropriately.
#'
#' The original code was released under an MIT license, and remains copyright of
#' the ggplot2 authors.
NULL

#' @noRd
stat_summary_by_z <- function(mapping = NULL, data = NULL,
                              geom = "contour", position = "identity",
                              ...,
                              bins = NULL,
                              binwidth = NULL,
                              breaks = NULL,
                              na.rm = FALSE,
                              show.legend = NA,
                              inherit.aes = TRUE) {
  ggplot2::layer(
    data = data,
    mapping = mapping,
    stat = StatSummarybyZ,
    geom = geom,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      bins = bins,
      binwidth = binwidth,
      breaks = breaks,
      na.rm = na.rm,
      ...
    )
  )
}

StatSummarybyZ <- ggplot2::ggproto("StatSummaryByZ", Stat,

                       required_aes = c("x", "y", "z"),
                       default_aes = aes(order = after_stat(level)),

                       setup_params = function(data, params) {

                         params$z.range <- range(data$z,
                                                 na.rm = TRUE,
                                                 finite = TRUE)
                         params
                       },

                       compute_group = function(data, scales,
                                                z.range, bins = NULL,
                                                binwidth = NULL,
                                                breaks = NULL, na.rm = FALSE) {

                         data <-
                           aggregate(z ~ x + y,
                                     data = data,
                                     FUN = mean,
                                     na.rm = na.rm,
                                     na.action = na.pass)

                         breaks <- contour_breaks(z.range, bins,
                                                  binwidth, breaks)

                         isolines <- xyz_to_isolines(data, breaks)

                         path_df <- iso_to_path(isolines,
                                                data$group[1])

                         path_df$level <- as.numeric(path_df$level)
                         path_df$nlevel <- scales::rescale_max(path_df$level)

                         path_df
                       }
)

#' Create an interpolated scalp surface
#'
#' `stat_scalpcontours` creates an interpolated surface for an irregular set of
#' x-y coordinates, as is typically required for a topographical EEG plot, and
#' then calculates contours.
#'
#' @inheritParams ggplot2::geom_contour
#' @author Matt Craddock \email{matt@@mattcraddock.com}
#' @param grid_res Resolution of the interpolation grid. (Defaults to 200
#'   points).
#' @param interp_limit Interpolate to the "skirt" or inside the "head".
#' @param method "biharmonic" or "gam" smooth to create interpolated surface
#' @param r Radius of interpolated surface
#' @param bins Number of contour bins. Overridden by binwidth.
#' @param binwidth The width of the contour bins. Overridden by breaks.
#' @param breaks Numeric vector to set the contour breaks. Overrides binwidth
#'   and bins. By default, this is a vector of length ten with pretty() breaks.
#' @family topoplot functions
#' @export
stat_scalpcontours <- function(mapping = NULL,
                               data = NULL,
                               position = "identity",
                               na.rm = FALSE,
                               show.legend = FALSE,
                               inherit.aes = TRUE,
                               grid_res = 200,
                               interp_limit = c("skirt", "head"),
                               method = "biharmonic",
                               r = NULL,
                               bins = 6,
                               binwidth = NULL,
                               breaks = NULL,
                               ...) {
  ggplot2::layer(
    stat = StatScalpContours,
    data = data,
    mapping = mapping,
    geom = ggplot2::GeomContour,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(na.rm = na.rm,
                  grid_res = grid_res,
                  interp_limit = interp_limit,
                  method = method,
                  r = r,
                  bins = bins,
                  breaks = breaks,
                  binwidth = binwidth,
                  ...)
  )
}

#' StatScalpcontours
#'
#' @export
#' @keywords internal

StatScalpContours <-
  ggplot2::ggproto("StatScalpContours",
                   Stat,
                   required_aes = c("x",
                                    "y",
                                    "z"),
                   # if no z provided, can add fill as default, but then get a
                   # different warning
                   default_aes = aes(order = after_stat(level),
                                     linetype = ggplot2::after_stat(level) < 0),

                   setup_params = function(data, params) {

                     params$z.range <- range(data$z,
                                             na.rm = TRUE,
                                             finite = TRUE)
                     params
                   },

                   compute_group = function(data,
                                            scales,
                                            method,
                                            z.range,
                                            bins = NULL,
                                            interp_limit,
                                            grid_res,
                                            r,
                                            na.rm = FALSE,
                                            breaks = NULL,
                                            binwidth = NULL) {

                     interp_limit <- match.arg(interp_limit,
                                               c("skirt", "head"))

                     data <- aggregate(fill ~ x + y,
                                       data = data,
                                       FUN = mean,
                                       na.rm = TRUE,
                                       na.action = na.pass)

                     if (is.null(r)) {
                       r <- update_r(95,
                                     data = data,
                                     interp_limit = interp_limit)
                     }

                     if (identical(method, "biharmonic")) {
                       data <- biharmonic(data,
                                          grid_res = grid_res,
                                          interp_limit = interp_limit,
                                          r = r)
                     } else {
                       data <- fit_gam_topo(data,
                                            grid_res = grid_res,
                                            interp_limit = interp_limit,
                                            r = r)
                     }

                     data <- dplyr::rename(data,
                                           z = fill)

                     z.range <- range(data$z,
                                      na.rm = TRUE,
                                      finite = TRUE)

                     breaks <- contour_breaks(z.range, bins,
                                              binwidth, breaks)

                     isolines <- xyz_to_isolines(data,
                                                 breaks)
                     path_df <- iso_to_path(isolines)#,
                                            #data$group[1])

                     path_df$level <- as.numeric(path_df$level)
                     path_df$nlevel <- scales::rescale_max(path_df$level)

                     path_df
                   }
  )

#' Calculate the breaks used for contouring
#'
#' @inheritParams geom_contour
#' @param z_range Range of values within which breaks should be calculated
#'
#' @return A vector of breaks
#' @noRd
#'
contour_breaks <- function(z_range, bins = NULL, binwidth = NULL, breaks = NULL) {
  if (!is.null(breaks)) {
    return(breaks)
  }

  # If no parameters set, use pretty bins
  if (is.null(bins) && is.null(binwidth)) {
    breaks <- pretty(z_range, 10)
    return(breaks)
  }

  # If provided, use bins to calculate binwidth
  if (!is.null(bins)) {
    # round lower limit down and upper limit up to make sure
    # we generate bins that span the data range nicely
    accuracy <- signif(diff(z_range), 1)/10
    z_range[1] <- floor(z_range[1]/accuracy)*accuracy
    z_range[2] <- ceiling(z_range[2]/accuracy)*accuracy

    if (bins == 1) {
      return(z_range)
    }

    binwidth <- diff(z_range) / (bins - 1)
    breaks <- scales::fullseq(z_range, binwidth)

    # Sometimes the above sequence yields one bin too few.
    # If this happens, try again.
    if (length(breaks) < bins + 1) {
      binwidth <- diff(z_range) / bins
      breaks <- scales::fullseq(z_range, binwidth)
    }

    return(breaks)
  }

  # if we haven't returned yet, compute breaks from binwidth
  scales::fullseq(z_range, binwidth)
}

#' Compute isoband objects
#'
#' @param data A data frame with columns `x`, `y`, and `z`.
#' @param breaks A vector of breaks. These are the values for
#'   which contour lines will be computed.
#'
#' @return An S3 "iso" object, which is a `list()` of `list(x, y, id)`s.
#' @noRd
#'
xyz_to_isolines <- function(data, breaks) {
  isoband::isolines(
    x = sort(unique(data$x)),
    y = sort(unique(data$y)),
    z = isoband_z_matrix(data),
    levels = breaks
  )
}

xyz_to_isobands <- function(data, breaks) {
  isoband::isobands(
    x = sort(unique(data$x)),
    y = sort(unique(data$y)),
    z = isoband_z_matrix(data),
    levels_low = breaks[-length(breaks)],
    levels_high = breaks[-1]
  )
}

#' Compute input matrix for isoband functions
#'
#' Note that [grDevices::contourLines()] needs transposed
#' output to the matrix returned by this function.
#'
#' @param data A data frame with columns `x`, `y`, and `z`.
#'
#' @return A [matrix()]
#' @noRd
#'
isoband_z_matrix <- function(data) {
  # Convert vector of data to raster
  x_pos <- as.integer(factor(data$x, levels = sort(unique(data$x))))
  y_pos <- as.integer(factor(data$y, levels = sort(unique(data$y))))

  nrow <- max(y_pos)
  ncol <- max(x_pos)

  raster <- matrix(NA_real_, nrow = nrow, ncol = ncol)
  raster[cbind(y_pos, x_pos)] <- data$z

  raster
}

#' Convert the output of isolines functions
#'
#' @param iso the output of [isoband::isolines()]
#' @param group the name of the group
#'
#' @return A data frame that can be passed to [geom_path()].
#' @noRd
#'
iso_to_path <- function(iso, group = 1) {
  lengths <- vapply(iso, function(x) length(x$x), integer(1))

  if (all(lengths == 0)) {
    rlang::warn("stat_contour(): Zero contours were generated")
    return(new_data_frame())
  }

  levels <- names(iso)
  xs <- unlist(lapply(iso, "[[", "x"), use.names = FALSE)
  ys <- unlist(lapply(iso, "[[", "y"), use.names = FALSE)
  ids <- unlist(lapply(iso, "[[", "id"), use.names = FALSE)
  item_id <- rep(seq_along(iso), lengths)

  # Add leading zeros so that groups can be properly sorted
  groups <- paste(group, sprintf("%03d", item_id), sprintf("%03d", ids), sep = "-")
  groups <- factor(groups)

  new_data_frame(
    list(
      level = rep(levels, lengths),
      x = xs,
      y = ys,
      piece = as.integer(groups),
      group = groups
    ),
    n = length(xs)
  )
}

#' Convert the output of isoband functions
#'
#' @param iso the output of [isoband::isobands()]
#' @param group the name of the group
#'
#' @return A data frame that can be passed to [geom_polygon()].
#' @noRd
#'
iso_to_polygon <- function(iso, group = 1) {
  lengths <- vapply(iso, function(x) length(x$x), integer(1))

  if (all(lengths == 0)) {
    rlang::warn("stat_contour(): Zero contours were generated")
    return(new_data_frame())
  }

  levels <- names(iso)
  xs <- unlist(lapply(iso, "[[", "x"), use.names = FALSE)
  ys <- unlist(lapply(iso, "[[", "y"), use.names = FALSE)
  ids <- unlist(lapply(iso, "[[", "id"), use.names = FALSE)
  item_id <- rep(seq_along(iso), lengths)

  # Add leading zeros so that groups can be properly sorted
  groups <- paste(group, sprintf("%03d", item_id), sep = "-")
  groups <- factor(groups)

  new_data_frame(
    list(
      level = rep(levels, lengths),
      x = xs,
      y = ys,
      piece = as.integer(groups),
      group = groups,
      subgroup = ids
    ),
    n = length(xs)
  )
}

#' Pretty isoband level names
#'
#' @param isoband_levels `names()` of an [isoband::isobands()] object.
#'
#' @return A vector of labels like those used in
#'   [cut()] and [cut_inverval()].
#' @noRd
#'
pretty_isoband_levels <- function(isoband_levels, dig.lab = 3) {
  interval_low <- gsub(":.*$", "", isoband_levels)
  interval_high <- gsub("^[^:]*:", "", isoband_levels)

  label_low <- format(as.numeric(interval_low), digits = dig.lab, trim = TRUE)
  label_high <- format(as.numeric(interval_high), digits = dig.lab, trim = TRUE)

  # from the isoband::isobands() docs:
  # the intervals specifying isobands are closed at their lower boundary
  # and open at their upper boundary
  sprintf("(%s, %s]", label_low, label_high)
}

# Fast data.frame constructor and indexing
# No checking, recycling etc. unless asked for
new_data_frame <- function(x = list(), n = NULL) {
  if (length(x) != 0 && is.null(names(x))) {
    rlang::abort("Elements must be named")
  }
  lengths <- vapply(x, length, integer(1))
  if (is.null(n)) {
    n <- if (length(x) == 0 || min(lengths) == 0) 0 else max(lengths)
  }
  for (i in seq_along(x)) {
    if (lengths[i] == n) next
    if (lengths[i] != 1) {
      rlang::abort("Elements must equal the number of rows or 1")
    }
    x[[i]] <- rep(x[[i]], n)
  }

  class(x) <- "data.frame"

  attr(x, "row.names") <- .set_row_names(n)
  x
}


