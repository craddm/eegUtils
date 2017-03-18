#' Topographical Plotting Function for EEG
#'
#' Allows simple plotting of functional data
#' @param df An EEG dataset. Must have columns x, y, and amplitude at present. x and y are (Cartesian) electrode co-ordinates), amplitude is amplitude. Only supply one timepoint, or it will crash.
#' @param electrodeLocs Not yet implemented.
#' @param method Interpolation method. Defaults to Biharmonic spline interpolation, which is the only one currently available.
#' @param rmax Maximum head radius if x is too small.
#' @param gridRes Resolution of the interpolated grid. Higher = smoother but slower.
#' @import ggplot2
#' @import dplyr

topoplot <- function(df, electrodeLocs = NULL, method = "Biharmonic", rmax = .75, gridRes = 67) {

  # Generate a grid to interpolate data to ---------
  tmp <- df

  df <- data.frame(x = tmp$x, y = tmp$y, z = tmp$amplitude)

  if("x" %in% colnames(df)){
    xo <- seq(min(-rmax, df$x), max(rmax, df$x), length = gridRes)
  } else {
    df <- left_join(df,electrodeLocs,by = "electrode")
    xo <- seq(min(-rmax, df$x), max(rmax, df$x), length = gridRes)
  }
  if("y" %in% colnames(df)){
    yo <- seq(max(rmax, df$y), min(-rmax, df$y), length = gridRes)
  }

  r <- round(max(df$x)) / 2;
  tt <- seq(0, 2 * pi, length.out = 100)
  headShape <- data.frame(x = r * cos(tt),
                           y = r * sin(tt))
  maskRing <- data.frame(x = (1.42/2) * cos(tt),
                         y = (1.42/2) * sin(tt)
                         )

  nose <- data.frame(x = c(-0.05, 0, .05), y = c(.495, .55, .495))

  xo <- matrix(rep(xo, length(yo)), nrow = length(xo), ncol = length(yo))
  yo <- t(matrix(rep(yo, length(xo)), nrow = length(yo), ncol = length(xo)))
  xy <- df$x + df$y * sqrt(as.complex(-1))
  d <- matrix(rep(xy, length(xy)), nrow = length(xy), ncol = length(xy))
  d <- abs(d - t(d))
  diag(d) <- 1
  g <- (d ^ 2) * (log(d)-1)   # Green's function.
  diag(g) <- 0
  weights <- qr.solve(g, df$z)
  xy <- t(xy)
  outmat <- matrix(nrow = gridRes, ncol = gridRes)
  for (i in 1:gridRes){
    for (j in 1:gridRes) {
      test4 <- abs((xo[i,j] + sqrt(as.complex(-1))*yo[i,j]) - xy)
      zCheck <- which(test4 == 0)
      if (length(zCheck) > 0){
        test4[zCheck] <- 1
      }
      g <- (test4^2) * (log(test4)-1)
      if (length(zCheck) > 0){
        g[zCheck] <- 0
      }
      outmat[i,j] <- g %*% weights
    }
  }
  outDf <- data.frame(x = xo[,1],outmat)
  names(outDf)[1:length(yo[1,])+1] <- yo[1,]

  outDf <- gather(outDf,
                key = y,
                value = amplitude,
                -x,
                convert = TRUE)

  outDf$incircle <- sqrt(outDf$x^2 + outDf$y^2) < .7

  ggplot(outDf[outDf$incircle, ], aes(x, y, fill = amplitude ))+
    geom_raster(interpolate = TRUE)+
    stat_contour(aes(z = amplitude,linetype = ..level..<0),
                 bins = 6,
                 colour = "black",
                 size = 1.2,
                 show.legend = FALSE
    )+
    scale_fill_distiller(palette = "RdBu",
                         guide = "colourbar",
                         oob = scales::squish)+
    annotate("path", x = maskRing$x, y = maskRing$y, colour = "white", size = 6)+
    geom_point(data = df,
               aes(x, y, z = NULL, fill = NULL),
               size = 1)+
    annotate("path", x = headShape$x, y = headShape$y,
              size = 1)+
    annotate("path", x = nose$x, y = nose$y,
              size = 1)+
    coord_equal()+
    theme_bw() +
    theme(rect = element_blank(),
      line = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank()
    )
}
