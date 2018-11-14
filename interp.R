# interp.R
#
# Set of 2D spatial interpolators.  Intended for use on SpatialPointsDataFrame data sets.
# Require input data, output grid (specified by a SpatialPoints object), and some
# routine-dependent parameters (e.g., linear or bicubic).

require(geometry)
require(Matrix)
require(akima)
require(sp)
require(SpatialTools)

akima_interp <- function (x, y = NULL, z, xo = seq(min(x), max(x), length = nx),
                          yo = seq(min(y), max(y), length = ny), linear = TRUE, extrap = FALSE,
                          duplicate = "error", dupfun = NULL, nx = 40, ny = 40, jitter = 10^-12,
                          jitter.iter = 6, jitter.random = FALSE)
{
  is.sp <- FALSE
  sp.coord <- NULL
  sp.z <- NULL
  sp.proj4string <- NULL
  if (is.null(y) && is.character(z)) {
    if (class(x) == "SpatialPointsDataFrame") {
      sp.coord <- dimnames(coordinates(x))[[2]]
      sp.z <- z
      sp.proj4string <- x@proj4string
      z <- x@data[, z]
      y <- coordinates(x)[, 2]
      x <- coordinates(x)[, 1]
      is.sp <- TRUE
      xo = seq(min(x), max(x), length = nx)
      yo = seq(min(y), max(y), length = ny)
    }
    else stop("either x,y,z are numerical or x is SpatialPointsDataFrame and z a name of a data column in x")
  }
  if (linear)
    ret <- interp.old(x, y, z, xo, yo, ncp = 0, extrap = FALSE,
                      duplicate = duplicate, dupfun = dupfun)
  else {
    if (!(all(is.finite(x)) && all(is.finite(y)) && all(is.finite(z))))
      stop("missing values and Infs not allowed")
    drx <- diff(range(x))
    dry <- diff(range(y))
    if (drx == 0 || dry == 0)
      stop("all data collinear")
    if (drx/dry > 10000 || drx/dry < 1e-04)
      stop("scales of x and y are too dissimilar")
    n <- length(x)
    nx <- length(xo)
    ny <- length(yo)
    if (length(y) != n || length(z) != n)
      stop("Lengths of x, y, and z do not match")
    xy <- paste(x, y, sep = ",")
    if (duplicate == "error") {
      if (any(duplicated(xy)))
        stop("duplicate data points: need to set 'duplicate = ..' ")
    }
    else {
      i <- match(xy, xy)
      if (duplicate == "user")
        dupfun <- match.fun(dupfun)
      ord <- !duplicated(xy)
      if (duplicate != "strip") {
        centre <- function(x) switch(duplicate, mean = mean(x),
                                     median = median(x), user = dupfun(x))
        z <- unlist(lapply(split(z, i), centre))
      }
      else {
        z <- z[ord]
      }
      x <- x[ord]
      y <- y[ord]
      n <- length(x)
    }
    miss <- !extrap
    ans <- .Fortran("sdsf3p", md = as.integer(1), ndp = as.integer(n),
                    xd = as.double(x), yd = as.double(y), zd = as.double(z),
                    nx = as.integer(nx), x = as.double(xo), ny = as.integer(ny),
                    y = as.double(yo), z = as.double(matrix(0, nx, ny)),
                    ier = integer(1), wk = double(36 * n), iwk = integer(25 *
                                                                           n), extrap = as.logical(matrix(extrap, nx, ny)),
                    near = integer(n), nxt = integer(n), dist = double(n),
                    linear = as.logical(linear), PACKAGE = "akima")
    if (miss)
      ans$z[ans$extrap] <- NA
    if (ans$ier == 10) {
      warning("collinear points, trying to add some jitter to avoid colinearities!")
      jitter.trials <- 1
      success <- FALSE
      while (jitter.trials < jitter.iter & !success) {
        if (jitter.random) {
          j <- list()
          j[[1]] <- rep(c(-1, 0, 1), length.out = length(x))
          j[[2]] <- rep(c(0, 1, -1), length.out = length(x))
          j[[3]] <- rep(c(1, -1, 0), length.out = length(x))
          jx <- sample(1:3, 1)
          jy <- sample(1:3, 1)
          xj <- x + j[[jx]] * diff(range(x)) * jitter *
            jitter.trials^1.5
          yj <- y + j[[jy]] * diff(range(y)) * jitter *
            jitter.trials^1.5
        }
        else {
          xj <- x + rep(c(-1, 0, 1), length.out = length(x)) *
            diff(range(x)) * jitter * jitter.trials^1.5
          yj <- y + rep(c(0, 1, -1), length.out = length(y)) *
            diff(range(y)) * jitter * jitter.trials^1.5
        }
        ans <- .Fortran("sdsf3p", as.integer(1), as.integer(n),
                        xd = as.double(xj), yd = as.double(yj), as.double(z),
                        as.integer(nx), x = as.double(xo), as.integer(ny),
                        y = as.double(yo), z = as.double(matrix(0,
                                                                nx, ny)), ier = integer(1), double(36 * n),
                        integer(25 * n), extrap = as.logical(matrix(extrap,
                                                                    nx, ny)), near = integer(n), nxt = integer(n),
                        dist = double(n), linear = as.logical(linear),
                        PACKAGE = "akima")
        if (miss)
          ans$z[ans$extrap] <- NA
        if (linear)
          ans$z[ans$extrap] <- NA
        success <- (ans$ier == 0)
        if (success)
          warning("success: collinearities reduced through jitter")
        jitter.trials <- jitter.trials + 1
      }
    }
    if (is.sp) {
      zm <- ny
      zn <- nx
      zvec <- c(ans$z)
      xvec <- c(matrix(rep(ans$x, zn), nrow = zm, ncol = zn,
                       byrow = FALSE))
      yvec <- c(matrix(rep(ans$y, zm), nrow = zm, ncol = zn,
                       byrow = TRUE))
      nona <- !is.na(zvec)
      ret <- data.frame(xvec[nona], yvec[nona], zvec[nona])
      names(ret) <- c(sp.coord[1], sp.coord[2], sp.z)
      coordinates(ret) <- sp.coord
      ret@proj4string <- sp.proj4string
      gridded(ret) <- TRUE
    }
    else {
      ret <- list(x = ans$x, y = ans$y, z = matrix(ans$z,
                                                   nx, ny))
    }
  }
  ret
}


# Surface fit to enclosing Delaunay Triangle.  Does not extrapolate, points outside the convex
# hull have value NA.
interp.bary <- function( x, y, v = 1 ) {
    if ( class( y ) != "SpatialPointsDataFrame" ) stop( "Input must be SpatialPointsDataFrame." )
    if ( class( x ) != "SpatialPoints" ) stop( "Output grid must be SpatialPoints." )
    f <- y@data[,v]
    dn <- delaunayn( y@coords )
    tri <- tsearch( y@coords[,1], y@coords[,2], dn, x@coords[,1], x@coords[,2], bary = T )
    ll <- which( !is.na( tri$idx ) )
    m <- length( ll )
    active <- dn[tri$idx[ll],]
    M <- sparseMatrix(i=rep(1:m,each=3),j=as.numeric(t(active)),x=as.numeric(t(tri$p[ll,])),dims=c(m,length(f)))
    bi <- as.numeric( M %*% f )
    bn <- rep( NA, length( x@coords[,1] ) )
    bn[ll] <- bi
    return( bn )
}

# Bicubic or Bilinear Spline Interpolation
interp.bcspl <- function( x, y, v = 1, linear = FALSE ) {
    if ( class( y ) != "SpatialPointsDataFrame" ) stop( "Input must be SpatialPointsDataFrame." )
    if ( class( x ) != "SpatialPoints" ) stop( "Output grid must be SpatialPoints." )
    x0 <- unique( x@coords[,1] )
    y0 <- unique( x@coords[,2] )

    bcs <- akima_interp( x = y, z = v, xo = x0, yo = y0,
                   linear = linear, extrap = TRUE,
                   nx = length( x0 ), ny = length( y0 ) )

    return( bcs@data[,v] )
}


# Inverse Distance Interpolation
interp.invdist <- function( x, y, v = 1 ) {
    eps <- 0.000001
    if ( class( y ) != "SpatialPointsDataFrame" ) stop( "Input must be SpatialPointsDataFrame." )
    if ( class( x ) != "SpatialPoints" ) stop( "Output grid must be SpatialPoints." )
    dM <- dist2( x@coords, y@coords )
    ll <- which( dM < eps )
    dM[ll] <- eps
    W <- 1 / dM
    Wsum <- apply( W, 1, sum )
    val <- W %*% y@data[,v] / Wsum
    return( val )
}

# Inverse Exponential Weighting Interpolation
interp.invexp <- function( x, y, v = 1, sigma = 1 ) {
    if ( class( y ) != "SpatialPointsDataFrame" ) stop( "Input must be SpatialPointsDataFrame." )
    if ( class( x ) != "SpatialPoints" ) stop( "Output grid must be SpatialPoints." )
    dM <- dist2( x@coords, y@coords )
    W <- exp( -dM^2 / ( 2 * sigma^2 ) )
    Wsum <- apply( W, 1, sum )
    val <- W %*% y@data[,v] / Wsum
    return( val )
}

# Local Mean Value Interpolation (return mean value of n-closest points. If there is a tie
# then include all values with distances <= n-th largest distance)
interp.nn <- function( x, y, v = 1, n = 1 ) {
    if ( class( y ) != "SpatialPointsDataFrame" ) stop( "Input must be SpatialPointsDataFrame." )
    if ( class( x ) != "SpatialPoints" ) stop( "Output grid must be SpatialPoints." )
    dM <- dist2( x@coords, y@coords )
    dS <- apply( dM, 1, sort )
    lst <- t( dS[n,] )
    m <- length( x@coords[,1] )
    val <- vector( mode = "double", length = m )
    for ( i in 1:m ) {
        ll <- which( dM[i,] <= lst[i] )
        val[i] <- sum( y@data[ll,v] ) / length( ll )
    }
    return( val )
}
