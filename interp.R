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
    bcs <- interp( x = y, z = v, xo = x0, yo = y0, 
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