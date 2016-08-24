pointWulff <- function( theta, phi, R = 1 ) {
  r = R * tan( pi / 4 - theta / 2 )
  x = r * sin( phi )
  y = r * cos( phi )
  points( x, y, pch = 16 )  
  invisible( )
}

pointSchmidt <- function( theta, phi, R = 1 ) {
  theta <- theta * pi / 180.0
  phi <- phi * pi / 180.0
  r = R * sqrt( 2 ) * cos( theta / 2 + pi / 4 )
  x = r * sin( phi )
  y = r * cos( phi )
  points( x, y, pch = 16 )    
  invisible( )
}

gcWulff <- function( strike, dip, R = 1 ) {
    theta = strike - 90
    theta <- theta * pi / 180.0
    dip <- dip * pi / 180.0
    c <- cos( -theta )
    s <- sin( -theta )
    rot <- matrix( c( c, s, -s, c ), nrow = 2 )
    beta <- seq( pi, 0, length.out = 90 )
    x = R / cos( dip ) * sin( beta )
    y = R * tan( dip ) + R / cos( dip ) * cos( beta )
    X <- cbind( x, y ) %*% t( rot )
    l <- x^2 + y^2
    ll <- which( l <= R )
    print( lines( X[ll,] ) )
    
    beta <- seq( pi, 2 * pi, length.out = 90 )
    x = R / cos( dip ) * sin( beta )
    y = R * tan( dip ) + R / cos( dip ) * cos( beta )
    l <- x^2 + y^2
    ll <- which( l <= R )
    X <- cbind( x, y ) %*% t( rot )
    print( lines( X[ll,] )  )
}


gcSchmidt <- function( strike, dip, R = 1 ) {
    theta = strike
    theta <- theta * pi / 180.0
    dip <- dip * pi / 180.0
    c <- cos( -theta )
    s <- sin( -theta )
    rot <- matrix( c( c, s, -s, c ), nrow = 2 )
    phi = seq( from = -90, to = 90, by = 5 ) * pi / 180
    lam0 <- pi / 2
    R = R * sqrt( 2 )/2
    kp = sqrt( 2 / ( 1 + cos( phi ) * cos( dip - lam0 ) ) )
    x = -R * kp * cos( phi ) * sin( dip - lam0 )
    y = R * kp * sin( phi )
    X <- cbind( x, y ) %*% t( rot )
    lines( X )
    invisible()
}