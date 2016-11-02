# Tools for applying PCA to sets of time series.  Uses output of prcomp.


#pcaVar - Tabulate variance and cumulative variance explained by the principal components in a
#         prcomp analysis.  If plot = T, make a nice plot of the result.

pcaVar <- function( res.pca, plot = T, n = NULL ) {
    if ( is.null( n ) ) n <- length( res.pca$sdev )
    eig <- res.pca$sdev^2
    variance <- eig * 100 / sum( eig )
    cumvar <- cumsum( variance )
    var <- data.frame( PC = seq( 1, n ), Eigen = eig, Var = variance, CumVar = cumvar )
    if ( plot ) {
        g <- ggplot( data = var[1:n,], aes( PC, Var ) ) + geom_bar( stat = "identity", fill = "steelblue", color = "black" ) + 
             labs( x = "Eigenvector", y = "Percentage of Variance Explained", 
             title = "Variance Distribution" )
        print( g )
    }
    return( var )
}

# pcaPPC - Plots the PCs scaled by st.dev.

pcaPPC <- function( res.pca, which = 1:length( res.pca$sdev )  ) {
    if( is.null( which ) ) which <- 1:length( res.pca$sdev )
    n <- nrow( res.pca$x )
    t <- seq( 1: n )
    x.df <- data.frame( Time = t )
    for ( i in 1:length( which ) ) {
        s <- sprintf( "PC #%d", which[i] )
        t.df <- data.frame( res.pca$x[ , which[i]])
        colnames( t.df ) <- c( s )
        x.df <- cbind( x.df, t.df )
    }
    x.gr <- melt( x.df, id.var = c( "Time" ))
    g <- ggplot( data = x.gr, aes( x = Time ) ) + 
          geom_line( aes( y = value ) ) + 
          facet_wrap( "variable" ) +
          labs( x = "Time Index", y = "Value", title = "Principal Components" )
    return( g )
}
