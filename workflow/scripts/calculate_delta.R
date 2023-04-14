library( ape )
library( Biostrings )
library( readr )
library( R.utils )
library("expm")
source( "workflow/scripts/delta.R" )

args <- commandArgs(asValue=TRUE, excludeReserved=TRUE, defaults=list( shuffle=FALSE ) )[-1]

ma <- readDNAMultipleAlignment( filepath=args$alignment, format="fasta" )
ma <- as.matrix( ma, use.names=TRUE )

tree <- read.tree( args$tree )
tree$edge.length[tree$edge.length==0] <- quantile(tree$edge.length,0.1)*0.1
browser()
positions = read_table( args$positions, col_names=FALSE )

delta_wrapper <- function( position ) {
    char_vector <- as.vector( ma[tree$tip.label,position] )
    suppressWarnings( delta_result <- delta( char_vector, tree, 0.1, 0.0589, 10000, 10, 100 ) )
    print( paste( "Finished with", position ) )
    return( delta_result )
}

#print( positions )
#apply( positions, 2, print )
results <- apply( positions, 1, delta_wrapper )
results <- data.frame( positions=positions, delta=results )
colnames( results ) <- c( "position", "delta" )
write_csv( results, args$output )