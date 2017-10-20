# run with
# Rscript local.gradient.distance.r -i input_distmat.txt -o input_distmat.txt

library('optparse')

# make option list and parse command line
option_list <- list(
    make_option(c("-i","--input"), type="character",
        help="Path to input distance matrix [required]."),
    make_option(c("-n","--neighborhood_size"), type="numeric",
        default=NULL,
        help="Neighborhood size in which to trust distances. Smaller is better for detrending. [default: find minimum neighborhood size for connected graph]."),
    make_option(c("-o","--output"), type="character",
        help="Path to output distance matrix [required].")
)
opts <- parse_args(OptionParser(option_list=option_list), 
    args=commandArgs(trailing=TRUE))


# Note: if neighborhood size results in multiple connected components,
# greedily adds bridges between connected components until graph 
# is fully connected.
# 
"lg.dist" <- function(d,neighborhood.size=NULL, weighted=TRUE) {
    if(is.null(neighborhood.size)){
        neighborhood.sizes <- 3:10
    } else {
        neighborhood.sizes <- neighborhood.size
    }

    ix <- 1
    is.valid <- FALSE

    # Calculate number of eigenvalues == 0 in graph;
    # This is the number of connected components
    # If > 1, greedily bridge components till all connected
    # note: needs to be implemented;
    # right now this just increases neighborhood size
    # note: also need to throw a warning if a bridge added
    # is more than twice the length of a local edge
    while(ix <= length(neighborhood.sizes) && !is.valid){
        ns <- neighborhood.sizes[ix]
        g <- lg.graph(d, neighborhood.size=ns, weighted=weighted)
        eigs <- eigen(graph.laplacian(g),only.values=TRUE)$values
        is.valid <- sum(eigs < 10 * .Machine$double.eps) == 1
        # is.valid <- all(is.finite(lgd))
        ix <- ix + 1
    }
    lgd <- shortest.paths(g)

    if(!is.valid) return(NULL)
    return(lgd)
}


"lg.graph" <- function(d,neighborhood.size=4,weighted=TRUE) {
    require('igraph', warn.conflicts=FALSE, quietly=TRUE)
    d <- as.matrix(d)
    dd <- matrix(0,nrow(d), nrow(d))
    for(i in 1:nrow(dd)) dd[i,order(d[i,])[1:(neighborhood.size+1)]] <- 1
    if(weighted) dd[dd>0] <- d[dd>0] else weighted <- NULL

    g <- graph.adjacency(dd,weighted=weighted, mode='undirected')
    return(g)
}


d <- as.matrix(read.table(opts$input),sep='\t',head=T,row=1)
cat('Calculating local gradient distance...\n')
lgd <- lg.dist(d,neighborhood.size=opts$neighborhood_size,weighted=TRUE)

cat('Calculating PCoA of original distances...\n')
pc.d <- cmdscale(d)
cat('Calculating PCoA of transformed distances...\n')

if(is.null(lgd)){
    stop("Error: graph was not connected with given neighborhood size. Try a larger neighborhood size or automatic selection.")
}

pc.lgd <- cmdscale(lgd)
pdf('d-v-lgd.pdf',width=9,height=5)
par(mfrow=c(1,2))
plot(pc.d, xlim=range(pc.d), ylim=range(pc.d))
plot(pc.lgd, xlim=range(pc.lgd), ylim=range(pc.lgd))
dev.off()

sink(opts$output)
write.table(lgd, sep='\t',quote=F)
sink(NULL)