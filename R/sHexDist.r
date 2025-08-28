#' Function to calculate distances between hexagons/rectangles in a 2D grid
#'
#' \code{sHexDist} is supposed to calculate euclidian distances between each pair of hexagons/rectangles in a 2D grid of input "sTopol" or "sMap" object. It returns a symmetric matrix containing pairwise distances. 
#'
#' @param sObj an object of class "sTopol" or "sInit" or "sMap"
#' @return 
#' \itemize{
#'  \item{\code{dist}: a symmetric matrix of nHex x nHex, containing pairwise distances, where nHex is the total number of hexagons/rectanges in the grid}
#' }
#' @note The return matrix has rows/columns ordered in the same order as the "coord" matrix of the input object does.
#' @export
#' @seealso \code{\link{sTopology}}, \code{\link{sInitial}}
#' @include sHexDist.r
#' @examples
#' # 1) generate an iid normal random matrix of 100x10 
#' data <- matrix( rnorm(100*10,mean=0,sd=1), nrow=100, ncol=10) 
#'
#' # 2) from this input matrix, determine nHex=5*sqrt(nrow(data))=50, 
#' # but it returns nHex=61, via "sHexGrid(nHex=50)", to make sure a supra-hexagonal grid
#' sTopol <- sTopology(data=data, lattice="hexa", shape="suprahex") 
#'
#' # 3) initialise the codebook matrix using "uniform" method
#' sI <- sInitial(data=data, sTopol=sTopol, init="uniform") 
#'
#' # 4) calculate distances between hexagons/rectangles in a 2D grid based on different objects
#' # 4a) based on an object of class "sTopol"
#' dist <- sHexDist(sObj=sTopol) 
#' # 4b) based on an object of class "sMap"
#' dist <- sHexDist(sObj=sI) 

sHexDist <- function(sObj)
{

    if (!is(sObj,"sTopol") & !is(sObj,"sInit") & !is(sObj,"sMap")){
        stop("The funciton must apply to either 'sTopol' or 'sInit' or 'sMap' object.\n")
    }
    
    coord <- sObj$coord
    shape <- sObj$shape
    lattice <- sObj$lattice
    
    if(shape == "sheet" | shape != "sheet"){
        if(lattice == "hexa"){
            dist <- as.matrix(stats::dist(coord))
        }else if(lattice == "rect"){
            dist <- as.matrix(stats::dist(coord, method = "maximum"))
        }
    }
    
    invisible(dist)
}