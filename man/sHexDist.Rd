% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sHexDist.r
\name{sHexDist}
\alias{sHexDist}
\title{Function to calculate distances between hexagons/rectangles in a 2D grid}
\usage{
sHexDist(sObj)
}
\arguments{
\item{sObj}{an object of class "sTopol" or "sInit" or "sMap"}
}
\value{
\itemize{
\item{\code{dist}: a symmetric matrix of nHex x nHex, containing
pairwise distances, where nHex is the total number of
hexagons/rectanges in the grid}
}
}
\description{
\code{sHexDist} is supposed to calculate euclidian distances between
each pair of hexagons/rectangles in a 2D grid of input "sTopol" or
"sMap" object. It returns a symmetric matrix containing pairwise
distances.
}
\note{
The return matrix has rows/columns ordered in the same order as the
"coord" matrix of the input object does.
}
\examples{
# 1) generate an iid normal random matrix of 100x10 
data <- matrix( rnorm(100*10,mean=0,sd=1), nrow=100, ncol=10)

# 2) from this input matrix, determine nHex=5*sqrt(nrow(data))=50, 
# but it returns nHex=61, via "sHexGrid(nHex=50)", to make sure a supra-hexagonal grid
sTopol <- sTopology(data=data, lattice="hexa", shape="suprahex")

# 3) initialise the codebook matrix using "uniform" method
sI <- sInitial(data=data, sTopol=sTopol, init="uniform")

# 4) calculate distances between hexagons/rectangles in a 2D grid based on different objects
# 4a) based on an object of class "sTopol"
dist <- sHexDist(sObj=sTopol)
# 4b) based on an object of class "sMap"
dist <- sHexDist(sObj=sI)
}
\seealso{
\code{\link{sTopology}}, \code{\link{sInitial}}
}
