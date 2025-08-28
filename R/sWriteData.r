#' Function to write out the best-matching hexagons and/or cluster bases in terms of data
#'
#' \code{sWriteData} is supposed to write out the best-matching hexagons and/or cluster bases in terms of data. 
#'
#' @param sMap an object of class "sMap" or a codebook matrix
#' @param data a data frame or matrix of input data
#' @param sBase an object of class "sBase"
#' @param filename a character string naming a filename
#' @param keep.data logical to indicate whether or not to also write out the input data. By default, it sets to false for not keeping it. It is highly expensive to keep the large data sets
#' @return 
#' a data frame with following components:
#' \itemize{
#'  \item{\code{ID}: ID for data. It inherits the rownames of data (if exists). Otherwise, it is sequential integer values starting with 1 and ending with dlen, the total number of rows of the input data}
#'  \item{\code{Hexagon_index}: the index for best-matching hexagons}
#'  \item{\code{Qerr_distance}: the quantification error (distance) for best-matching hexagons}
#'  \item{\code{Cluster_base}: optional, it is only appended when sBase is given. It stores the cluster memberships/bases}
#'  \item{\code{data}: optional, it is only appended when keep.data is true}
#' }
#' @note If "filename" is not NULL, a tab-delimited text file will be also written out. If "sBase" is not NULL and comes from the "sMap" partition, then cluster bases are also appended. if "keep.data" is true, the data will be part of output.
#' @export
#' @seealso \code{\link{sBMH}}
#' @include sWriteData.r
#' @examples
#' # 1) generate an iid normal random matrix of 100x10 
#' data <- matrix( rnorm(100*10,mean=0,sd=1), nrow=100, ncol=10) 
#'
#' # 2) get trained using by default setup 
#' sMap <- sPipeline(data=data)
#'
#' # 3) write data's BMH hitting the trained map
#' output <- sWriteData(sMap=sMap, data=data, filename="sData_output.txt") 
#'
#' # 4) partition the grid map into cluster bases
#' sBase <- sDmatCluster(sMap=sMap, which_neigh=1,
#' distMeasure="median", clusterLinkage="average") 
#'
#' # 5) write data's BMH and cluster bases
#' output <- sWriteData(sMap=sMap, data=data, sBase=sBase, filename="sData_base_output.txt")

sWriteData <- function(sMap, data, sBase=NULL, filename=NULL, keep.data=FALSE)
{
    
    response <- sBMH(sMap=sMap, data=data, which_bmh="best")    
    bmh <- as.vector(response$bmh)
    
    if(is.null(rownames(data))){
        rownames(data) <- seq(1,nrow(data))
    }
    if(is.null(colnames(data))){
        colnames(data) <- seq(1,ncol(data))
    }
        
    ## Prepare output data
    ## 1st column for "ID"
    ## 2nd column for "Hexagon_index"
    ## 3rd column for "Qerr_distance"
    data_output <- tibble::tibble(ID=rownames(data), Hexagon_index=bmh, Qerr_distance=as.vector(response$qerr))
    
    ## The column for "Cluster_base" (if sBase is given)
    if(!is.null(sBase)){
        if(is(sBase,"sBase")){
            if(sMap$nHex == length(sBase$bases)){
            	data_output <- data_output %>% dplyr::mutate(Cluster_base=sBase$bases[bmh])
            }
        }
    }
    
    ## The next columns for data itself (if keep.data is true)
    if(keep.data){
        output <- dplyr::bind_cols(data_output, tibble::as_tibble(data))
    }else{
    	output <- data_output
    }
    
    ## convert into a data frame called 'output'
    #output <- as.data.frame(data_output, stringsAsFactors=FALSE)
    
    ## If the filename is given, output data is written into a tab-delimited text file
    if(!is.null(filename)){
    	readr::write_delim(output, filename, delim="'t")
    }

    invisible(output)
    
}