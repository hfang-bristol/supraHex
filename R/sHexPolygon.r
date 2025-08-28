#' Function to extract polygon location per hexagon within a supra-hexagonal grid
#'
#' \code{sHexPolygon} is supposed to extract polygon location per hexagon within a supra-hexagonal grid
#'
#' @param sObj an object of class "sMap" or "sInit" or "sTopol" or "sHex"
#' @param area.size an integer or a vector specifying the area size of each hexagon
#' @return 
#' a tibble of 7 columns ('index','x','y','node','edge','stepCentroid','angleCentroid') storing polygon location per hexagon. 'node' for nodes (including n1,n2,n3,n4,n5,n6), and 'edge' for a list-column where each is a tibble with a single column 'edge' containing two rows (such as edges 'e12' and 'e16' for the node 'n1').
#' @note None
#' @export
#' @seealso \code{\link{sHexGridVariant}}, \code{\link{sPipeline}}
#' @include sHexPolygon.r
#' @examples
#' sObj <- sTopology(xdim=4, ydim=4, lattice="hexa", shape="suprahex")
#' df_polygon <- sHexPolygon(sObj, area.size=1)

sHexPolygon <- function (sObj, area.size=1)
{
    
    if (!is(sObj,"sHex") & !is(sObj,"sTopol") & !is(sObj,"sInit") & !is(sObj,"sMap")){
        stop("The funciton must apply to either 'sHex' or 'sTopol' or 'sInit' or 'sMap' object.\n")
    }
    
    dat <- data.frame(sObj$coord)
    #xdim <- sObj$xdim
    #ydim <- sObj$ydim
    nHex <- sObj$nHex
    xnew <- dat$x
    ynew <- dat$y
    
    ####################################
    stepCentroid <- ''
    angleCentroid <- ''
    sHex <- ''
    
    if(is(sObj, "sHex")){
    	sHex <- sObj
    }else{
    	shape <- sObj$shape
    	if(shape != 'sheet'){
    		sHex <- sHexGridVariant(nHex=nHex, shape=shape)
    	}
    }
    
    if(is(sHex, 'sHex')){
    	stepCentroid <- sHex$stepCentroid
    	angleCentroid <- sHex$angleCentroid
    }
    plts <- rep(stepCentroid, each=6)
    plta <- rep(angleCentroid, each=6)
    ####################################
    
	###################################
    hexC <- list()
    hexC$x <- c(0.5,0.5,0,-0.5,-0.5,0)
    hexC$y <- c(-1*sqrt(3)/6, sqrt(3)/6, sqrt(3)/3, sqrt(3)/6, -1*sqrt(3)/6, -1*sqrt(3)/3)
	###################################
    
    minSize <- 0.25
    if(length(area.size) == nHex){
        scaled <- (area.size - min(area.size))/(max(area.size)-min(area.size))
        scaled[scaled <= minSize] <- minSize ## replace those <= minSize with minSize
        area.size <- scaled
    }else if(length(area.size) == 1){
        if(area.size > 1 | area.size <= minSize){
            area.size <- 1
        }
    }else{
        area.size <- 1
    }
	###################################
	
    radius <- rep.int(1,nHex) * area.size
    n6 <- rep.int(6:6, nHex)
    
    pltx <- rep.int(hexC$x, nHex) * rep.int(radius, n6) + rep.int(xnew, n6)
    plty <- rep.int(hexC$y, nHex) * rep.int(radius, n6) + rep.int(ynew, n6)
    
	## return location per hexagon
    index <- rep(1:nHex, each=6)
    #df_polygon <- data.frame(x=pltx, y=plty, index=index, stepCentroid=plts, angleCentroid=plta, stringsAsFactors=FALSE)
    df_polygon <- data.frame(x=round(pltx,digits=2), y=round(plty,digits=3), index=index, stepCentroid=plts, angleCentroid=plta, stringsAsFactors=FALSE)

    ## append "node" and "edge" (list-column)
    edge <- node <- group <- x <- y <- index <- tmp <- data <- NULL
    tibble::tibble(node=1:6, edge=stringr::str_c(c(2:6,1),',',c(6,1:5))) %>% tidyr::separate_rows(edge) %>% dplyr::mutate(group=ifelse(node<edge,stringr::str_c(node,edge),stringr::str_c(edge,node))) %>% dplyr::transmute(node=stringr::str_c('n',node), edge=stringr::str_c('e',group)) %>% tidyr::nest(edge=-node) -> df_node_edge
    df_polygon %>% tidyr::nest(data=-index) %>% dplyr::mutate(tmp=purrr::map(data, ~.x %>% dplyr::mutate(node=df_node_edge %>% dplyr::pull(node), edge=df_node_edge %>% dplyr::pull(edge)))) %>% tidyr::unnest(tmp) %>% dplyr::select(index,x,y,node,edge,stepCentroid,angleCentroid) -> df_polygon

    invisible(df_polygon)
}