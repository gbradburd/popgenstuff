#' Function for making standardized color range
#' 
#' @param x A vector of numeric values for which 
#'			you want to generate colors that indicate 
#'			value.
#' @param cols A palette or vector of colors that you want 
#'			to plot in
#' @param nCols The total number of colors you want in the absolute 
#'			(standardized) scale
#' @param valRange The range of values for which you want to 
#'			generate colors. If unspecified, the function will 
#'			use the range of \code{x}. Default is \code{NULL}.
#' @return A vector of colors that relate the value of each 
#'			element of \code{x} within the range specified by 
#'			\code{valRange}.
#' @export
colFunc <- function (x, cols, nCols, valRange) {
	if(is.null(valRange)){
		valRange <- c(min(x),max(x))
	}
    cols <- grDevices::colorRampPalette(cols)(nCols)[findInterval(x, seq(valRange[1], valRange[2], length.out = nCols))]
    return(cols)
}