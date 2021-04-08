#' Data Normalization
#'
#' Normalize data to have mean zero and std 1
#'
#' @param dataset The input list of data sets matrix
#'
#' @keywords normalize
#'
#' @examples
#' dataset = list(
#' matrix(runif(5000, 1, 2), nrow = 100, ncol = 50),
#' matrix(runif(5000, 1, 2), nrow = 100, ncol = 50)
#' )
#' normalizeData(dataset)
#'
#' @export

normalizeData <- function(dataset){
    dataset = lapply(dataset, FUN = function(x){
        scale(t(scale(t(x), scale = FALSE)))
    })
    return(dataset)
}
