#' Single Data Set Decomposition with Independent Component Analysis
#'
#' Apply ICA (Independent Component Analysis) to a single data set
#'
#' @param dataset A dataframe/matrix to be decomposed
#' @param comp_num Number of ICs to be extracted
#'
#'
#' @importFrom fastICA fastICA
#'
#' @return A list of scores and component
#'
#' @keywords separate analysis, ICA
#'
#' @examples
#' dataset = matrix(runif(5000, 1, 2), nrow = 100, ncol = 50)
#' comp_num = 2
#' res_sepICA = sepICA(dataset, comp_num)
#'
#' @export

sepICA <- function(dataset, comp_num){
    ica_temp = fastICA(dataset, comp_num)
    component = ica_temp$S
    score = ica_temp$A

    return(list(component = component, score = score))
}
