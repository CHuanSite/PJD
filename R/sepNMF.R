#' Single Data Set Decomposition with Nonnegative Matrix Factorization
#'
#' Apply NMF (Nonnegative Matrix Factorization) to a single data set
#'
#' @param dataset A dataframe/matrix to be decomposed
#' @param comp_num Number of NMFs to be extracted
#' @param perturbation A small perturbation to ensure nmf works well
#'
#' @importFrom NMF nmf
#'
#' @return A list of scores and component
#'
#' @keywords separate analysis, NMF
#'
#' @examples
#' dataset = matrix(runif(5000, 1, 2), nrow = 100, ncol = 50)
#' comp_num = 2
#' res_sepNMF = sepNMF(dataset, comp_num)
#'
#' @export

sepNMF <- function(dataset, comp_num, perturbation = 0.0001){
    nmf_temp = nmf(dataset + perturbation, comp_num)
    component = nmf_temp@fit@W
    score = nmf_temp@fit@H

    return(list(component = component, score = score))
}
