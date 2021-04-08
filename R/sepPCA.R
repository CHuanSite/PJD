#' Single Data Set Decomposition with Principal Component Analysis
#'
#' Apply PCA (Principal Component Analysis) to a single data set
#'
#' @param dataset A dataframe/matrix to be decomposed
#' @param comp_num Number of PCs to be extracted
#'
#'
#' @importFrom RSpectra svds
#'
#' @return A list of scores and component
#'
#' @keywords separate analysis, PCA
#'
#' @examples
#' dataset = list(matrix(runif(5000, 1, 2), nrow = 100, ncol = 50))
#' comp_num = 2
#' res_sepPCA = sepPCA(dataset, comp_num)
#'
#' @export

sepPCA <- function(dataset, comp_num){
    dataset = frameToMatrix(dataset)

    N = length(dataset)

    list_component = list()
    list_score = list()

    for(i in 1 : N){
        svd_temp = svds(dataset[[i]], comp_num[i])
        component = svd_temp$u
        score = diag(svd_temp$d) %*% t(svd_temp$v)
        list_component[[i]] = component
        list_score[[i]] = score
    }

    return(list(linked_component_list = list_component, score_list = list_score))
}
