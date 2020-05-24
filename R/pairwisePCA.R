#' Pairwise Decomposition with Principal Component Analysis
#'
#' Pairwise decomposition of several matrices with Principal Component Analysis (PCA)
#'
#' @keywords pairwise, PCA
#'
#' @export

pairwisePCA <- function(dataset, group, comp_num, max_ite = 100, max_err = 0.0001){
    ## Parameters to be initialized
    N = length(dataset)
    K = length(group)
    M = sum(comp_num)
    p = nrow(dataset[[1]])

    N_dataset = unlist(lapply(dataset, ncol))


    ## Output the component and scores
    list_component = list()
    list_score = list()
    for(j in 1 : N){
        list_score[[j]] = list()
    }

    for(i in 1 : M){
        list_component[[i]] = matrix(0, nrow = p, ncol = comp_num[i])
        for(j in 1 : N){
            list_score[[j]][[i]] = matrix(0, nrow = comp_num[i], ncol = N_dataset[j])
        }
    }




}
