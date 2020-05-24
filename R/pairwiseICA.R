#' Pairwise decomposition with Independent Component Analysis
#'
#' Pairwise decomposition of several linked matrices with Independent Component Analysis (ICA)
#'
#' @param dataset A list of dataset to be analyzed
#' @param group A list of grouping of the datasets, indicating the relationship between datasets
#' @param comp_num A vector indicates the dimension of each compoent
#'
#' @return A list contains the component and the score of each dataset on every component after pairwisePCA algorithm
#'
#' @keywords pairwise, ICA
#'
#' @export

pairwiseICA <- function(dataset, group, comp_num){
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

    for(i in 1 : K){
        list_component[[i]] = matrix(0, nrow = p, ncol = comp_num[i])
        for(j in 1 : N){
            list_score[[j]][[i]] = matrix(0, nrow = comp_num[i], ncol = N_dataset[j])
        }
    }

    ## Extract pairwise PCA from the datasets
    for(i in 1 : K){
        temp_dat = c()
        temp_sample_n = c()
        for(j in group[[i]]){
            temp_dat = cbind(temp_dat, dataset[[j]])
            temp_sample_n = c(temp_sample_n, ncol(datasets[[j]]))
        }
        ica_temp = fastICA(temp_dat, comp_num[i])
        list_component[[i]] = ica_temp$S
        for(j in 1 : length(group[[i]])){
            list_score[[group[[i]][j]]][[i]] = ica_temp$A[, ifelse(j == 1, 1, cumsum(temp_sample_n[j - 1])) : cumsum(temp_sample_n[j])]
        }
    }

    return(list(linked_component_list = list_component, score_list = list_score))
}
