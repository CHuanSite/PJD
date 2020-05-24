#' Joint decomposition with Independent Component Analysis
#'
#' Joint decomposition of several linked matrices with Independent Component Analysis (ICA)
#'
#' @param dataset A list of dataset to be analyzed
#' @param group A list of grouping of the datasets, indicating the relationship between datasets
#' @param comp_num A vector indicates the dimension of each compoent
#' @param max_ite The maximum number of iterations for the jointPCA algorithms to run, default value is set to 100
#' @param max_err The maximum error of loss between two iterations, or the program will terminate and return, default value is set to be 0.001
#'
#' @importFrom fastICA fastICA
#'
#' @return A list contains the component and the score of each dataset on every component after jointPCA algorithm
#'
#' @keywords Joint, ICA
#'
#' @export
#'

jointICA <- function(dataset, group, comp_num, max_ite = 100, max_err = 0.0001){
    ## Parameters to be initialized
    N = length(dataset)
    K = length(group)
    M = sum(comp_num)
    p = nrow(dataset[[1]])
    N_dataset = unlist(lapply(dataset, ncol))

    ## PCA preprocessing on the data
    PCA_preprocess = jointPCA(dataset, group, comp_num, max_ite, max_err)

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

    ## Applying joint ICA
    for(i in 1 : K){
        dat_temp = c()
        n_sample_temp = c()
        for(j in group[[i]]){
            dat_temp = cbind(dat_temp, PCA_preprocess$linked_component_list[[i]] %*% list_score[[j]][[i]])
            n_sample_temp = c(n_sample_temp, N_dataset[j])
        }
        ica_temp = fastICA(dat_temp, comp_num[i])
        list_component[[i]] = ica_temp$S
        for(j in 1 : length(group[[i]])){
            list_score[[group[[i]][j]]][[i]] = ica_temp$A[, ifelse(j == 1, 1, cumsum(n_sample_temp[j - 1]) + 1) : cumsum(n_sample_temp[j])]
        }
    }

    return(list(linked_component_list = list_component, score_list = list_score))
}
