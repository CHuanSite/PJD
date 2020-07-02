#' Sequential Decomposition with PCA
#'
#' Sequential decomposition of several matrices with PCA
#'
#' @param dataset A list of dataset to be analyzed
#' @param group A list of grouping of the datasets, indicating the relationship between datasets
#' @param comp_num A vector indicates the dimension of each compoent
#'
#' @importFrom RSpectra svds
#'
#' @return A list contains the component and the score of each dataset on every component after seqPCA algorithm
#'
#' @keywords sequential, PCA
#'
#' @examples
#' dataset = list(matrix(runif(5000, 1, 2), nrow = 100, ncol = 50),
#' matrix(runif(5000, 1, 2), nrow = 100, ncol = 50),
#' matrix(runif(5000, 1, 2), nrow = 100, ncol = 50),
#' matrix(runif(5000, 1, 2), nrow = 100, ncol = 50))
#' group = list(c(1,2,3,4), c(1,2), c(3,4), c(1,3), c(2,4), c(1), c(2), c(3), c(4))
#' comp_num = c(2,2,2,2,2,2,2,2,2)
#' res_seqPCA = seqPCA(dataset, group, comp_num)
#'
#' @export

seqPCA <- function(dataset, group, comp_num){

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

    ## compute the total number of component each dataset has
    data_comp_num = rep(0, N)
    for(i in 1 : K){
        for(j in group[[i]]){
            data_comp_num[group[[i]][j]] = data_comp_num[group[[i]][j]] + 1
        }
    }

    ## compute the component for each dataset
    data_comp_total = list()
    for(i in 1 : N){
        data_comp_total[[i]] = svds(dataset[[i]], data_comp_num[i])$u
    }

    ## compute the components sequentially
    temp_comp = c()
    for(i in 1 : K){
        for(j in group[[i]]){
            temp_comp = cbind(temp_comp, data_comp_total[[j]])
        }
        temp_comp_svd = svds(temp_comp, comp_num[i])
        list_component[[i]] = temp_comp_svd$u
        for(i in 1 : N){
            data_comp_total[[i]] = data_comp_total[[i]] - temp_comp_svd$u %*% (t(temp_comp_svd$u) %*% data_comp_total[[i]])
        }
    }

    ## compute the score for each dataset
    for(i in 1 : K){
        for(j in 1 : N){
            if (j %in% group[[i]]){
                list_score[[j]][[i]] = t(list_component[[i]]) %*% dataset[[j]]
            }
        }
    }

    return(list(linked_component_list = list_component, score_list = list_score))
}
