#' Joint Decomposition with Nonnegative Matrix Factorization
#'
#' Joint decomposition of several linked matrices with Nonnegative Matrix Factorization (NMF)
#' It is based on the MSE loss, proposed by Lee, Daniel D., and H. Sebastian Seung. "Learning the parts of objects by non-negative matrix factorization." Nature 401.6755 (1999): 788-791.
#'
#' @param dataset A list of dataset to be analyzed
#' @param group A list of grouping of the datasets, indicating the relationship between datasets
#' @param comp_num A vector indicates the dimension of each compoent
#' @param max_ite The maximum number of iterations for the jointNMF algorithms to run, default value is set to 100
#' @param max_err The maximum error of loss between two iterations, or the program will terminate and return, default value is set to be 0.001
#'
#' @return A list contains the component and the score of each dataset on every component after jointNMF algorithm
#'
#' @keywords joint, NMF
#'
#' @export

jointNMF <- function(dataset, group, comp_num, max_ite = 100, max_err = 0.0001){
    ## Initialize values for the algorithm
    N = length(dataset)
    K = length(group)
    M = sum(comp_num)
    p = nrow(dataset[[1]])

    N_dataset = unlist(lapply(dataset, ncol))
    max_element = -Inf
    min_element = Inf
    for(i in 1 : N){
        max_element = max(max_element, max(dataset[[i]]))
        min_element = min(min_element, min(dataset[[i]]))
    }

    ## Initialize the values of W and H
    X = c()
    for(i in 1 : N){
        X = cbind(X, dataset[[i]])
    }
    W = matrix(runif(p * M, min_element, max_element), nrow = p)
    H = c()

    for(i in 1 : K){
        H_temp = c()
        for(j in 1 : N){
            if(j %in% group[[i]]){
                H_temp = cbind(H_temp, matrix(runif(N_dataset[j] * comp_num[i], min_element, max_element), nrow = comp_num[i], ncol = N_dataset[j]))
            }else{
                H_temp = cbind(H_temp, matrix(0, nrow = comp_num[i], ncol = N_dataset[j]))
            }
        }
        H = rbind(H, H_temp)
    }

    ## Iteratively estimate the NMF with Euclidean distance
    error_out = c()

    for(ite in 1 : max_ite){
        H = H * t(W) %*% X / (t(W) %*% W %*% H)
        W = W * (X %*% t(H)) / (W %*% H %*% t(H))
        error_out = c(error_out, sum((X - W %*% H)^2))
        ## Break when the error difference is small
        if(abs(error_out[length(error_out)] - error_out[length(error_out) - 1]) <= max_err){
            break
        }
        print(ite)
        print(abs(error_out[length(error_out)] - error_out[length(error_out) - 1]))
    }

    ## Output component and scores
    list_component = list()
    list_score = list()
    for(j in 1 : N){
        list_score[[j]] = list()
    }

    for(i in 1 : M){
        list_component[[i]] = W[, ifelse(i == 1, 1, cumsum(comp_num)[i - 1] + 1) : cumsum(comp_num)[i]]
        for(j in 1 : N){
            list_score[[j]][[i]] = H[ifelse(i == 1, 1, cumsum(comp_num)[i - 1] + 1) : cumsum(comp_num)[i], ifelse(j == 1, 1, cumsum(N_dataset)[j - 1] + 1) : cumsum(N_dataset)[j]]
        }
    }

    return(list(linked_component_list = list_component, score_list = list_score))
}
