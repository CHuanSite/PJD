#' Joint Decomposition with Principal Component Analysis
#'
#' Joint decomposition of several linked matrices with Principal Component Analysis (PCA)
#'
#' @param dataset A list of dataset to be analyzed
#' @param group A list of grouping of the datasets, indicating the relationship between datasets
#' @param comp_num A vector indicates the dimension of each compoent
#' @param max_ite The maximum number of iterations for the jointPCA algorithms to run, default value is set to 100
#' @param max_err The maximum error of loss between two iterations, or the program will terminate and return, default value is set to be 0.001
#'
#' @importFrom stats runif
#' @return A list contains the component and the score of each dataset on every component after jointPCA algorithm
#'
#' @keywords joint, PCA
#'
#' @export

jointPCA <- function(dataset, group, comp_num, max_ite = 100, max_err = 0.0001){
    ## Parameters to be initialized
    N = length(dataset)
    K = length(group)
    M = sum(comp_num)
    p = nrow(dataset[[1]])

    ## Combine the dataset into a huge one
    combine_data <- c()
    for(i in 1 : N){
        combine_data = cbind(combine_data, dataset[[i]])
    }

    ## List to store the random scores and initialize the scores
    score_list = list()
    for(i in 1 : N){
        score_list[[i]] = list()
    }
    for(i in 1 : N){
        for(j in 1 : K){
            if(i %in% group[[j]]){
                score_list[[i]][[j]] = matrix(runif(comp_num[j] * ncol(dataset[[i]])), nrow = comp_num[j])
            }else{
                score_list[[i]][[j]] = matrix(0, nrow = comp_num[j], ncol = ncol(dataset[[i]]))
            }
        }
    }

    ## Initialize the loss
    loss = 0

    ## Start the Alternative Projection
    for(t in 1 : max_ite){
        matrix_score = c()
        for(i in 1 : N){
            temp_score = c()
            for(j in 1 : K){
                temp_score = rbind(temp_score, score_list[[i]][[j]])
            }
            matrix_score = cbind(matrix_score, temp_score)
        }

        ## Apply Procrustes projection to obtain the linkedin component
        linked_component = Procrustes(combine_data, matrix_score)
        linked_component_list = list()
        index = 1
        for(j in 1 : K){
            linked_component_list[[j]] = linked_component[, index : (index + comp_num[j] - 1)]
            index = index + comp_num[j]
        }

        ## Compute the random scores for every dataset
        for(i in 1 : N){
            for(j in 1 : K){
                if(i %in% group[[j]]){
                    score_list[[i]][[j]] = t(linked_component_list[[j]]) %*% dataset[[i]]
                }else{
                    score_list[[i]][[j]] = matrix(0, nrow = comp_num[j], ncol = ncol(dataset[[i]]))
                }
            }
        }
        loss_current = sum((combine_data - linked_component %*% matrix_score)^2)
        if(abs(loss[length(loss)] - loss_current) / loss[length(loss)] < max_err){
            return(list(linked_component_list = linked_component_list, score_list = score_list))
        }
        loss = c(loss, loss_current)
        print(t)
        print(loss_current)
    }

    return(list(linked_component_list = linked_component_list, score_list = score_list))
}

#' Procrustes Projection
#'
#' Procrustes projection function to solve the procrustes problem ||A - UB||^2 U^T U = I
#'
#' @param A The input matrix A as target
#' @param B the input matrix B as basis
#'
#' @return The procrustes matrix U
#'
#' @keywords procrustes
#'
#' @export

Procrustes <- function(A, B){
    C = A %*% t(B)
    svd.temp = svd(C)
    return(svd.temp$u %*% t(svd.temp$v))
}

