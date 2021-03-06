#' Sequential Decomposition with PCA and automatic rank selection
#'
#' Sequential decomposition of several matrices with PCA, the rank selection procedure is automatic based on BEMA
#'
#' @param dataset A list of dataset to be analyzed
#' @param group A list of grouping of the datasets, indicating the relationship between datasets
#' @param threshold The threshold used to cutoff the eigenvalues
#' @param backup A backup variable, which permits the overselection of the components by BEMA
#' @param total_number Total number of components will be extracted, if default value is set to NA, then BEMA will be used.
#'
#' @importFrom RSpectra svds
#'
#' @return A list contains the component and the score of each dataset on every component after seqPCA algorithm
#'
#' @keywords sequential, rank, PCA
#'
#' @examples
#' dataset = list(matrix(runif(5000, 1, 2), nrow = 100, ncol = 50),
#' matrix(runif(5000, 1, 2), nrow = 100, ncol = 50),
#' matrix(runif(5000, 1, 2), nrow = 100, ncol = 50),
#' matrix(runif(5000, 1, 2), nrow = 100, ncol = 50))
#' group = list(c(1, 2, 3, 4), c(1, 2), c(3, 4), c(1, 3), c(2, 4), c(1), c(2), c(3), c(4))
#' threshold = c(3, 1.5, 1.5, 1.5, 1.5, 0.5, 0.5, 0.5, 0.5)
#' res_seqPCA = seqPCA.rank(dataset, group, threshold = threshold)
#'
#' @export

seqPCA.rank <- function(dataset, group, total_number = NULL, threshold, backup = 0){

    ## Parameters to be initialized
    N = length(dataset)
    K = length(group)
    p = nrow(dataset[[1]])
    N_dataset = unlist(lapply(dataset, ncol))

    ## Determine the number of components

    ## compute the component for each dataset
    data_comp_total = list()

    ## compute the total number of component each dataset has
    data_comp_num = rep(0, N)

    ## Using BEMA to extract the number of components
    for(i in 1 : N){
        svd_temp = svd(dataset[[i]])
        if(is.null(total_number)){
            data_comp_num[i] = BEMA(svd_temp$d^2 / ncol(dataset[[i]]), p = nrow(dataset[[i]]), n = ncol(dataset[[i]])) + backup
        }else{
            data_comp_num[i] = total_number[i] + backup
        }
        data_comp_total[[i]] = svd_temp$u[, 1 : data_comp_num[i]]
    }

    ## Output the component and scores
    list_component = list()
    list_score = list()
    for(j in 1 : N){
        list_score[[j]] = list()
    }

    for(i in 1 : K){
        list_component[[i]] = matrix(0, nrow = p, ncol = 1)
        for(j in 1 : N){
            list_score[[j]][[i]] = matrix(0, nrow = 1, ncol = N_dataset[j])
        }
    }

    ## compute the components sequentially
    for(i in 1 : K){
        temp_comp = c()
        for(j in group[[i]]){
            temp_comp = cbind(temp_comp, data_comp_total[[j]])
        }

        ## Orthogonalize the component
        if(i >= 2){
            for(j in 1 : (i - 1)){
                temp_comp = temp_comp - list_component[[j]] %*% (t(list_component[[j]]) %*% temp_comp)
            }
        }

        temp_comp_svd = svd(temp_comp)
        plot(temp_comp_svd$d^2, xlab = "index of eigenvalue", ylab = "eigenvalue", main = paste0("Group of Dataset: ", toString(group[[i]])), type = "o")

        index = which(temp_comp_svd$d^2 > threshold[i])
        if(length(index) > 0){
            list_component[[i]] = temp_comp_svd$u[, 1 : length(index)]
            for(j in 1 : N){
                data_comp_total[[j]] = data_comp_total[[j]] - list_component[[i]] %*% (t(list_component[[i]]) %*% data_comp_total[[j]])
                # data_comp_total[[j]] = svd(data_comp_total[[j]])$u[, 1 : (ncol(data_comp_total[[j]]) - comp_num[i])]
            }
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
