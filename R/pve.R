#' Percentage of Variance Explained for separate data set
#'
#' Compute the PVE (percentage of variance explained) for each data set
#'
#' @param dataset A list of data sets for input
#' @param list_score A list of extracted scores by the corresponding algorithm
#' @param list_component A list of components comptuted by the corresponding algorithm
#'
#' @return The list of scores
#'
#' @keywords PVE, Separate analysis
#'
#' dataset = list(matrix(runif(5000, 1, 2), nrow = 100, ncol = 50))
#' comp_num = 2
#' res_sepPCA = sepPCA(dataset, comp_num)
#' pveSep(dataset, res_sepPCA$score_list, res_sepPCA$linked_component_list)
#'
#' @export

pveSep <- function(dataset, list_score, list_component){
    ## Compute total variance for each data set
    total_variance = c()
    for(i in 1 : length(dataset)){
        total_variance = c(total_variance, sum(dataset[[i]] * dataset[[i]]))
    }

    ## Computer percentage of variance explained by each component for every data set
    pve_store = c()
    for(i in 1 : length(list_component)){
        u = svd(list_component[[i]])$u
        pve_store = c(pve_store, sum(diag((t(u) %*% dataset[[i]]) %*% t(t(u) %*% dataset[[i]]))) / total_variance[i])
    }

    ## Assign new name for data set
    for(i in 1 : length(list_score)){
        names(list_score)[i] = paste0(names(list_score)[i], ", PVE: ", formatC(pve_store[i], digits = 6, format = "f"))
    }

    return(list_score)
}




#' Percentage of Variance Explained for Multiple Data sets
#'
#' Compute the PVE (Percentage of Variance Explained) for multiple data sets on multiple components
#'
#' @param dataset A list of dataset to be analyzed
#' @param group A list of grouping of the datasets, indicating the relationship between datasets
#' @param comp_num A vector indicates the dimension of each compoent
#' @param list_score A list of extracted scores by the corresponding algorithm
#' @param list_component A list of components comptuted by the corresponding algorithm
#'
#' @return The list of scores
#'
#' @keywords PVE, Multiple analysis
#'
#' dataset = list(matrix(runif(5000, 1, 2), nrow = 100, ncol = 50),
#' matrix(runif(5000, 1, 2), nrow = 100, ncol = 50),
#' matrix(runif(5000, 1, 2), nrow = 100, ncol = 50),
#' matrix(runif(5000, 1, 2), nrow = 100, ncol = 50))
#' group = list(c(1,2,3,4), c(1,2), c(3,4), c(1,3), c(2,4), c(1), c(2), c(3), c(4))
#' comp_num = c(2,2,2,2,2,2,2,2,2)
#' res_concatPCA = concatPCA(dataset, group, comp_num)
#' pveMultiple(dataset, group, comp_num, res_concatPCA$score_list, res_concatPCA$linked_component_list)
#'
#'@export

pveMultiple <- function(dataset, group, comp_num, list_score, list_component){
    ## Compute total variance for each data set
    total_variance = c()
    for(i in 1 : length(dataset)){
        total_variance = c(total_variance, sum(dataset[[i]] * dataset[[i]]))
    }

    ## Computer percentage of variance explained by each component for every score list
    pve_store = list()
    for(i in 1 : length(group)){
        u = svd(list_component[[i]])$u
        pve_store[[i]] = rep(0, length(dataset))
        for(j in group[[i]]){
            pve_store[[i]][j] =  sum(diag((t(u) %*% dataset[[j]]) %*% t(t(u) %*% dataset[[j]]))) / total_variance[j]
        }
    }

    ## Assign new name for each element in score list
    for(i in 1 : length(group)){
        for(j in group[[i]]){
            names(list_score[[j]])[i] = paste0(names(list_score[[j]])[i], ", PVE: ", formatC(pve_store[[i]][j], digits = 6, format = "f"))
        }
    }

    return(list_score)
}
