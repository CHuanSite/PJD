#' Concatenated Decomposition with Principal Component Analysis
#'
#' Concatenated decomposition of several matrices with Principal Component Analysis (PCA)
#'
#' @param dataset A list of dataset to be analyzed
#' @param group A list of grouping of the datasets, indicating the relationship between datasets
#' @param comp_num A vector indicates the dimension of each compoent
#'
#' @importFrom RSpectra svds
#'
#' @return A list contains the component and the score of each dataset on every component after concatPCA algorithm
#'
#' @keywords pairwise, PCA
#'
#' @examples
#' dataset = list(matrix(runif(5000, 1, 2), nrow = 100, ncol = 50),
#' matrix(runif(5000, 1, 2), nrow = 100, ncol = 50),
#' matrix(runif(5000, 1, 2), nrow = 100, ncol = 50),
#' matrix(runif(5000, 1, 2), nrow = 100, ncol = 50))
#' group = list(c(1,2,3,4), c(1,2), c(3,4), c(1,3), c(2,4), c(1), c(2), c(3), c(4))
#' comp_num = c(2,2,2,2,2,2,2,2,2)
#' res_concatPCA = concatPCA(dataset, group, comp_num)
#'
#' @export

concatPCA <- function(dataset, group, comp_num){

    ## Obtain names for dataset, gene and samples
    dataset_name = datasetNameExtractor(dataset)
    gene_name = geneNameExtractor(dataset)
    sample_name = sampleNameExtractor(dataset)
    group_name = groupNameExtractor(group)

    ## Normalize and Preprocess dataset
    dataset = frameToMatrix(dataset)
    dataset = normalizeData(dataset)
    dataset = balanceData(dataset)

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
            temp_sample_n = c(temp_sample_n, ncol(dataset[[j]]))
        }
        svd_temp = svds(temp_dat, comp_num[i])
        list_component[[i]] = svd_temp$u
        for(j in 1 : length(group[[i]])){
            list_score[[group[[i]][j]]][[i]] = diag(svd_temp$d) %*% t(svd_temp$v)[, ifelse(j == 1, 1, sum(temp_sample_n[1 : (j - 1)]) + 1) : sum(temp_sample_n[1 : j])]
        }
    }

    ## Assign name for components
    list_component = compNameAssign(list_component, group_name)
    list_component = geneNameAssign(list_component, gene_name)
    list_score = scoreNameAssign(list_score, dataset_name, group_name)
    list_score = sampleNameAssign(list_score, sample_name)
    list_score = filterNAValue(list_score, dataset, group)
    list_score = rebalanceData(list_score, group, dataset)

    return(list(linked_component_list = list_component, score_list = list_score))
}
