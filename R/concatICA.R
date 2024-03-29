#' Concatenated decomposition with Independent Component Analysis
#'
#' Concatenated decomposition of several linked matrices with Independent Component Analysis (ICA)
#'
#' @param dataset A list of dataset to be analyzed
#' @param group A list of grouping of the datasets, indicating the relationship between datasets
#' @param comp_num A vector indicates the dimension of each compoent
#' @param weighting Weighting of each dataset, initialized to be NULL
#' @param proj_dataset The datasets to be projected on
#' @param proj_group The grouping of projected data sets
#'
#' @importFrom fastICA fastICA
#'
#' @return A list contains the component and the score of each dataset on every component after concatPCA algorithm
#'
#' @keywords pairwise, ICA
#'
#' @examples
#' dataset = list(matrix(runif(5000, 1, 2), nrow = 100, ncol = 50),
#' matrix(runif(5000, 1, 2), nrow = 100, ncol = 50),
#' matrix(runif(5000, 1, 2), nrow = 100, ncol = 50),
#' matrix(runif(5000, 1, 2), nrow = 100, ncol = 50))
#' group = list(c(1,2,3,4), c(1,2), c(3,4), c(1,3), c(2,4), c(1), c(2), c(3), c(4))
#' comp_num = c(2,2,2,2,2,2,2,2,2)
#' proj_dataset = list(matrix(runif(5000, 1, 2), nrow = 100, ncol = 50))
#' proj_group = list(c(TRUE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE))
#' res_concatICA = concatICA(dataset, group, comp_num, proj_dataset = proj_dataset, proj_group = proj_group)
#'
#' @export

concatICA <- function(dataset, group, comp_num, weighting = NULL, proj_dataset = NULL, proj_group = NULL){
    concatPCA_out = concatPCA(dataset, group, comp_num, weighting)

    ## Obtain names for dataset, gene and samples
    dataset_name = datasetNameExtractor(dataset)
    gene_name = geneNameExtractor(dataset)
    sample_name = sampleNameExtractor(dataset)
    group_name = groupNameExtractor(group)

    dataset = frameToMatrix(dataset)
    dataset = normalizeData(dataset)
    dataset = balanceData(dataset)
    dataset = weightData(dataset, weighting)

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

    ## Extract pairwise ICA from the datasets
    for(i in 1 : K){
        temp_dat = c()
        temp_sample_n = c()
        for(j in group[[i]]){
            temp_dat = cbind(temp_dat, concatPCA_out$score_list[[j]][[i]] / sqrt(ncol(dataset[[j]])))
            temp_sample_n = c(temp_sample_n, ncol(dataset[[j]]))
        }
        ica_temp = fastICA(t(temp_dat), comp_num[i])
        list_component[[i]] = concatPCA_out$linked_component_list[[i]] %*% t(ica_temp$A)
        for(j in 1 : length(group[[i]])){
            list_score[[group[[i]][j]]][[i]] = t(ica_temp$S)[, ifelse(j == 1, 1, sum(temp_sample_n[1 : (j - 1)]) + 1) : sum(temp_sample_n[1 : j])]
        }
    }

    ## Assign name for components
    list_component = compNameAssign(list_component, group_name)
    list_component = geneNameAssign(list_component, gene_name)
    list_score = scoreNameAssign(list_score, dataset_name, group_name)
    list_score = sampleNameAssign(list_score, sample_name)
    list_score = filterNAValue(list_score, dataset, group)
    list_score = rebalanceData(list_score, group, dataset)
    list_score = pveMultiple(dataset, group, comp_num, list_score, list_component)

    ## Project score
    proj_list_score = list()
    if(!is.null(proj_dataset)){
        proj_sample_name = sampleNameExtractor(proj_dataset)
        proj_dataset_name = datasetNameExtractor(proj_dataset)

        proj_dataset = frameToMatrix(proj_dataset)
        proj_dataset = normalizeData(proj_dataset)
        proj_dataset = balanceData(proj_dataset)

        for(i in 1 : length(proj_dataset)){
            proj_list_score[[i]] = list()
            for(j in 1 : length(proj_group[[i]])){
                if(proj_group[[i]][j]){
                    proj_list_score[[i]][[j]] =  t(list_component[[j]]) %*% proj_dataset[[i]]
                }else{
                    proj_list_score[[i]][[j]] = NA
                }
            }
        }


        proj_list_score = scoreNameAssign(proj_list_score, proj_dataset_name, group_name)
        proj_list_score = sampleNameAssign(proj_list_score, proj_sample_name)
    }

    return(list(linked_component_list = list_component, score_list = list_score, proj_score_list = proj_list_score))
}
