#' Two-staged Independent Linked Component Analysis
#'
#' Two-staged Independent Linked Component Analysis, a generalization based on the Two-staged Independent Linked Component Analysis
#'
#' @param dataset A list of dataset to be analyzed
#' @param group A list of grouping of the datasets, indicating the relationship between datasets
#' @param comp_num A vector indicates the dimension of each compoent
#' @param backup A positive scalar to determine how many ICs to over select
#'
#' @importFrom RSpectra svds
#' @importFrom fastICA fastICA
#'
#' @return A list contains the component and the score of each dataset on every component after 2siLCA algorithm
#'
#' @keywords two-staged, independent LCA
#'
#' @examples
#' dataset = list(matrix(runif(5000, 1, 2), nrow = 100, ncol = 50),
#' matrix(runif(5000, 1, 2), nrow = 100, ncol = 50),
#' matrix(runif(5000, 1, 2), nrow = 100, ncol = 50),
#' matrix(runif(5000, 1, 2), nrow = 100, ncol = 50))
#' group = list(c(1, 2, 3, 4), c(1, 2), c(3, 4), c(1, 3), c(2, 4), c(1), c(2), c(3), c(4))
#' comp_num = c(2, 2, 2, 2, 2, 2, 2, 2, 2)
#' res_twoStageiLCA = twoStageiLCA(dataset, group, comp_num)
#'
#' @export


twoStageiLCA <- function(dataset, group, comp_num, backup = 0){

    ## Obtain names for dataset, gene and samples
    dataset_name = datasetNameExtractor(dataset)
    gene_name = geneNameExtractor(dataset)
    sample_name = sampleNameExtractor(dataset)
    group_name = groupNameExtractor(group)

    dataset = frameToMatrix(dataset)
    dataset = normalizeData(dataset)

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
            data_comp_num[j] = data_comp_num[j] + comp_num[i] + backup
        }
    }

    ## compute the component for each dataset
    data_comp_total = list()
    for(i in 1 : N){
        data_comp_total[[i]] = svds(dataset[[i]], data_comp_num[i])$u
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

        temp_comp_svd = svds(temp_comp, comp_num[i])
        # print(svd(temp_comp)$d^2)
        plot(svd(temp_comp)$d^2, xlab = "index of eigenvalue", ylab = "eigenvalue", main = paste0("Group of Dataset: ", toString(group[[i]])), type = "o")
        list_component[[i]] = temp_comp_svd$u

        # ## Orthogonalize the component
        # if(i >= 2){
        #     for(j in 1 : (i - 1)){
        #         list_component[[i]] = list_component[[i]] - list_component[[j]] %*% (t(list_component[[j]]) %*% list_component[[i]])
        #     }
        # }

        ## Project the extracted space onto orthogonal space
        # for(j in group[[i]]){
        for(j in 1 : N){
            data_comp_total[[j]] = data_comp_total[[j]] - list_component[[i]] %*% (t(list_component[[i]]) %*% data_comp_total[[j]])
            # data_comp_total[[j]] = svd(data_comp_total[[j]])$u[, 1 : (ncol(data_comp_total[[j]]) - comp_num[i])]
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

    ## Conduct ICA on the extracted Scores
    for(i in 1 : K){
        for(j in 1 : N){
            if(j %in% group[[i]] & nrow(list_score[[j]][[i]]) >= 2){
                temp_ica_res = fastICA(t(list_score[[j]][[i]]), n.comp = nrow(list_score[[j]][[i]]))
                list_score[[j]][[i]] = t(temp_ica_res$S)
            }
        }
    }

    ## Assign name for components
    list_component = compNameAssign(list_component, group_name)
    list_component = geneNameAssign(list_component, gene_name)
    list_score = scoreNameAssign(list_score, dataset_name, group_name)
    list_score = sampleNameAssign(list_score, sample_name)


    return(list(linked_component_list = list_component, score_list = list_score))
}
