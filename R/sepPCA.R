#' Single Data Set Decomposition with Principal Component Analysis
#'
#' Apply PCA (Principal Component Analysis) to a single data set
#'
#' @param dataset A dataframe/matrix to be decomposed
#' @param comp_num Number of PCs to be extracted
#' @param weighting Weighting of each dataset, initialized to be NULL
#'
#' @importFrom RSpectra svds
#'
#' @return A list of scores and component
#'
#' @keywords separate analysis, PCA
#'
#' @examples
#' dataset = list(matrix(runif(5000, 1, 2), nrow = 100, ncol = 50))
#' comp_num = 2
#' res_sepPCA = sepPCA(dataset, comp_num)
#'
#' @export

sepPCA <- function(dataset, comp_num, weighting = NULL){

    ## Obtain names for dataset, gene and samples
    dataset_name = datasetNameExtractor(dataset)
    gene_name = geneNameExtractor(dataset)
    sample_name = sampleNameExtractor(dataset)

    ## Prepare dataset
    dataset = frameToMatrix(dataset)
    dataset = normalizeData(dataset)
    dataset = weightData(dataset, weighting)

    N = length(dataset)

    list_component = list()
    list_score = list()

    for(i in 1 : N){
        svd_temp = svds(dataset[[i]], comp_num[i])
        component = svd_temp$u
        score = diag(svd_temp$d) %*% t(svd_temp$v)
        list_component[[i]] = component
        list_score[[i]] = score
    }

    ## Assign name for components
    list_component = compNameAssignSep(list_component, dataset_name)
    list_component = geneNameAssign(list_component, gene_name)
    list_score = scoreNameAssignSep(list_score, dataset_name)
    list_score = sampleNameAssignSep(list_score, sample_name)
    list_score = pveSep(dataset, list_score, list_component)

    return(list(linked_component_list = list_component, score_list = list_score))
}
