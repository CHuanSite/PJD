for(i in 1 : 1){
    set.seed(i)

    ################################################################
    ##
    ## Simulate data
    ##
    ################################################################

    ## Simulation part of the algorithm
    configuration_setting = configuration_setting_generation(featureNum = 1000,
                                                             DataNum = c(30, 200, 30, 200),
                                                             noiseVariance = c(1, 1, 1, 1))
    data_list = simulated_data_generation(configuration_setting, 10)
    cov_list = list(cov(t(data_list[[1]])), cov(t(data_list[[2]])), cov(t(data_list[[3]])), cov(t(data_list[[4]])))
    eigen_space = eigen(cov(t(data_list[[1]])) + cov(t(data_list[[2]])) + cov(t(data_list[[3]])) + cov(t(data_list[[4]])))$vectors[, 1 : 18]
    group = list(c(1, 2, 3, 4), c(1, 2), c(3, 4), c(1, 3), c(2, 4), c(1), c(2), c(3), c(4))
    factor_num = c(2, 2, 2, 2, 2, 2, 2, 2, 2)

    ################################################################
    ##
    ## Compare the result of the three LCA algorithm
    ##
    ################################################################
    concatPCA.result <- concatPCA(data_list, group, factor_num)
    jointPCA.result <- jointPCA(data_list, group, factor_num)
    linkedPCA.result <- linkedPCA(data_list, cov_list, eigen_space, group, factor_num)

    ############################################################
    ##
    ##
    ## LCA.1 result
    ##
    ##
    ############################################################
    comp.1 = cbind(concatPCA.result$linked_component_list[[1]], concatPCA.result$linked_component_list[[2]], concatPCA.result$linked_component_list[[4]], concatPCA.result$linked_component_list[[6]])
    comp.2 = cbind(concatPCA.result$linked_component_list[[1]], concatPCA.result$linked_component_list[[2]], concatPCA.result$linked_component_list[[5]], concatPCA.result$linked_component_list[[7]])
    comp.3 = cbind(concatPCA.result$linked_component_list[[1]], concatPCA.result$linked_component_list[[3]], concatPCA.result$linked_component_list[[4]], concatPCA.result$linked_component_list[[8]])
    comp.4 = cbind(concatPCA.result$linked_component_list[[1]], concatPCA.result$linked_component_list[[3]], concatPCA.result$linked_component_list[[5]], concatPCA.result$linked_component_list[[9]])
    comp.1 = svd(comp.1)$u
    comp.2 = svd(comp.2)$u
    comp.3 = svd(comp.3)$u
    comp.4 = svd(comp.4)$u

    pct.1 = sum((svd(comp.1 %*% (t(comp.1) %*% data_list[[1]]))$d[1 : 8])^2) / sum(svd(data_list[[1]])$d[1 : 8]^2)
    pct.2 = sum((svd(comp.2 %*% (t(comp.2) %*% data_list[[2]]))$d[1 : 8])^2) / sum(svd(data_list[[2]])$d[1 : 8]^2)
    pct.3 = sum((svd(comp.3 %*% (t(comp.3) %*% data_list[[3]]))$d[1 : 8])^2) / sum(svd(data_list[[3]])$d[1 : 8]^2)
    pct.4 = sum((svd(comp.4 %*% (t(comp.4) %*% data_list[[4]]))$d[1 : 8])^2) / sum(svd(data_list[[4]])$d[1 : 8]^2)
    out.1 = c(pct.1, pct.2, pct.3, pct.4)

    distance.space.1 = mean(svd(t(concatPCA.result$linked_component_list[[1]]) %*% configuration_setting$commonComponent)$d)
    distance.space.2 = mean(svd(t(concatPCA.result$linked_component_list[[2]]) %*% configuration_setting$partialComponent[[1]])$d)
    distance.space.3 = mean(svd(t(concatPCA.result$linked_component_list[[3]]) %*% configuration_setting$partialComponent[[2]])$d)
    distance.space.4 = mean(svd(t(concatPCA.result$linked_component_list[[4]]) %*% configuration_setting$partialComponent[[3]])$d)
    distance.space.5 = mean(svd(t(concatPCA.result$linked_component_list[[5]]) %*% configuration_setting$partialComponent[[4]])$d)
    distance.space.6 = mean(svd(t(concatPCA.result$linked_component_list[[6]]) %*% configuration_setting$individualComponent[[1]])$d)
    distance.space.7 = mean(svd(t(concatPCA.result$linked_component_list[[7]]) %*% configuration_setting$individualComponent[[2]])$d)
    distance.space.8 = mean(svd(t(concatPCA.result$linked_component_list[[8]]) %*% configuration_setting$individualComponent[[3]])$d)
    distance.space.9 = mean(svd(t(concatPCA.result$linked_component_list[[9]]) %*% configuration_setting$individualComponent[[4]])$d)
    out.distance.1 = c(distance.space.1, distance.space.2, distance.space.3, distance.space.4, distance.space.5, distance.space.6, distance.space.7, distance.space.8, distance.space.9)
    print("Concatenate Finished!")
    ############################################################
    ##
    ##
    ## LCA.2 result
    ##
    ##
    ############################################################
    comp.1 = svd(cbind(jointPCA.result$linked_component_list[[1]], jointPCA.result$linked_component_list[[2]], jointPCA.result$linked_component_list[[4]], jointPCA.result$linked_component_list[[6]]))$u
    comp.2 = svd(cbind(jointPCA.result$linked_component_list[[1]], jointPCA.result$linked_component_list[[2]], jointPCA.result$linked_component_list[[5]], jointPCA.result$linked_component_list[[7]]))$u
    comp.3 = svd(cbind(jointPCA.result$linked_component_list[[1]], jointPCA.result$linked_component_list[[3]], jointPCA.result$linked_component_list[[4]], jointPCA.result$linked_component_list[[8]]))$u
    comp.4 = svd(cbind(jointPCA.result$linked_component_list[[1]], jointPCA.result$linked_component_list[[3]], jointPCA.result$linked_component_list[[5]], jointPCA.result$linked_component_list[[9]]))$u
    pct.1 = sum((svd(comp.1 %*% (t(comp.1) %*% data_list[[1]]))$d[1 : 8])^2) / sum(svd(data_list[[1]])$d[1 : 8]^2)
    pct.2 = sum((svd(comp.2 %*% (t(comp.2) %*% data_list[[2]]))$d[1 : 8])^2) / sum(svd(data_list[[2]])$d[1 : 8]^2)
    pct.3 = sum((svd(comp.3 %*% (t(comp.3) %*% data_list[[3]]))$d[1 : 8])^2) / sum(svd(data_list[[3]])$d[1 : 8]^2)
    pct.4 = sum((svd(comp.4 %*% (t(comp.4) %*% data_list[[4]]))$d[1 : 8])^2) / sum(svd(data_list[[4]])$d[1 : 8]^2)
    out.2 = c(pct.1, pct.2, pct.3, pct.4)

    distance.space.1 = mean(svd(t(jointPCA.result$linked_component_list[[1]]) %*% configuration_setting$commonComponent)$d)
    distance.space.2 = mean(svd(t(jointPCA.result$linked_component_list[[2]]) %*% configuration_setting$partialComponent[[1]])$d)
    distance.space.3 = mean(svd(t(jointPCA.result$linked_component_list[[3]]) %*% configuration_setting$partialComponent[[2]])$d)
    distance.space.4 = mean(svd(t(jointPCA.result$linked_component_list[[4]]) %*% configuration_setting$partialComponent[[3]])$d)
    distance.space.5 = mean(svd(t(jointPCA.result$linked_component_list[[5]]) %*% configuration_setting$partialComponent[[4]])$d)
    distance.space.6 = mean(svd(t(jointPCA.result$linked_component_list[[6]]) %*% configuration_setting$individualComponent[[1]])$d)
    distance.space.7 = mean(svd(t(jointPCA.result$linked_component_list[[7]]) %*% configuration_setting$individualComponent[[2]])$d)
    distance.space.8 = mean(svd(t(jointPCA.result$linked_component_list[[8]]) %*% configuration_setting$individualComponent[[3]])$d)
    distance.space.9 = mean(svd(t(jointPCA.result$linked_component_list[[9]]) %*% configuration_setting$individualComponent[[4]])$d)
    out.distance.2 = c(distance.space.1, distance.space.2, distance.space.3, distance.space.4, distance.space.5, distance.space.6, distance.space.7, distance.space.8, distance.space.9)
    print("Joint Finished!")

    ############################################################
    ##
    ##
    ## LCA.3 result
    ##
    ##
    ############################################################
    comp.1 = svd(cbind(linkedPCA.result$linked_component_list[[1]], linkedPCA.result$linked_component_list[[2]], linkedPCA.result$linked_component_list[[4]], linkedPCA.result$linked_component_list[[6]]))$u
    comp.2 = svd(cbind(linkedPCA.result$linked_component_list[[1]], linkedPCA.result$linked_component_list[[2]], linkedPCA.result$linked_component_list[[5]], linkedPCA.result$linked_component_list[[7]]))$u
    comp.3 = svd(cbind(linkedPCA.result$linked_component_list[[1]], linkedPCA.result$linked_component_list[[3]], linkedPCA.result$linked_component_list[[4]], linkedPCA.result$linked_component_list[[8]]))$u
    comp.4 = svd(cbind(linkedPCA.result$linked_component_list[[1]], linkedPCA.result$linked_component_list[[3]], linkedPCA.result$linked_component_list[[5]], linkedPCA.result$linked_component_list[[9]]))$u
    pct.1 = sum((svd(comp.1 %*% (t(comp.1) %*% data_list[[1]]))$d[1 : 8])^2) / sum(svd(data_list[[1]])$d[1 : 8]^2)
    pct.2 = sum((svd(comp.2 %*% (t(comp.2) %*% data_list[[2]]))$d[1 : 8])^2) / sum(svd(data_list[[2]])$d[1 : 8]^2)
    pct.3 = sum((svd(comp.3 %*% (t(comp.3) %*% data_list[[3]]))$d[1 : 8])^2) / sum(svd(data_list[[3]])$d[1 : 8]^2)
    pct.4 = sum((svd(comp.4 %*% (t(comp.4) %*% data_list[[4]]))$d[1 : 8])^2) / sum(svd(data_list[[4]])$d[1 : 8]^2)
    out.3 = c(pct.1, pct.2, pct.3, pct.4)


    distance.space.1 = mean(svd(t(linkedPCA.result$linked_component_list[[1]]) %*% configuration_setting$commonComponent)$d)
    distance.space.2 = mean(svd(t(linkedPCA.result$linked_component_list[[2]]) %*% configuration_setting$partialComponent[[1]])$d)
    distance.space.3 = mean(svd(t(linkedPCA.result$linked_component_list[[3]]) %*% configuration_setting$partialComponent[[2]])$d)
    distance.space.4 = mean(svd(t(linkedPCA.result$linked_component_list[[4]]) %*% configuration_setting$partialComponent[[3]])$d)
    distance.space.5 = mean(svd(t(linkedPCA.result$linked_component_list[[5]]) %*% configuration_setting$partialComponent[[4]])$d)
    distance.space.6 = mean(svd(t(linkedPCA.result$linked_component_list[[6]]) %*% configuration_setting$individualComponent[[1]])$d)
    distance.space.7 = mean(svd(t(linkedPCA.result$linked_component_list[[7]]) %*% configuration_setting$individualComponent[[2]])$d)
    distance.space.8 = mean(svd(t(linkedPCA.result$linked_component_list[[8]]) %*% configuration_setting$individualComponent[[3]])$d)
    distance.space.9 = mean(svd(t(linkedPCA.result$linked_component_list[[9]]) %*% configuration_setting$individualComponent[[4]])$d)
    out.distance.3 = c(distance.space.1, distance.space.2, distance.space.3, distance.space.4, distance.space.5, distance.space.6, distance.space.7, distance.space.8, distance.space.9)

    print("Linked Finished!")
    ############################################################
    ##
    ##
    ## Output result
    ##
    ##
    ############################################################
    out.dat = data.frame(out.1, out.2, out.3)
    out.distance.dat = data.frame(out.distance.1, out.distance.2, out.distance.3)
    # names(out.dat) = c("lca.biconvex", "lca.concat", "lca.covariance")
    # save(out.dat, file = paste0('./simulation/output1/output_', i, ".RData"))
    # print(i)
    print(list(pctExplained = out.dat, angleDistance = out.distance.dat))
}



