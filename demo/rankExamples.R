## Simulation part of the algorithm
configuration_setting = configuration_setting_generation(featureNum = 500,
                                                         DataNum = c(50, 50, 50, 50),
                                                         noiseVariance = c(1, 1, 1, 1))
data_list = simulated_data_generation(configuration_setting,
                                      amplitude =  12,
                                      heterogeneousNoise = FALSE)

BEMA(svd(data_list[[2]])$d^2 / 50, p = 500, n = 50, alpha = 0.2, beta = 0.1)
