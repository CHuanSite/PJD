## Human Data Set
data_human = data.frame(matrix(rnorm(50), nrow = 10, ncol = 5))
rownames(data_human) = c("Gene-01", "Gene-02", "Gene-03", "Gene-04","Gene-05",
                    "Gene-06", "Gene-07", "Gene-08", "Gene-09", "Gene-10")

colnames(data_human) = c("Human-01", "Human-02", "Human-03", "Human-04", "Human-05")

## Mouse Data Set
data_mouse = data.frame(matrix(rnorm(50), nrow = 10, ncol = 5))
rownames(data_mouse) = c("")
