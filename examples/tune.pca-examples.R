# load data
data(liver.toxicity)

# run tuning
tune <- tune.pca(liver.toxicity$gene, center = TRUE, scale = TRUE)
plot(tune)

# set up multilevel dataset
repeat.indiv <- c(1, 2, 1, 2, 1, 2, 1, 2, 3, 3, 4, 3, 4, 3, 4, 4, 5, 6, 5, 5,
                  6, 5, 6, 7, 7, 8, 6, 7, 8, 7, 8, 8, 9, 10, 9, 10, 11, 9, 9,
                  10, 11, 12, 12, 10, 11, 12, 11, 12, 13, 14, 13, 14, 13, 14,
                  13, 14, 15, 16, 15, 16, 15, 16, 15, 16)
design <- data.frame(sample = repeat.indiv)

# run tuning
tune <- tune.pca(liver.toxicity$gene, center = TRUE, scale = TRUE, multilevel = design)
plot(tune)