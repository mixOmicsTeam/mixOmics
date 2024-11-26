#load data
data(nutrimouse)
X <- nutrimouse$lipid
Y <- nutrimouse$gene

# run tuning
tune_res <- tune.rcc(X, Y, validation = "Mfold")

# plot output
plot(tune_res)