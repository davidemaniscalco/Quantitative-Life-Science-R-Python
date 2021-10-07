Items_List <- c(8,16,32,64,128,256,512)
Repetitions <- 10000
starting_item <- 1
#Recalled_Matrix <- matrix(0, nrow=Repetitions, ncol=length(Items_List))
ff = 0.01
NN = 1000000
verbose = TRUE
each = 50
symmetric = FALSE
single = TRUE
with_ff = TRUE                # if false, random matrices with runif are created, either symm or asymm
filename = './files/19giu_asymmetric_simple_001.csv'
