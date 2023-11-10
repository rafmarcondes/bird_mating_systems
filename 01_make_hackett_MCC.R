library(ape)
library(phangorn)
tr <- read.tree("/Users/rs155/Dropbox/my_research/mating-systems/HackettStage1Full_3.tre")

bla <- sample(x = tr, size = 100) # subset to only 100 trees 

mcctr <- maxCladeCred(x = bla) # mcc tree
write.tree(phy = mcctr, file = "/Users/rs155/Dropbox/my_research/mating-systems/MCC_of_100_hackett_trees.tre")
