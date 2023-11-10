library(ape)
library(phangorn)
setwd("/Users/rs155/Dropbox/my_research/mating-systems/")

tr <- read.tree("/Users/rs155/Dropbox/my_research/mating-systems/HackettStage1Full_3.tre")

bla <- sample(x = tr, size = 100) # subset to only 100

mcctr <- read.tree("MCC_of_100_hackett_trees.tre")


rez <- setNames(object = rep(NA, 100), nm = 1:100)
for (tr in 1:100) {
  rez[tr] <- wRF.dist(tree1 = mcctr, tree2 = bla[[tr]], rooted = T)
  cat(tr, "\n")
}

most_different_tree <- bla[which.max(rez)][[1]]
midpoint_tree <- bla[which.min(abs(rez - quantile(x = rez, probs = 0.50)))][[1]]

write.tree(phy = most_different_tree, file = "most_different_tree.tre")
write.tree(phy = midpoint_tree, file = "midpoint_tree.tre")
