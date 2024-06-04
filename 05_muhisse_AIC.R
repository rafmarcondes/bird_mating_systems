library(hisse)

setwd("/Users/rs155/Dropbox/my_research/mating-systems/REVISION/unrandomized_analyses/muhisse_models_same_divers_across_hidden_states")
files <- list.files()
mods <- list()
for (f in files) {
  mods[length(mods) + 1] <- readRDS(f)$AICc
  names(mods)[length(mods)] <- gsub(".rds", "", f)
}

setwd("/Users/rs155/Dropbox/my_research/mating-systems/REVISION/unrandomized_analyses/muhisse_models_divers_dep_on_hidden_states")
files <- list.files()
for (f in files) {
  aic <- readRDS(f)$AICc
  mods[gsub(".rds", "", f)] <- aic
}

mods <- as.data.frame(unlist(mods))

write.csv(x=mods,file="/Users/rs155/Dropbox/my_research/mating-systems/REVISION/unrandomized_analyses/AIC.csv")

