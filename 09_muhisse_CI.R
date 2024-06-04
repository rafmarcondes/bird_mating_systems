library(hisse)

setwd('/Users/rs155/Dropbox/my_research/mating-systems/REVISION/unrandomized_analyses')

#top model
mod=readRDS('muhisse_models_divers_dep_on_hidden_states/FD2.rds')

Sys.time()
CI=SupportRegionMuHiSSE(mod,scale.int=0.005,desired.delta=4,verbose=T,n.points=20000)
Sys.time()

saveRDS(object=CI, file='muhisse_CI.rds')