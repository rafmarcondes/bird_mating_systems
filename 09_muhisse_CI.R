library(hisse)

#top model
mod=readRDS('/Users/rs155/Dropbox/my_research/mating-systems/muhisse_models_divers_dep_on_hidden_states_02jun2023/FD2.rds')

Sys.time()
CI=SupportRegionMuHiSSE(mod,scale.int=0.005,desired.delta=4,verbose=T,n.points=20000)
Sys.time()

saveRDS(object=CI, file='/Users/rs155/Dropbox/my_research/mating-systems/muhisse_CI_09oct2023.rds')