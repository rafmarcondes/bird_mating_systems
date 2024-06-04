library(hisse)
library(RColorBrewer)
setwd("/Users/rs155/Dropbox/my_research/mating-systems/REVISION/unrandomized_analyses")

# top model
mod <- readRDS("muhisse_models_divers_dep_on_hidden_states/FD2.rds")

#### reconstruction####

Sys.time()
reco <- MarginReconMuHiSSE(
  phy = mod$phy,
  data = mod$data,
  f = mod$f,
  pars = mod$solution,
  hidden.states = 2,
  n.cores = 5
)
Sys.time()
saveRDS(object = reco, file = "muhisse_reco.rds")


reco <- readRDS("/Users/rs155/Dropbox/my_research/mating-systems/muhisse_reco_07jun2023.rds")
reco$node.mat <- reco$node.mat[, c(2:4, 6:8)]


### plot combined states###
pal <- setNames(nm = colnames(reco$node.mat), object = brewer.pal(6, "Accent"))
pdf("combined_states_muhisse.pdf", width = 20, height = 60)
plot.phylo(
  x = mod$phy, cex = 0.1, edge.width = 0.5,
  edge.lty = 1, bg = "grey", direction = "leftwards"
)
nodelabels(pie = reco$node.mat, cex = 0.1, piecol = pal)
legend(x = "topright", legend = names(pal), col = pal, fill = pal)

dev.off()


### plot hidden states###
hidden <- as.data.frame(reco$node.mat)
hidden$A <- hidden$`(00A)` + hidden$`(01A)` + hidden$`(10A)`
hidden$B <- hidden$`(00B)` + hidden$`(01B)` + hidden$`(10B)`
hidden <- hidden[, c("A", "B")]

pdf("hidden_states_muhisse.pdf", width = 20, height = 60)
pal <- setNames(nm = c("A", "B"), object = c("red", "blue"))
plot.phylo(
  x = mod$phy, cex = 0.1, edge.width = 0.5,
  edge.lty = 1, bg = "grey", direction = "leftwards"
)
nodelabels(pie = hidden, cex = 0.1, piecol = pal)
legend(x = "topright", legend = names(pal), col = pal, fill = pal)
dev.off()

# plot trait states
traits <- as.data.frame(reco$node.mat)
traits$m <- traits$`(00A)` + traits$`(00B)`
traits$p <- traits$`(01A)` <- traits$`(01B)`
traits$l <- traits$`(10A)` <- traits$`(10B)`

traits <- traits[, c("m", "p", "l")]

pdf("trait_states_muhisse.pdf", width = 20, height = 60)
pal <- setNames(nm = c("m", "p", "l"), object = c("red", "blue", "white"))
plot.phylo(
  x = mod$phy, cex = 0.1, edge.width = 0.5,
  edge.lty = 1, bg = "grey", direction = "leftwards"
)
nodelabels(pie = traits, cex = 0.1, piecol = pal)
legend(x = "topright", legend = names(pal), col = pal, fill = pal)
dev.off()
