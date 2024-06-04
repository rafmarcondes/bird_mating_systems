library(phytools)
setwd('/Users/rs155/Dropbox/my_research/mating-systems/REVISION/unrandomized_analyses')

s <- readRDS("simmap_muhisse_500sims.rds")
ds <- readRDS("dsimmap_muhisse_500sims.rds")

#### avg n of mating system transitions####
ds$count <- as.data.frame(ds$count)

ds$count$MP <- ds$count$`00A,01A` + ds$count$`00B,01B`
ds$count$ML <- ds$count$`00A,10A` + ds$count$`00B,10B`
ds$count$PM <- ds$count$`01A,00A` + ds$count$`01B,00B`
ds$count$PL <- ds$count$`01A,10A` + ds$count$`01B,10B`
ds$count$LP <- ds$count$`10A,01A` + ds$count$`10B,01B`
ds$count$LM <- ds$count$`10A,00A` + ds$count$`10B,00B`

avg <- as.data.frame(colMeans(ds$count[, 32:37]))
avg$sd <- apply(X = (ds$count[, 32:37]), FUN = sd, MARGIN = 2)


######## plot traits w colored edges#########
ace <- as.data.frame(ds$ace)
ace$M <- ace$`00A` + ace$`00B`
ace$P <- ace$`01A` + ace$`01B`
ace$L <- ace$`10A` + ace$`10B`
ace <- ace[, 7:9]

#### FIG 1####
#pdf("simmap_muhisse_with_colored_edges.pdf", width = 100, height = 100)
png("simmap_muhisse_with_colored_edges.png",units='px',width=1000,height=1000)

# code below adapted from http://blog.phytools.org/2023/04/coloring-edges-of-plotted-tree-by.html
tip.liks <- ace[(6619 + 1):13239, ]
node.liks <- ace[1:6619, ]
allStates <- rbind(tip.liks, node.liks)
cols <- setNames(
  object = c("black", "yellow", "red"),
  nm = colnames(ace)
)
edge.cols <- cols[apply(allStates, 1, function(x) {
  which(x == max(x))
})][s[[1]]$edge[, 2]]
plot.phylo(
  x = s[[1]], cex = 0.2, edge.width = 0.01, type = "fan",
  edge.lty = 3, edge.color = edge.cols
)
legend(
  x = "bottomleft", legend = c("Monogamous", "Resource-defense polygamous", "Lek"),
  col = c("black", "yellow", "red"), cex = 10, pch = 15
)
dev.off()

######## plot traits w pies#########

pdf("muhisse_simmap_with_pie_charts.pdf", width = 100, height = 100)
plot.phylo(
  x = s[[1]], cex = 0.2, edge.width = 0.01, type = "fan",
  edge.lty = 3, edge.color = edge.cols
)
cols <- make.transparent(cols, 0.5)
nodelabels(pie = ace, cex = 0.4, piecol = cols, bg = cols)
tiplabels(pie = ace, cex = 0.05, offset = 1, pch = 20, piecol = cols, bg = cols)
dev.off()


######## plot hidden states with colored edges#########

ace <- as.data.frame(ds$ace)
ace$A <- ace$`00A` + ace$`01A` + ace$`10A`
ace$B <- ace$`00B` + ace$`01B` + ace$`10B`

ace <- ace[, 7:8] 
colnames(ace)=c('B','A')#rename regimes so that A is the background

pdf("simmap_muhisse_hidden_regimes_with_colored_edges.pdf", width = 100, height = 100)
# par(bg='grey')
tip.liks <- ace[(6619 + 1):13239, ]
node.liks <- ace[1:6619, ]
allStates <- rbind(tip.liks, node.liks)
# allStates=as.matrix(allStates)
allStates$A <- unlist(allStates$A)
allStates$B <- unlist(allStates$B)
cols <- setNames(
  object = c("green", "purple"),
  nm = colnames(ace)
)
dummyfunction <- function(x) {
  return(which(x == max(x)))
}
bla <- apply(allStates, 1, dummyfunction)
edge.cols <- cols[unlist(bla)][s[[1]]$edge[, 2]]
plot.phylo(
  x = s[[1]], cex = 0.2, edge.width = 0.01, type = "fan",
  edge.lty = 3, edge.color = edge.cols, show.tip.label = F
)
legend(
  x = "bottomleft", legend = c("Regime A", "Regime B"),
  col = c("purple", "green"), cex = 20, pch = 15
)
dev.off()
