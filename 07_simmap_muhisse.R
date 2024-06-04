library(hisse)
library(phytools)
setwd("/Users/rs155/Dropbox/my_research/mating-systems/REVISION/unrandomized_analyses")

#### set up Q-matrix based on muhisse FD2 model estimates####

muhisse <- readRDS("muhisse_models_divers_dep_on_hidden_states/FD2.rds")$solution
muhisse <- muhisse[grep(x = names(muhisse), pattern = "q")]
muhisse <- muhisse[muhisse > 0.0000000001]
names(muhisse) <- gsub(x = names(muhisse), pattern = "q", replacement = "")
q <- matrix(ncol = 6, nrow = 6, data = rep(0, 36))
colnames(q) <- c("00A", "01A", "10A", "00B", "01B", "10B")
rownames(q) <- colnames(q)
q <- as.data.frame(q)

for (p in names(muhisse)) {
  from <- substr(x = p, start = 1, stop = 3)
  to <- substr(x = p, start = 5, stop = 7)
  q[from, to] <- muhisse[p]
}

q <- as.matrix(q)

diag(q) <- (-rowSums(q))



#### set up matrix of tip probabilities based on muhisse recosntruction####

reco <- readRDS("muhisse_reco.rds")
dat <- reco$tip.mat
dat <- dat[, c(2:4, 6:8)]
colnames(dat) <- gsub(x = colnames(dat), pattern = "\\(", replacement = "")
colnames(dat) <- gsub(x = colnames(dat), pattern = "\\)", replacement = "")
rownames(dat) <- reco$phy$tip.label
dat <- as.matrix(dat)

#### run simmap####
Sys.time()
sm <- make.simmap(tree = reco$phy, x = dat, Q = q, nsim = 500, pi = "estimated")
Sys.time()
saveRDS(object = sm, file = "simmap_muhisse_500sims.rds")
smd <- describe.simmap(sm)
Sys.time()
saveRDS(object = smd, file = "dsimmap_muhisse_500sims.rds")
