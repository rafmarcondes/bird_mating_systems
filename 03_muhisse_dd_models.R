library(hisse)

setwd("/Users/rs155/Dropbox/my_research/mating-systems/")
library(corHMM)
library(geiger)

tr <- read.tree("MCC_of_100_hackett_trees.tre")
dat <- read.csv("imputed_data_binned_lek_15may2023.csv", row.names = 1)

td <- treedata(phy = tr, data = dat, sort = F)
tr <- td$phy
dat <- as.data.frame(td$data)
dat[, 2] <- dat[, 1]
dat[, 1] <- rownames(dat)

# transform dat so that M=00, P=01 and L=10

dat$T1 <- rep(NA, nrow(dat))
dat$T2 <- rep(NA, nrow(dat))

for (sp in 1:nrow(dat)) {
  if (dat[sp, 2] == "M") {
    dat[sp, 3] <- 0
    dat[sp, 4] <- 0
  }
  if (dat[sp, 2] == "P") {
    dat[sp, 3] <- 0
    dat[sp, 4] <- 1
  }
  if (grepl(x = dat[sp, 2], pattern = "L") == T) {
    dat[sp, 3] <- 1
    dat[sp, 4] <- 0
  }
}
dat <- dat[, c(1, 3, 4)]
colnames(dat)[1] <- "sp"


setwd("/Users/rs155/Dropbox/my_research/mating-systems/muhisse_models_divers_dep_on_hidden_states_02jun2023")

#### FD1####

# diversific params will be different for each trait state
turnover <- c(1, 2, 3, 0) # last has to be zero bc it corresponds to the 11 state
extinction.fraction <- c(1, 2, 3, 0)
f <- rep(0.66, 4) # sampling fraction

trans.rate <- TransMatMakerMuHiSSE(hidden.traits = 0, make.null = T)
trans.rate[!is.na(trans.rate)] <- 1:12 # make dual transitions possible

# drop 11 from the matrix
trans.rate <- ParDrop(trans.rate, c(3, 6, 9, 10, 11, 12))

Sys.time()
FD1 <- MuHiSSE(
  phy = tr, data = dat, f = f, turnover = turnover,
  eps = extinction.fraction, hidden.states = F,
  trans.rate = trans.rate
)
Sys.time()
saveRDS(FD1, "FD1.rds")


#### C1D1####
# any transitions back into monogamy forbidden

trans.rate <- ParDrop(trans.rate, c(1, 2))

Sys.time()
C1D1 <- MuHiSSE(
  phy = tr, data = dat, f = f, turnover = turnover,
  eps = extinction.fraction, hidden.states = F,
  trans.rate = trans.rate
)
Sys.time()
saveRDS(C1D1, "C1D1.rds")


#### C2D1####
# same as C1, but  transitions from L into P also forbidden

trans.rate <- ParDrop(trans.rate, 2)

Sys.time()
C2D1 <- MuHiSSE(
  phy = tr, data = dat, f = f, turnover = turnover,
  eps = extinction.fraction, hidden.states = F,
  trans.rate = trans.rate
)
Sys.time()
saveRDS(C2D1, "C2D1.rds")



#### C3D1####
# same as C2, but can not go straight from M to  L

trans.rate <- ParDrop(trans.rate, 2)

Sys.time()
C3D1 <- MuHiSSE(
  phy = tr, data = dat, f = f, turnover = turnover,
  eps = extinction.fraction, hidden.states = F,
  trans.rate = trans.rate
)
Sys.time()
saveRDS(C3D1, "C3D1.rds")

#### FD2####
trans.rate <- TransMatMakerMuHiSSE(hidden.traits = 1, make.null = T, include.diagonals = T)
trans.rate[!is.na(trans.rate)] <- 1:32 # make dual transitions possible
# the way this is set up, it doesn't allow simultanous transitions in the observed trait
# and in the hidden trait

# drop 11 from the matrix
trans.rate <- ParDrop(
  rate.mat = trans.rate,
  drop.par = c(3, 7, 11, 13:16, 20, 24, 28:32)
)

turnover <- c(1, 2, 3, 0, 4, 5, 6, 0) # zeroes correspond to the 11 state
extinction.fraction <- c(1, 2, 3, 0, 4, 5, 6, 0)
f <- rep(0.66, 4) # sampling fraction

Sys.time()
FD2 <- MuHiSSE(
  phy = tr, data = dat, f = f, turnover = turnover,
  eps = extinction.fraction, hidden.states = T,
  trans.rate = trans.rate
)
Sys.time()
saveRDS(FD2, "FD2.rds")


#### C1D2####
# transitions to monogamy forbidden + 2 hidden states

trans.rate <- ParDrop(trans.rate, c(1, 2, 11, 12))

Sys.time()
C1D2 <- MuHiSSE(
  phy = tr, data = dat, f = f, turnover = turnover,
  eps = extinction.fraction, hidden.states = T,
  trans.rate = trans.rate
)
Sys.time()
saveRDS(C1D2, "C1D2.rds")


#### C2D2####
# same as C2 + transitions from L to P forbidden
#+ 2 hidden states

trans.rate <- ParDrop(trans.rate, c(3, 11))

Sys.time()
C2D2 <- MuHiSSE(
  phy = tr, data = dat, f = f, turnover = turnover,
  eps = extinction.fraction, hidden.states = T,
  trans.rate = trans.rate
)
Sys.time()
saveRDS(C2D2, "C2D2.rds")

#### C3D2####
# same as C2 + cannot go straight from M to  L
#+ 2 hidden states

trans.rate <- ParDrop(trans.rate, c(4, 11))

Sys.time()
C3D2 <- MuHiSSE(
  phy = tr, data = dat, f = f, turnover = turnover,
  eps = extinction.fraction, hidden.states = T,
  trans.rate = trans.rate
)
Sys.time()
saveRDS(C3D2, "C3D2.rds")



#### FD3####
trans.rate <- TransMatMakerMuHiSSE(hidden.traits = 2, make.null = T, include.diagonals = T)
trans.rate[!is.na(trans.rate)] <- 1:60 # make dual transitions possible
# the way this is set up, it doesn't allow simultaneous transitions in the observed trait
# and in the hidden trait

# drop 11 from the matrix
trans.rate <- ParDrop(
  rate.mat = trans.rate,
  drop.par = c(
    16:20, 37:40, 56:60, 3, 8, 13, 36, 56, 19,
    24, 29, 34, 57, 20, 40, 45, 50, 55
  )
)

# diversification params will be different for each trait state, but same across hidden states
turnover <- c(1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0) # zeroes correspond to the 11 state
extinction.fraction <- c(1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0)
f <- rep(0.66, 4) # sampling fraction

Sys.time()
FD3 <- MuHiSSE(
  phy = tr, data = dat, f = f, turnover = turnover,
  eps = extinction.fraction, hidden.states = T,
  trans.rate = trans.rate
)
Sys.time()
saveRDS(FD3, "FD3.rds")


#### C1D3####
# transitions to monogamy forbidden + 3 hidden states

trans.rate <- ParDrop(trans.rate, c(1:4, 13:16, 25:28))

Sys.time()
C1D3 <- MuHiSSE(
  phy = tr, data = dat, f = f, turnover = turnover,
  eps = extinction.fraction, hidden.states = T,
  trans.rate = trans.rate
)
Sys.time()
saveRDS(C1D3, "C1D3.rds")


#### C2D3####
# same as C2 + transitions from L to P forbidden
#+ 3 hidden states

trans.rate <- ParDrop(trans.rate, c(2, 11, 20))

Sys.time()
C2D3 <- MuHiSSE(
  phy = tr, data = dat, f = f, turnover = turnover,
  eps = extinction.fraction, hidden.states = T,
  trans.rate = trans.rate
)
Sys.time()
saveRDS(C2D3, "C2D3.rds")

#### C3D3####
# same as C2 + cannot go straight from M to  L
#+ 3 hidden states

trans.rate <- ParDrop(trans.rate, c(4, 12, 20))

Sys.time()
C3D3 <- MuHiSSE(
  phy = tr, data = dat, f = f, turnover = turnover,
  eps = extinction.fraction, hidden.states = T,
  trans.rate = trans.rate
)
Sys.time()
saveRDS(C3D3, "C3D3.rds")
