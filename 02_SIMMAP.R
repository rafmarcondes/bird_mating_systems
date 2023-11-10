#### IN THIS VERSION OF THE SCRIPT THE STATES FOR SL, EL AND CL
#### ARE ALL BINNED INTO CL, SO EL AND SL ARE EXCLUDED (LINE 159).
# THE PRIOR PROBABILITIES ARE SET
#### ACCORDINGLY


setwd("/Users/rs155/Dropbox/my_research/mating-systems/")
library(phytools)
library(geiger)

tr <- read.tree("MCC_of_100_hackett_trees.tre")
# tr=keep.tip(phy=tr,tip=sample(tr$tip.label,100)) #prune to X random tips

dat <- read.csv("/Users/rs155/Dropbox/my_research/mating-systems/MS_DATASET.csv")
dat$Mating_system <- gsub(x = dat$Mating_system, pattern = " ", replacement = "")
rownames(dat) <- dat$Species
dat$Species <- NULL
rownames(dat) <- gsub(x = rownames(dat), pattern = " ", replacement = "_")

td <- treedata(phy = tr, data = dat, sort = T)
tr <- td$phy
dat <- td$data

pdf("test.pdf", height = 21, widt = 21)
plot(tr, type = "fan")
dev.off()

# convert dat to matrix of prior probabilities
# order of states: M, P, SL, EL, CL

priormatrix <- as.data.frame(dat)
priormatrix <- priormatrix[c("Mating_system", "Data_quality")]

# replace NA data quality (U mating) with -1
priormatrix$Data_quality <- replace(x = priormatrix$Data_quality, list = is.na(priormatrix$Data_quality), values = -1)


priormatrix$M <- rep(NA, nrow(dat))
priormatrix$P <- rep(NA, nrow(dat))
priormatrix$SL <- rep(NA, nrow(dat))
priormatrix$EL <- rep(NA, nrow(dat))
priormatrix$CL <- rep(NA, nrow(dat))
priormatrix$Data_quality <- as.numeric(priormatrix$Data_quality)


# priors for quality 3
mprior <- c(1, 0, 0, 0, 0)
pprior <- c(0, 1, 0, 0, 0)
# slprior=c(0,0,1,0,0)
# elprior=c(0,0,0,1,0)
clprior <- c(0, 0, 0, 0, 1)
slprior <- clprior
elprior <- clprior

# add priors to sp with quality 3
for (sp in rownames(priormatrix)) {
  cat(sp, priormatrix[sp, "Mating_system"], priormatrix[sp, "Data_quality"], "\n")
  if (priormatrix[sp, "Data_quality"] == 3) {
    if (priormatrix[sp, "Mating_system"] == "M") {
      priormatrix[sp, c("M", "P", "SL", "EL", "CL")] <- mprior
    }
    if (priormatrix[sp, "Mating_system"] == "P") {
      priormatrix[sp, c("M", "P", "SL", "EL", "CL")] <- pprior
    }
    if (priormatrix[sp, "Mating_system"] == "SL") {
      priormatrix[sp, c("M", "P", "SL", "EL", "CL")] <- slprior
    }
    if (priormatrix[sp, "Mating_system"] == "EL") {
      priormatrix[sp, c("M", "P", "SL", "EL", "CL")] <- elprior
    }
    if (priormatrix[sp, "Mating_system"] == "CL") {
      priormatrix[sp, c("M", "P", "SL", "EL", "CL")] <- clprior
    }
  }
}

# priors for quality 2
mprior <- c(0.95, 0.025, 0.025, 0.025, 0.025)
pprior <- c(0.025, 0.95, 0.025, 0, 0.025)
# slprior=c(0.0125,0.0125,0.95,0.0125,0.0125)
# elprior=c(0.0125,0.0125,0.0125,0.95,0.0125)
clprior <- c(0.025, 0.025, 0.025, 0.025, 0.95)
slprior <- clprior
elprior <- clprior

# add priors to sp with quality 2
for (sp in rownames(priormatrix)) {
  cat(sp, priormatrix[sp, "Mating_system"], priormatrix[sp, "Data_quality"], "\n")
  if (priormatrix[sp, "Data_quality"] == 2) {
    if (priormatrix[sp, "Mating_system"] == "M") {
      priormatrix[sp, c("M", "P", "SL", "EL", "CL")] <- mprior
    }
    if (priormatrix[sp, "Mating_system"] == "P") {
      priormatrix[sp, c("M", "P", "SL", "EL", "CL")] <- pprior
    }
    if (priormatrix[sp, "Mating_system"] == "SL") {
      priormatrix[sp, c("M", "P", "SL", "EL", "CL")] <- slprior
    }
    if (priormatrix[sp, "Mating_system"] == "EL") {
      priormatrix[sp, c("M", "P", "SL", "EL", "CL")] <- elprior
    }
    if (priormatrix[sp, "Mating_system"] == "CL") {
      priormatrix[sp, c("M", "P", "SL", "EL", "CL")] <- clprior
    }
  }
}


# priors for quality 1
mprior <- c(0.9, 0.05, 0.025, 0.025, 0.05)
pprior <- c(0.05, 0.9, 0.5, 0.5, 0.05)
# slprior=c(0.025,0.025,0.9,0.025,0.025)
# elprior=c(0.025,0.025,0.025,0.9,0.025)
clprior <- c(0.05, 0.05, 0.05, 0.05, 0.9)
slprior <- clprior
elprior <- clprior

# add priors to sp with quality 1
for (sp in rownames(priormatrix)) {
  cat(sp, priormatrix[sp, "Mating_system"], priormatrix[sp, "Data_quality"], "\n")
  if (priormatrix[sp, "Data_quality"] == 1) {
    if (priormatrix[sp, "Mating_system"] == "M") {
      priormatrix[sp, c("M", "P", "SL", "EL", "CL")] <- mprior
    }
    if (priormatrix[sp, "Mating_system"] == "P") {
      priormatrix[sp, c("M", "P", "SL", "EL", "CL")] <- pprior
    }
    if (priormatrix[sp, "Mating_system"] == "SL") {
      priormatrix[sp, c("M", "P", "SL", "EL", "CL")] <- slprior
    }
    if (priormatrix[sp, "Mating_system"] == "EL") {
      priormatrix[sp, c("M", "P", "SL", "EL", "CL")] <- elprior
    }
    if (priormatrix[sp, "Mating_system"] == "CL") {
      priormatrix[sp, c("M", "P", "SL", "EL", "CL")] <- clprior
    }
  }
}


# priors for quality 0
mprior <- c(0.85, 0.075, 0.075, 0.075, 0.075)
pprior <- c(0.075, 0.85, 0.075, 0.075, 0.075)
# slprior=c(0.0375,0.0375,0.85,0.0375,0.0375)
# elprior=c(0.0375,0.0375,0.0375,0.85,0.0375)
clprior <- c(0.075, 0.075, 0.075, 0.075, 0.85)
slprior <- clprior
elprior <- clprior

# add priors to sp with quality 0
for (sp in rownames(priormatrix)) {
  cat(sp, priormatrix[sp, "Mating_system"], priormatrix[sp, "Data_quality"], "\n")
  if (priormatrix[sp, "Data_quality"] == 0) {
    if (priormatrix[sp, "Mating_system"] == "M") {
      priormatrix[sp, c("M", "P", "SL", "EL", "CL")] <- mprior
    }
    if (priormatrix[sp, "Mating_system"] == "P") {
      priormatrix[sp, c("M", "P", "SL", "EL", "CL")] <- pprior
    }
    if (priormatrix[sp, "Mating_system"] == "SL") {
      priormatrix[sp, c("M", "P", "SL", "EL", "CL")] <- slprior
    }
    if (priormatrix[sp, "Mating_system"] == "EL") {
      priormatrix[sp, c("M", "P", "SL", "EL", "CL")] <- elprior
    }
    if (priormatrix[sp, "Mating_system"] == "CL") {
      priormatrix[sp, c("M", "P", "SL", "EL", "CL")] <- clprior
    }
  }
}

# priors for completely U will be flat

for (sp in rownames(priormatrix)) {
  cat(sp, priormatrix[sp, "Mating_system"], priormatrix[sp, "Data_quality"], "\n")
  if (priormatrix[sp, "Mating_system"] == "U") {
    priormatrix[sp, c("M", "P", "SL", "EL", "CL")] <- rep(1 / 3, 5)
  }
}


# priors for PU will be flat for P/L and 0 for M

puprior <- c(0, 0.5, 0.5, 0.5, 0.5)

for (sp in rownames(priormatrix)) {
  cat(sp, priormatrix[sp, "Mating_system"], priormatrix[sp, "Data_quality"], "\n")
  if (priormatrix[sp, "Mating_system"] == "PU") {
    priormatrix[sp, c("M", "P", "SL", "EL", "CL")] <- puprior
  }
}


# priors for LU will be flat for lekking and 0 for M and P

luprior <- c(0, 0, 1, 1, 1)

for (sp in rownames(priormatrix)) {
  cat(sp, priormatrix[sp, "Mating_system"], priormatrix[sp, "Data_quality"], "\n")
  if (priormatrix[sp, "Mating_system"] == "LU") {
    priormatrix[sp, c("M", "P", "SL", "EL", "CL")] <- luprior
  }
}



priormatrix <- priormatrix[, c("M", "P", "CL")]


### simmap####
# priormatrix=priormatrix[,3:7]
td <- treedata(phy = tr, data = priormatrix, sort = T)
tr <- td$phy
priormatrix <- td$data

Sys.time()
m.ard <- make.simmap(tree = tr, x = priormatrix, model = "ARD", nsim = 1000, pi = "estimated")
Sys.time()
md.ard <- describe.simmap(m.ard)

c.ard <- describe.simmap(m.ard) # n of each type of transition across all stochastic maps
means.c.ard <- colMeans(c.ard$count)
means.c.ard <- round(x = means.c.ard, digits = 1)

saveRDS(object = m.ard, file = "simmap_binned_lek_15may2023.rds")

###### compare tip state likelihoods to priors
colnames(md.ard$tips) <- paste0("post_", colnames(md.ard$tips))
colnames(priormatrix) <- paste0("prior_", colnames(priormatrix))
comparison <- merge(x = priormatrix, y = md.ard$tips, by = 0)
rownames(comparison) <- comparison$Row.names


### DF of imputed states####

imputed <- vector(length = nrow(tip.liks))
names(imputed) <- rownames(tip.liks)
for (t in 1:nrow(tip.liks)) {
  s <- which.max(tip.liks[t, ])
  imputed[t] <- names(s)
}
imputed <- as.data.frame(imputed)

write.csv(x = imputed, file = "imputed_data_binned_lek_15may2023.csv")
