library(hisse)

setwd("/Users/rs155/Dropbox/my_research/mating-systems/REVISION/unrandomized_analyses")
mod <- readRDS("muhisse_models_divers_dep_on_hidden_states/FD2.rds")
ci <- readRDS("muhisse_CI.rds")
ci <- ci$ci

params <- mod$solution
params <- params[params > 10**-100]
names(params) <- gsub(x = names(params), pattern = "00", replacement = "M")
names(params) <- gsub(x = names(params), pattern = "01", replacement = "P")
names(params) <- gsub(x = names(params), pattern = "10", replacement = "L")
names(params) <- gsub(x = names(params), pattern = "q", replacement = "")

colnames(ci) <- gsub(x = colnames(ci), pattern = "00", replacement = "M")
colnames(ci) <- gsub(x = colnames(ci), pattern = "01", replacement = "P")
colnames(ci) <- gsub(x = colnames(ci), pattern = "10", replacement = "L")
colnames(ci) <- gsub(x = colnames(ci), pattern = "q", replacement = "")
ci <- ci[, -grep(x = colnames(ci), pattern = "11")]

#### transitions A####
transitionsA <- params[grep(x = names(params), pattern = "A_[A-Z]A")]
names(transitionsA) <- gsub(x = names(transitionsA), pattern = "q", replacement = "")
names(transitionsA) <- gsub(x = names(transitionsA), pattern = "A", replacement = "")
CItransitionsA <- ci[, grep(x = colnames(ci), pattern = "A_[A-Z]A")]
colnames(CItransitionsA) <- gsub(x = colnames(CItransitionsA), pattern = "A", replacement = "")
CItransitionsA <- CItransitionsA[c(1, 5), ]
bp <- barplot(transitionsA, ylim = c(0, 0.1), main = "Transition rates, regime A")
lines(lwd = 1, x = rep(bp[1, 1], 2), y = CItransitionsA[, 1], col = "blue")
lines(lwd = 1, x = rep(bp[2, 1], 2), y = CItransitionsA[, 2], col = "blue")
lines(lwd = 1, x = rep(bp[3, 1], 2), y = CItransitionsA[, 3], col = "blue")
lines(lwd = 1, x = rep(bp[4, 1], 2), y = CItransitionsA[, 4], col = "blue")
lines(lwd = 1, x = rep(bp[5, 1], 2), y = CItransitionsA[, 5], col = "blue")
lines(lwd = 1, x = rep(bp[6, 1], 2), y = CItransitionsA[, 6], col = "blue")


#### transitions B####

transitionsB <- params[grep(x = names(params), pattern = "B_[A-Z]B")]
names(transitionsB) <- gsub(x = names(transitionsB), pattern = "q", replacement = "")
names(transitionsB) <- gsub(x = names(transitionsB), pattern = "B", replacement = "")
CItransitionsB <- ci[, grep(x = colnames(ci), pattern = "B_[B-Z]B")]
colnames(CItransitionsB) <- gsub(x = colnames(CItransitionsB), pattern = "B", replacement = "")
CItransitionsB <- CItransitionsB[c(1, 5), ]
bp <- barplot(transitionsB, ylim = c(0, 0.1), main = "Transition rates, regime B")
lines(lwd = 1, x = rep(bp[1, 1], 2), y = CItransitionsB[, 1], col = "blue")
lines(lwd = 1, x = rep(bp[2, 1], 2), y = CItransitionsB[, 2], col = "blue")
lines(lwd = 1, x = rep(bp[3, 1], 2), y = CItransitionsB[, 3], col = "blue")
lines(lwd = 1, x = rep(bp[4, 1], 2), y = CItransitionsB[, 4], col = "blue")
lines(lwd = 1, x = rep(bp[5, 1], 2), y = CItransitionsB[, 5], col = "blue")
lines(lwd = 1, x = rep(bp[6, 1], 2), y = CItransitionsB[, 6], col = "blue")


#### turnover A####
turnA <- params[grep(x = names(params), pattern = "turnover[A-Z]A")]
ciTurnA <- ci[c(1, 5), grep(x = colnames(ci), pattern = "turnover[A-Z]A")]
colnames(ciTurnA) <- gsub(x = colnames(ciTurnA), pattern = "A", replacement = "")

bp <- barplot(turnA, ylim = c(0, 0.1), main = "turnover, regime A")
lines(lwd = 1, x = rep(bp[1, 1], 2), y = ciTurnA[, 1], col = "blue")
lines(lwd = 1, x = rep(bp[2, 1], 2), y = ciTurnA[, 2], col = "blue")
lines(lwd = 1, x = rep(bp[3, 1], 2), y = ciTurnA[, 3], col = "blue")


#### turnover B####
turnB <- params[grep(x = names(params), pattern = "turnover[B-Z]B")]
ciTurnB <- ci[c(1, 5), grep(x = colnames(ci), pattern = "turnover[B-Z]B")]
colnames(ciTurnB) <- gsub(x = colnames(ciTurnB), pattern = "B", replacement = "")

bp <- barplot(turnB, ylim = c(0, 1), main = "turnover, regime B")
lines(lwd = 1, x = rep(bp[1, 1], 2), y = ciTurnB[, 1], col = "blue")
lines(lwd = 1, x = rep(bp[2, 1], 2), y = ciTurnB[, 2], col = "blue")
lines(lwd = 1, x = rep(bp[3, 1], 2), y = ciTurnB[, 3], col = "blue")


#### eps A####
epsA <- params[grep(x = names(params), pattern = "eps[A-Z]A")]
ciEpsA <- ci[c(1, 5), grep(x = colnames(ci), pattern = "eps[A-Z]A")]
colnames(ciEpsA) <- gsub(x = colnames(ciEpsA), pattern = "A", replacement = "")

bp <- barplot(epsA, ylim = c(0, 0.07), main = "eps, regime A")
lines(lwd = 1, x = rep(bp[1, 1], 2), y = ciEpsA[, 1], col = "blue")
lines(lwd = 1, x = rep(bp[2, 1], 2), y = ciEpsA[, 2], col = "blue")
lines(lwd = 1, x = rep(bp[3, 1], 2), y = ciEpsA[, 3], col = "blue")

#### eps B####
epsB <- params[grep(x = names(params), pattern = "eps[B-Z]B")]
ciEpsB <- ci[c(1, 5), grep(x = colnames(ci), pattern = "eps[B-Z]B")]
colnames(ciEpsB) <- gsub(x = colnames(ciEpsB), pattern = "B", replacement = "")

bp <- barplot(epsB, ylim = c(0, 0.6), main = "eps, regime B")
lines(lwd = 1, x = rep(bp[1, 1], 2), y = ciEpsB[, 1], col = "blue")
lines(lwd = 1, x = rep(bp[2, 1], 2), y = ciEpsB[, 2], col = "blue")
lines(lwd = 1, x = rep(bp[3, 1], 2), y = ciEpsB[, 3], col = "blue")


#### FIG 3#####
pdf("fig3.pdf")
layout(mat = matrix(nrow = 3, ncol = 2, dat = 1:6, byrow = T))

# transitiosn, regime A
bp <- barplot(transitionsA,
  ylim = c(0, 0.1), main = "Transition rates\nRegime A",
  names.arg = c("M->P", "M->L", "P->M", "P->L", "L->M", "L->P"), border = NA
)
lines(lwd = 1, x = rep(bp[1, 1], 2), y = CItransitionsA[, 1], col = "blue")
lines(lwd = 1, x = rep(bp[2, 1], 2), y = CItransitionsA[, 2], col = "blue")
lines(lwd = 1, x = rep(bp[3, 1], 2), y = CItransitionsA[, 3], col = "blue")
lines(lwd = 1, x = rep(bp[4, 1], 2), y = CItransitionsA[, 4], col = "blue")
lines(lwd = 1, x = rep(bp[5, 1], 2), y = CItransitionsA[, 5], col = "blue")
lines(lwd = 1, x = rep(bp[6, 1], 2), y = CItransitionsA[, 6], col = "blue")


# transitiosn, regime B
bp <- barplot(transitionsB,
  ylim = c(0, 0.1), main = "Transition rates\nRegime B",
  names.arg = c("M->P", "M->L", "P->M", "P->L", "L->M", "L->P"), border = NA
)
lines(lwd = 1, x = rep(bp[1, 1], 2), y = CItransitionsB[, 1], col = "blue")
lines(lwd = 1, x = rep(bp[2, 1], 2), y = CItransitionsB[, 2], col = "blue")
lines(lwd = 1, x = rep(bp[3, 1], 2), y = CItransitionsB[, 3], col = "blue")
lines(lwd = 1, x = rep(bp[4, 1], 2), y = CItransitionsB[, 4], col = "blue")
lines(lwd = 1, x = rep(bp[5, 1], 2), y = CItransitionsB[, 5], col = "blue")
lines(lwd = 1, x = rep(bp[6, 1], 2), y = CItransitionsB[, 6], col = "blue")

# turnover, regime A
bp <- barplot(turnA,
  ylim = c(0, 1), main = "Turnover (speciation+extinction)\nRegime A",
  names.arg = c("Monogamous", "Polygamous", "Lek"), border = NA
)
lines(lwd = 1, x = rep(bp[1, 1], 2), y = ciTurnA[, 1], col = "blue")
lines(lwd = 1, x = rep(bp[2, 1], 2), y = ciTurnA[, 2], col = "blue")
lines(lwd = 1, x = rep(bp[3, 1], 2), y = ciTurnA[, 3], col = "blue")

# turnover, B
bp <- barplot(turnB,
  ylim = c(0, 1), main = "Turnover (speciation+extinction)\nRegime B",
  names.arg = c("Monogamous", "Polygamous", "Lek"), border = NA
)
lines(lwd = 1, x = rep(bp[1, 1], 2), y = ciTurnB[, 1], col = "blue")
lines(lwd = 1, x = rep(bp[2, 1], 2), y = ciTurnB[, 2], col = "blue")
lines(lwd = 1, x = rep(bp[3, 1], 2), y = ciTurnB[, 3], col = "blue")


# Eps A
bp <- barplot(epsA,
  ylim = c(0, 0.6), main = "Extinction fraction (extinction/speciation)\nRegime A",
  names.arg = c("Monogamous", "Polygamous", "Lek"), border = NA
)
lines(lwd = 1, x = rep(bp[1, 1], 2), y = ciEpsA[, 1], col = "blue")
lines(lwd = 1, x = rep(bp[2, 1], 2), y = ciEpsA[, 2], col = "blue")
lines(lwd = 1, x = rep(bp[3, 1], 2), y = ciEpsA[, 3], col = "blue")

# Eps B
bp <- barplot(epsB,
  ylim = c(0, 0.6), main = "Extinction fraction (extinction/speciation)\nRegime B",
  names.arg = c("Monogamous", "Polygamous", "Lek"), border = NA
)
lines(lwd = 1, x = rep(bp[1, 1], 2), y = ciEpsB[, 1], col = "blue")
lines(lwd = 1, x = rep(bp[2, 1], 2), y = ciEpsB[, 2], col = "blue")
lines(lwd = 1, x = rep(bp[3, 1], 2), y = ciEpsB[, 3], col = "blue")


dev.off()
