# Independent validation using R state enumeration and standard glm.
# Run from the repository root after Simulation 5.
source("R/revision_2026/00_parameters_uncentered01.R")
source("R/revision_2026/utils_uncentered01_model.R")
z <- readRDS("res/revision_2026/sim5/edge_estimates.rds")
metrics <- read.csv("res/revision_2026/sim5/replicate_metrics.csv")
states <- as.matrix(expand.grid(rep(list(0:1),N)))
ij <- z$ij
features <- states[,ij[,1]]*states[,ij[,2]]
field <- drop(states%*%gamma)

# Finite mixture identity, independently computed from all joint probabilities.
P <- c(-1,.2,1.3)
mass <- c(.2,.5,.3)
base <- drop(states%*%tau + features%*%omega[ij])
weights <- exp(outer(base,rep(1,3))+outer(field,P))
Z <- colSums(weights)
joint <- drop((weights/rep(Z,each=nrow(states)))%*%mass)
H <- function(t) log(sum(mass*exp(P*t)/Z))
max_identity_error <- 0
min_extra_association <- Inf
for(e in seq_len(nrow(ij))) {
  i <- ij[e,1]; j <- ij[e,2]
  zero <- which(states[,i]==0 & states[,j]==0)
  p10 <- zero+2^(i-1)
  p01 <- zero+2^(j-1)
  p11 <- zero+2^(i-1)+2^(j-1)
  actual <- log(joint[p11])+log(joint[zero])-log(joint[p10])-log(joint[p01])
  t <- field[zero]
  expected <- omega[i,j]+vapply(t,function(v)
    H(v+gamma[i]+gamma[j])-H(v+gamma[i])-H(v+gamma[j])+H(v),numeric(1))
  max_identity_error <- max(max_identity_error,max(abs(actual-expected)))
  min_extra_association <- min(min_extra_association,min(actual-omega[i,j]))
}
stopifnot(max_identity_error<1e-12,min_extra_association>=-1e-12)

# Replay a successful dataset from every scenario and sample size. Independent
# R sampling uses exactly the stored seed rule and the same uniforms.
design <- z$design
keys <- unique(design[,c("scenario","n")])
checks <- list()
RNGkind("Mersenne-Twister","Inversion","Rejection")
for(k in seq_len(nrow(keys))) {
  good <- metrics$scenario==keys$scenario[k] & metrics$n==keys$n[k] &
    metrics$method=="true_context" & metrics$ok
  chosen <- metrics[which(good)[1],]
  index <- which(design$scenario==chosen$scenario & design$n==chosen$n & design$rep==chosen$rep)
  par <- z$parameters[[chosen$scenario]]
  set.seed(chosen$seed)
  err <- 0
  for(g in 1:2) {
    context <- rnorm(chosen$n,par$mu[g],par$sigma[g])
    E <- rnorm(chosen$n)
    uniforms <- runif(chosen$n)
    b <- drop(states%*%tau+features%*%par$J[[g]][ij])
    Y <- t(vapply(seq_along(context),function(person) {
      eta <- b+context[person]*field
      probability <- exp(eta-max(eta))
      states[which(cumsum(probability)>=uniforms[person]*sum(probability))[1],]
    },numeric(N)))
    for(m in seq_along(z$methods)) {
      r <- c(NA,.2,.5,.8,1)[m]
      W <- if(is.na(r)) NULL else context+par$sigma[g]*sqrt((1-r)/r)*E
      fitted <- fit_edges(Y,W)[ij]
      err <- max(err,max(abs(fitted-z$estimates[[index]][g,m,])))
    }
  }
  stopifnot(err<1e-6)
  checks[[k]] <- data.frame(scenario=chosen$scenario,n=chosen$n,rep=chosen$rep,max_error=err)
}
write.csv(do.call(rbind,checks),"res/revision_2026/sim5/independent_replay.csv",row.names=FALSE)
writeLines(c(paste("Maximum mixture identity discrepancy",max_identity_error),
             paste("Minimum additional conditional association",min_extra_association),
             paste("Maximum replay discrepancy",max(vapply(checks,function(x)x$max_error,numeric(1))))),
           "res/revision_2026/sim5/independent_validation.txt")
cat("Independent mixture derivation and 12 dataset replays passed\n")
