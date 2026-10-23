# Simulation 5. Run from the repository root.
# Arguments are repetitions, workers, and optional R backend.
# Manuscript settings are 500 repetitions and four workers.
suppressPackageStartupMessages({
  library(Rcpp)
  library(parallel)
  library(lpSolve)
})
source("R/revision_2026/00_parameters_uncentered01.R")
source("R/revision_2026/utils_uncentered01_model.R")
args <- commandArgs(trailingOnly=TRUE)
use_cpp <- !(length(args) > 2 && args[3] == "R")
if(use_cpp) use_cpp <- tryCatch({
  sourceCpp("R/revision_2026/sim5_core.cpp"); TRUE
},error=function(e) {
  message("Compiled backend unavailable. Using standard R. ",conditionMessage(e)); FALSE
})
if(!use_cpp) {
  sim5_sample <- function(states,base,field,context,uniforms) {
    t(vapply(seq_along(context),function(k) {
      eta <- base + context[k]*field
      prob <- exp(eta-max(eta))
      states[which(cumsum(prob)>=uniforms[k]*sum(prob))[1],]
    },numeric(ncol(states))))
  }
  sim5_fit_cpp <- function(Y,context,adjust) {
    p <- ncol(Y)
    beta <- matrix(0,p,p)
    fail <- function(reason) list(J=matrix(NA_real_,p,p),ok=FALSE,max_score=NA_real_,
                                  iterations=NA_integer_,max_coefficient=NA_real_,reason=reason)
    for(i in seq_len(p-1)) for(j in (i+1):p) {
      if(any(tabulate(1L+as.integer(2*Y[,i]+Y[,j]),nbins=4L)==0L))
        return(fail("zero symptom-pair cell"))
    }
    worst <- largest <- iterations <- 0
    for(i in seq_len(p)) {
      others <- setdiff(seq_len(p),i)
      X <- cbind(1,Y[,others,drop=FALSE])
      if(adjust) X <- cbind(X,context)
      fit <- suppressWarnings(glm.fit(X,Y[,i],family=binomial(),
                                      control=glm.control(epsilon=1e-12,maxit=100)))
      if(!fit$converged || any(!is.finite(fit$coefficients))) return(fail("numerical failure"))
      score <- max(abs(crossprod(X,Y[,i]-fit$fitted.values)))
      if(score>=1e-7) return(fail("numerical failure"))
      beta[i,others] <- fit$coefficients[2:p]
      worst <- max(worst,score)
      largest <- max(largest,abs(fit$coefficients))
      iterations <- max(iterations,fit$iter)
    }
    list(J=(beta+t(beta))/2,ok=TRUE,max_score=worst,
         iterations=iterations,max_coefficient=largest,reason="success")
  }
}
has_separation <- function(X,y) {
  # A nonzero nonnegative vector of signed margins gives a separating
  # direction. Split free coefficients into bounded positive/negative parts.
  D <- X * (2*y-1)
  A <- cbind(D,-D)
  objective <- colSums(A)
  fit <- lpSolve::lp("max",objective,rbind(A,diag(ncol(A))),
    c(rep(">=",nrow(A)),rep("<=",ncol(A))),c(rep(0,nrow(A)),rep(1,ncol(A))))
  if(fit$status != 0) stop("Separation linear program failed")
  fit$objval > 1e-7
}
sim5_fit <- function(Y,W,adjust) {
  ans <- sim5_fit_cpp(Y,W,adjust)
  # Check for separation when coefficients exceed 8.
  # Large finite estimates are retained if no separating direction is found.
  if (ans$ok && ans$max_coefficient > 8) {
    for (i in seq_len(ncol(Y))) {
      X <- cbind(1,Y[,-i,drop=FALSE])
      if (adjust) X <- cbind(X,W)
      if(has_separation(X,Y[,i])) {
        ans$ok <- FALSE; ans$J[] <- NA_real_
        ans$reason <- "linear separation"; break
      }
    }
  }
  ans
}
B <- if (length(args)) as.integer(args[1]) else 500L
workers <- if (length(args) > 1) as.integer(args[2]) else 4L
stopifnot(B > 0, B <= 500, workers > 0)
out <- if (B == 500) "res/revision_2026/sim5" else "res/revision_2026/sim5_pilot"
dir.create(out, recursive=TRUE, showWarnings=FALSE)
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
states <- as.matrix(expand.grid(rep(list(0:1), N)))
ij <- which(upper.tri(omega), arr.ind=TRUE)
edge_names <- paste(symptoms[ij[,1]], symptoms[ij[,2]], sep="__")
features <- states[,ij[,1]] * states[,ij[,2]]
context_field <- drop(states %*% gamma)
scenarios <- c("neither", "context_only", "coupling_only", "both",
               "balanced", "context_mean_only")
methods <- c("omit", "reliability_02", "reliability_05", "reliability_08", "true_context")
reliabilities <- c(NA, .2, .5, .8, 1)

parameters <- function(scenario) {
  J1 <- omega
  if (scenario %in% c("coupling_only", "both")) {
    for (pair in list(c(1,2), c(3,4), c(4,7))) {
      J1[pair[1],pair[2]] <- J1[pair[2],pair[1]] <- omega[pair[1],pair[2]] + .30
    }
  }
  if (scenario == "balanced") {
    J1[1,2] <- J1[2,1] <- omega[1,2] + .30
    J1[3,4] <- J1[4,3] <- omega[3,4] - .30
  }
  if (scenario %in% c("context_only", "both", "balanced")) {
    mu <- c(-.5,.5); sigma <- c(.35,1)
  } else if (scenario == "context_mean_only") {
    mu <- c(-.5,.5); sigma <- c(.6,.6)
  } else { mu <- c(0,0); sigma <- c(.6,.6) }
  list(J=list(omega,J1), mu=mu, sigma=sigma,
       delta=(J1-omega)[ij])
}

# Check conditional log odds directly from enumerated joint probabilities.
pr <- parameters("both")
base <- drop(states %*% tau + features %*% pr$J[[2]][ij])
conditional_error <- 0
for (i in seq_len(N)) {
  zero <- which(states[,i] == 0)
  one <- zero + 2^(i-1)
  actual <- base[one] - base[zero] + .7*(context_field[one] - context_field[zero])
  expected <- tau[i] + drop(states[zero,,drop=FALSE] %*% pr$J[[2]][i,]) + .7*gamma[i]
  conditional_error <- max(conditional_error,max(abs(actual-expected)))
}
stopifnot(conditional_error < 1e-12)
set.seed(73401)
testP <- rnorm(2000)
testY <- sim5_sample(states, base, context_field, testP, runif(2000))
test_errors <- vapply(c(FALSE,TRUE), function(adjust) {
  cpp <- sim5_fit(testY,testP,adjust)
  reference <- fit_edges(testY,if(adjust) testP else NULL)
  stopifnot(cpp$ok)
  max(abs(cpp$J-reference))
},numeric(1))
stopifnot(max(test_errors) < 1e-6)
stopifnot(has_separation(cbind(1,c(-2,-1,1,2)),c(0,0,1,1)))
stopifnot(!has_separation(cbind(1,c(-2,-1,1,2)),c(0,1,0,1)))
separatedY <- testY
separatedY[separatedY[,1]==1,2] <- 0
stopifnot(!sim5_fit(separatedY,testP,FALSE)$ok)
# Independent multinomial frequency check at a fixed contextual value.
eta <- base + .2*context_field
prob <- exp(eta-max(eta)); prob <- prob/sum(prob)
testn <- 100000L
testS <- sim5_sample(states,base,context_field,rep(.2,testn),runif(testn))
moments <- colMeans(testS)
target <- drop(crossprod(prob,states))
standardized <- abs(moments-target)/sqrt(target*(1-target)/testn)
stopifnot(max(standardized)<6)
writeLines(c(
  paste("Conditional log-odds maximum error",format(conditional_error,digits=12)),
  paste("Nodewise fitter versus R glm maximum errors",paste(format(test_errors,digits=12),collapse=" ")),
  paste("Sampler maximum marginal z score",max(standardized)),
  "The sampler check uses 100000 independent profiles at a fixed context."
),file.path(out,"validation.txt"))

design <- expand.grid(scenario=scenarios,n=c(500L,2000L),
                      rep=seq_len(B),stringsAsFactors=FALSE)
run_one <- function(k) {
  cell <- design[k,]
  sc <- match(cell$scenario,scenarios)
  seed <- 20261020L + sc*100000L + match(cell$n,c(500L,2000L))*10000L + cell$rep
  set.seed(seed)
  par <- parameters(cell$scenario)
  estimates <- array(NA_real_,c(2,length(methods),nrow(ij)))
  okay <- matrix(FALSE,2,length(methods))
  scores <- matrix(NA_real_,2,length(methods))
  reasons <- matrix("",2,length(methods))
  for (g in 1:2) {
    P <- rnorm(cell$n,par$mu[g],par$sigma[g])
    E <- rnorm(cell$n)
    base_g <- drop(states %*% tau + features %*% par$J[[g]][ij])
    Y <- sim5_sample(states,base_g,context_field,P,runif(cell$n))
    for (m in seq_along(methods)) {
      r <- reliabilities[m]
      W <- if (is.na(r)) P else P + par$sigma[g]*sqrt((1-r)/r)*E
      fit <- sim5_fit(Y,W,!is.na(r))
      okay[g,m] <- fit$ok
      scores[g,m] <- fit$max_score
      reasons[g,m] <- fit$reason
      estimates[g,m,] <- fit$J[ij]
    }
  }
  delta <- estimates[2,,] - estimates[1,,]
  change <- abs(par$delta)>1e-12
  tab <- do.call(rbind,lapply(seq_along(methods),function(m) {
    err <- delta[m,]-par$delta
    data.frame(scenario=cell$scenario,n=cell$n,rep=cell$rep,seed=seed,
      method=methods[m],ok=all(okay[,m]),max_score=if(all(okay[,m])) max(scores[,m]) else NA_real_,
      failure=if(all(okay[,m])) "" else paste(unique(reasons[!okay[,m],m]),collapse=" | "),
      contrast_rmse=sqrt(mean(err^2)), signed_sum=sum(delta[m,]),
      signed_sum_error=sum(err), unchanged_mae=mean(abs(err[!change])),
      changed_oriented=if(any(change)) mean(sign(par$delta[change])*delta[m,change]) else NA_real_,
      changed_mae=if(any(change)) mean(abs(err[change])) else NA_real_,
      group_rmse=sqrt(mean(c(estimates[1,m,]-par$J[[1]][ij],
                             estimates[2,m,]-par$J[[2]][ij])^2)))
  }))
  list(metrics=tab,estimates=estimates)
}

started <- Sys.time()
cat("Design",nrow(design),"datasets and",length(methods),"methods\n")
all_results <- vector("list",nrow(design))
# Process 120 datasets at a time and save after the final batch.
chunks <- split(seq_len(nrow(design)),ceiling(seq_len(nrow(design))/120))
for (a in seq_along(chunks)) {
  ids <- chunks[[a]]
  chunk <- mclapply(ids,run_one,mc.cores=workers,mc.set.seed=FALSE)
  if (any(vapply(chunk,inherits,logical(1),"try-error"))) stop("Worker failed")
  all_results[ids] <- chunk
  cat(format(Sys.time()),"completed",max(ids),"of",nrow(design),"\n")
  flush.console()
}
metrics <- do.call(rbind,lapply(all_results,`[[`,"metrics"))
write.csv(metrics,file.path(out,"replicate_metrics.csv"),row.names=FALSE)
saveRDS(list(design=design,methods=methods,edges=edge_names,
             ij=ij,parameters=lapply(setNames(scenarios,scenarios),parameters),
             estimates=lapply(all_results,`[[`,"estimates"),
             seed_rule="20261020 + scenario_index*100000 + n_index*10000 + replicate"),
        file.path(out,"edge_estimates.rds"),compress="xz")
writeLines(capture.output(sessionInfo()),file.path(out,"sessionInfo.txt"))
writeLines(c(paste("Started",started),paste("Finished",Sys.time()),
             paste("Backend",if(use_cpp) "RcppArmadillo" else "standard R"),
             paste("Failed method fits",sum(!metrics$ok)),
             paste("Maximum successful score",max(metrics$max_score[metrics$ok])),
             paste("Repetitions per scenario and sample size",B)),
           file.path(out,"run_record.txt"))
cat("Finished. Failed method fits",sum(!metrics$ok),"\n")
