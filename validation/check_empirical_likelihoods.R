# Run from the package root. Requires authorized HELIUS data as described in README.md.
source('R/revision_2026/empirical/00_context_ising_core.R')
source('R/revision_2026/empirical/01_data.R')
dat <- sf_load_data()
st <- ci_setup(sf_symptoms)
X <- as.matrix(dat[,sf_symptoms]); Tobs <- ci_features(X,st$pair_index)
cv <- data.frame(srl=dat$srl_z,srl2=dat$srl_z^2-mean(dat$srl_z^2), high=dat$high_srl,age=dat$age_z,female=dat$female,nondutch=dat$non_dutch)
f <- readRDS('res/revision_2026/empirical/helius/step1_model_fits/fits.rds')
cat('Sample',nrow(dat),'Low',sum(!dat$high_srl),'High',sum(dat$high_srl),'Missing ethnicity',sum(is.na(dat$non_dutch)),'\n')
out <- list()
for(tag in names(f)) for(m in names(f[[tag]])) {
 z <- f[[tag]][[m]]; cols <- setdiff(z$Z_cols,'(Intercept)'); rows <- if(tag=='adjusted') which(complete.cases(cv[,c('age','female','nondutch')])) else seq_len(nrow(cv))
 Z <- cbind(1,as.matrix(cv[rows,cols,drop=FALSE])); attr(z$S,'npar') <- z$npar
 v <- ci_negll_grad(z$par,z$S,Tobs[rows,,drop=FALSE],Z,st)
 out[[paste(tag,m)]] <- data.frame(adjustment=tag,model=m,n=length(rows),stored_logLik=z$logLik,recomputed_logLik=-as.numeric(v),abs_difference=abs(z$logLik+as.numeric(v)),max_gradient=max(abs(attr(v,'gradient'))),AIC=2*z$npar+2*as.numeric(v),BIC=log(length(rows))*z$npar+2*as.numeric(v))
 cat(tag,m,'LL difference',abs(z$logLik+as.numeric(v)),'max gradient',max(abs(attr(v,'gradient'))),'\n')
}
out <- do.call(rbind,out)
write.csv(out,'validation/empirical_likelihood_verification.csv',row.names=FALSE)
stopifnot(max(out$abs_difference)<1e-6,max(out$max_gradient)<.01)
cat('PASS all 18 fitted likelihoods and gradients verified against local HELIUS data\n')
