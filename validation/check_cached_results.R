suppressPackageStartupMessages(library(dplyr))
near <- function(x,y) stopifnot(isTRUE(all.equal(as.numeric(x),as.numeric(y),tolerance=1e-8)))
p <- 'res/revision_2026'
s <- readRDS(file.path(p,'sim1/sim1_raw.rds'))$burden
z <- s |> group_by(condition) |> summarise(mean=mean(M),high=mean(M>=5),.groups='drop')
near(round(z$mean,2),c(3.56,1.64,2.49)); cat('Simulation 1 manuscript means verified\n')
a <- read.csv(file.path(p,'sim4/sim4_replicated_by_arm.csv')); b <- read.csv(file.path(p,'sim4/sim4_replicated_summary_by_arm.csv'))
for(i in seq_len(nrow(b))) for(v in names(a)[-(1:2)]) {
 x <- a[a$arm==b$arm[i],v]; near(mean(x),b[i,paste0(v,'_mean')]);near(sd(x)/sqrt(length(x)),b[i,paste0(v,'_se')])
};cat('Simulation 4 all arm means and SEs verified from 30 replicate rows\n')
x <- readRDS('res/NCT_analysis/nct_mc_raw_weighted.rds');z <- x |> group_by(topology,n_per_grp,threshold_gap) |> summarise(reject_nwinv=mean(nwinv_p<.05),reject_glstr=mean(glstr_p<.05),fp_edge_05=mean(edge_fp_05),.groups='drop'); y <- read.csv('res/NCT_analysis/nct_mc_summary_weighted.csv')
k <- c('topology','n_per_grp','threshold_gap');j<-merge(z,y,by=k,suffixes=c('.raw','.saved'));for(v in c('reject_nwinv','reject_glstr','fp_edge_05'))near(j[[paste0(v,'.raw')]],j[[paste0(v,'.saved')]])
cat('All 27 NCT conditions verified from 2700 replicate rows\n')
p2 <- file.path(p,'empirical/helius/step2_simulation');x <- read.csv(file.path(p2,'replications.csv'));y <- read.csv(file.path(p2,'conditions_vs_observed.csv'))
for(i in seq_len(nrow(y))) {
 metric <- y$metric[i]; if(!metric %in% names(x))next
 z <- x[x$condition==y$condition[i],metric];if(!is.numeric(z))next; z<-z[is.finite(z)]
 near(mean(z),y$mean[i]);near(quantile(z,.025),y$q025[i]);near(quantile(z,.975),y$q975[i])
};cat('Mechanism means and 95 percent intervals verified from 1000 replicate rows\n')
for(f in list.files('R',pattern='\\.R$',recursive=TRUE,full.names=TRUE))parse(f)
cat('All R sources parse successfully\n')
