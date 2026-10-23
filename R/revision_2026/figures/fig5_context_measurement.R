suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})
source("R/revision_2026/figures/theme_publication.R")
root <- "res/revision_2026/sim5"
dat <- read.csv(file.path(root,"replicate_metrics.csv"))
stopifnot(nrow(dat)==30000L)
failures <- dat |> group_by(scenario,n,method) |>
  summarise(attempted=n(),failed=sum(!ok),.groups="drop")
write.csv(failures,file.path(root,"failures.csv"),row.names=FALSE)
# Paired comparisons require the same valid datasets for every method.
dat <- dat |> group_by(scenario,n,rep) |> filter(all(ok)) |> ungroup()
stopifnot(all(is.finite(dat$contrast_rmse)))
metrics <- c("contrast_rmse","signed_sum","signed_sum_error","unchanged_mae",
             "changed_oriented","changed_mae","group_rmse")
se <- function(x) if (all(is.na(x))) NA_real_ else sd(x,na.rm=TRUE)/sqrt(sum(!is.na(x)))
avg <- function(x) if (all(is.na(x))) NA_real_ else mean(x,na.rm=TRUE)
summary <- dat |>
  group_by(scenario,n,method) |>
  summarise(repetitions=n(),across(all_of(metrics),list(mean=avg,mcse=se)),.groups="drop")
write.csv(summary,file.path(root,"summary.csv"),row.names=FALSE)
paired <- dat |>
  select(scenario,n,rep,method,all_of(metrics)) |>
  pivot_wider(names_from=method,values_from=all_of(metrics))
pairs <- lapply(setdiff(unique(dat$method),"omit"),function(m) {
  ans <- paired |> select(scenario,n,rep)
  for (metric in metrics) ans[[metric]] <- paired[[paste0(metric,"_",m)]]-paired[[paste0(metric,"_omit")]]
  ans$method <- m
  ans
}) |> bind_rows() |>
  group_by(scenario,n,method) |>
  summarise(across(all_of(metrics),list(mean=avg,mcse=se)),.groups="drop")
write.csv(pairs,file.path(root,"paired_vs_omit.csv"),row.names=FALSE)

# Check RMSE and signed sums against the saved edge estimates.
saved <- readRDS(file.path(root,"edge_estimates.rds"))
edge_rows <- vector("list",nrow(saved$design))
for (k in seq_len(nrow(saved$design))) {
  cell <- saved$design[k,]
  est <- saved$estimates[[k]]
  delta <- est[2,,]-est[1,,]
  truth <- saved$parameters[[cell$scenario]]$delta
  check <- dat[dat$scenario==cell$scenario & dat$n==cell$n & dat$rep==cell$rep,]
  if(!nrow(check)) next
  check <- check[match(saved$methods,check$method),]
  stopifnot(max(abs(sqrt(rowMeans(sweep(delta,2,truth)^2))-check$contrast_rmse))<1e-12)
  stopifnot(max(abs(rowSums(delta)-check$signed_sum))<1e-12)
  edge_rows[[k]] <- data.frame(scenario=cell$scenario,n=cell$n,rep=cell$rep,
    method=rep(saved$methods,each=length(truth)),edge=rep(saved$edges,length(saved$methods)),
    truth=rep(truth,length(saved$methods)),estimate=as.vector(t(delta)))
}
edge_summary <- bind_rows(edge_rows) |>
  group_by(scenario,n,method,edge,truth) |>
  summarise(mean=mean(estimate),bias=mean(estimate-truth),
    bias_mcse=sd(estimate)/sqrt(n()),rmse=sqrt(mean((estimate-truth)^2)),.groups="drop")
write.csv(edge_summary,file.path(root,"edge_summary.csv"),row.names=FALSE)
writeLines(paste("All",nrow(dat),"valid replicate RMSE and signed-sum values agree with independent reconstruction to 1e-12."),
           file.path(root,"metric_validation.txt"))

levels_method <- c("omit","reliability_02","reliability_05","reliability_08","true_context")
scenario_labels <- c(neither="Neither difference",context_only="Context only",
  coupling_only="Coupling only",both="Context and coupling",
  balanced="Opposing coupling changes",context_mean_only="Context mean only")
summary <- summary |>
  mutate(method=factor(method,levels=levels_method),
         size=factor(n,levels=c(500,2000),labels=c("500 per group","2000 per group")),
         scenario=factor(scenario,levels=names(scenario_labels),labels=scenario_labels))
base_plot <- function(d, metric, ylabel) {
  ggplot(d,aes(x=method,y=.data[[paste0(metric,"_mean")]],colour=size,shape=size,group=size))+
    geom_line(linewidth=.6)+geom_point(size=2)+
    geom_errorbar(aes(ymin=.data[[paste0(metric,"_mean")]]-1.96*.data[[paste0(metric,"_mcse")]],
                      ymax=.data[[paste0(metric,"_mean")]]+1.96*.data[[paste0(metric,"_mcse")]]),
                  width=.12,linewidth=.4)+
    scale_x_discrete(labels=c("Omitted","0.2","0.5","0.8","1.0"))+
    scale_colour_manual(values=c("#7A8DA6","#1B9E77"))+
    scale_shape_manual(values=c(17,16))+
    labs(x="Within-group reliability of the context measure",y=ylabel,colour=NULL,shape=NULL)+
    theme_pub()+theme(strip.text=element_text(size=10),legend.position="bottom")
}
main <- base_plot(filter(summary,scenario %in% scenario_labels[1:4]),
                  "contrast_rmse","Mean root mean squared error")+
  facet_wrap(~scenario,ncol=2)+
  scale_y_continuous(limits=c(0,NA))
dir.create("figs/revision_2026",showWarnings=FALSE,recursive=TRUE)
ggsave("figs/revision_2026/Figure_sim5_recovery.pdf",main,width=7.2,height=5.4)
ggsave("figs/revision_2026/Figure_sim5_recovery.png",main,width=7.2,height=5.4,dpi=200)

targets <- data.frame(scenario=factor(scenario_labels,levels=scenario_labels),
                     target=c(0,0,.9,.9,0,0))
signed <- base_plot(summary,"signed_sum","Mean signed sum of estimated differences")+
  geom_hline(data=targets,aes(yintercept=target),linetype="dashed",colour="grey35")+
  facet_wrap(~scenario,ncol=2,scales="free_y")
ggsave("figs/revision_2026/Figure_sim5_signed.pdf",signed,width=7.2,height=7.1)
ggsave("figs/revision_2026/Figure_sim5_signed.png",signed,width=7.2,height=7.1,dpi=180)
changed <- base_plot(filter(summary,!is.na(changed_oriented_mean)),
                     "changed_oriented","Mean recovery of a true coupling change")+
  geom_hline(yintercept=.3,linetype="dashed",colour="grey35")+
  facet_wrap(~scenario,ncol=3)
ggsave("figs/revision_2026/Figure_sim5_changed.pdf",changed,width=7.2,height=3.3)
ggsave("figs/revision_2026/Figure_sim5_changed.png",changed,width=7.2,height=3.3,dpi=180)
print(summary |> filter(size=="2000 per group") |>
        select(scenario,method,contrast_rmse_mean,signed_sum_mean,changed_oriented_mean),n=30)
