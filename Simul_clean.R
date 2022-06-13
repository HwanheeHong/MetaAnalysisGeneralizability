# Simulation code for OM2, SM2 with 4 RCTs under the random sampling selection model
# To perform simulation for other outcome models, change values of phi1 and phi2
# To perform simulation for other selection models, change values of gam0, gam1 and gam2
library(dplyr)
library(rlang)
library(meta)
library(survey)
library(twang) ## For GBM
library(BayesTree) ## For BART
library(SuperLearner)
library(gam)
library(tmle) ## For tmle
library(lme4) ## For logistic regression with random effects
library(betareg)
source("Functions.R")

# True parameters
N_pop_true <- 11500
N_pop_obs <- 1500
x1_mean <- 0; x1_sd <- 1 
x2_p <- 0.5 
b0 <- 0; b1 <- 1; b2 <- 1; tau <- 2; phi1 <- 1; phi2 <- 0; y_sd<-1
TRUTH <- tau+phi2*0.5
N_RCT <- 5
gam0 <-c(logit(0.0019),logit(0.003),logit(0.003),logit(0.003),logit(0.003)); gam1 <- c(1,0,0,0,0); gam2 <- c(0,0,0,0,0)

# Simulation
nsim <- 1000
RES <- NULL
for(h in 1:nsim) {
  cat(h, "\n")
  set.seed(h)
  
  out <- fun.single.simulation(n=N_pop_true, x1_mean=x1_mean, x1_sd=x1_sd, x2_p=x2_p, b0=b0, b1=b1, b2=b2, tau=tau, phi1=phi1, phi2=phi2, y_sd=y_sd, TRUTH=TRUTH)
  RES <- rbind(RES, out) 
}

# Summary
OUT <- RES
vars <- c( paste0("naive",1:5), "naive.FMA", "naive.RMA", 
           paste0("logistic",1:5), "logistic1.FMA", "logistic1.RMA", 
           paste0("GBM",1:5), "GBM1.FMA", "GBM1.RMA",
           paste0("SL",1:5), "SL1.FMA", "SL1.RMA",
           paste0("BARTo",1:5), "BARTo.FMA", "BARTo.RMA",
           paste0("TMLEpop",1:5), "TMLEpop.FMA", "TMLEpop.RMA" )
stat.m <- apply(OUT,2,mean)
stat.v <- apply(OUT,2,var)

truth <- TRUTH
bias <- stat.m[ paste0(vars,".est") ] - truth 
mse <- bias^2 + stat.v[ paste0(vars,".est") ] 
cp <- stat.m[ paste0(vars,".cpind") ]

OUT_est <- OUT[,paste0(vars,".est")]
B <- 100
boot <- NULL
for(i in 1:B) {
  set.seed(i)
  sam <- sample(1:nsim, nsim, replace=T)
  temp <- OUT_est[sam,]   
  m <- apply(temp,2,mean)
  boot <- rbind(boot, m)
}
MCE.est <- apply(boot, 2, function(X) sum((X-mean(X))^2)/B )

smd.x1.m <- c(stat.m[paste0("naive",1:5,".SMD.x1")],NA,NA,stat.m[paste0("logistic",1:5,".SMD.x1")],NA,NA,stat.m[paste0("GBM",1:5,".SMD.x1")],NA,NA,stat.m[paste0("SL",1:5,".SMD.x1")],rep(NA,16)) 
smd.x1.sd <- c(sqrt(stat.v[paste0("naive",1:5,".SMD.x1")]),NA,NA,sqrt(stat.v[paste0("logistic",1:5,".SMD.x1")]),NA,NA,sqrt(stat.v[paste0("GBM",1:5,".SMD.x1")]),NA,NA,sqrt(stat.v[paste0("SL",1:5,".SMD.x1")]),rep(NA,16))
smd.x2.m <- c(stat.m[paste0("naive",1:5,".SMD.x2")],NA,NA,stat.m[paste0("logistic",1:5,".SMD.x2")],NA,NA,stat.m[paste0("GBM",1:5,".SMD.x2")],NA,NA,stat.m[paste0("SL",1:5,".SMD.x2")],rep(NA,16)) 
smd.x2.sd <- c(sqrt(stat.v[paste0("naive",1:5,".SMD.x2")]),NA,NA,sqrt(stat.v[paste0("logistic",1:5,".SMD.x2")]),NA,NA,sqrt(stat.v[paste0("GBM",1:5,".SMD.x2")]),NA,NA,sqrt(stat.v[paste0("SL",1:5,".SMD.x2")]),rep(NA,16))
wgt.mean.m <- c(rep(NA,7),stat.m[paste0("logistic",1:5,".wgt.mean")],NA,NA,stat.m[paste0("GBM",1:5,".wgt.mean")],NA,NA,stat.m[paste0("SL",1:5,".wgt.mean")],rep(NA,16))  
wgt.sd.m <- c(rep(NA,7),stat.m[paste0("logistic",1:5,".wgt.sd")],NA,NA,stat.m[paste0("GBM",1:5,".wgt.sd")],NA,NA,stat.m[paste0("SL",1:5,".wgt.sd")],rep(NA,16)) 
wgt.min.m <- c(rep(NA,7),stat.m[paste0("logistic",1:5,".wgt.min")],NA,NA,stat.m[paste0("GBM",1:5,".wgt.min")],NA,NA,stat.m[paste0("SL",1:5,".wgt.min")],rep(NA,16)) 
wgt.max.m <- c(rep(NA,7),stat.m[paste0("logistic",1:5,".wgt.max")],NA,NA,stat.m[paste0("GBM",1:5,".wgt.max")],NA,NA,stat.m[paste0("SL",1:5,".wgt.max")],rep(NA,16)) 
wgt.low.m <- c(rep(NA,7),stat.m[paste0("logistic",1:5,".wgt.low")],NA,NA,stat.m[paste0("GBM",1:5,".wgt.low")],NA,NA,stat.m[paste0("SL",1:5,".wgt.low")],rep(NA,16)) 
wgt.up.m <- c(rep(NA,7),stat.m[paste0("logistic",1:5,".wgt.up")],NA,NA,stat.m[paste0("GBM",1:5,".wgt.up")],NA,NA,stat.m[paste0("SL",1:5,".wgt.up")],rep(NA,16)) 
wgt.q1.m <- c(rep(NA,7),stat.m[paste0("logistic",1:5,".wgt.q1")],NA,NA,stat.m[paste0("GBM",1:5,".wgt.q1")],NA,NA,stat.m[paste0("SL",1:5,".wgt.q1")],rep(NA,16)) 
wgt.q2.m <- c(rep(NA,7),stat.m[paste0("logistic",1:5,".wgt.q2")],NA,NA,stat.m[paste0("GBM",1:5,".wgt.q2")],NA,NA,stat.m[paste0("SL",1:5,".wgt.q2")],rep(NA,16)) 
wgt.q3.m <- c(rep(NA,7),stat.m[paste0("logistic",1:5,".wgt.q3")],NA,NA,stat.m[paste0("GBM",1:5,".wgt.q3")],NA,NA,stat.m[paste0("SL",1:5,".wgt.q3")],rep(NA,16)) 
wgt.mean.sd <- c(rep(NA,7),sqrt(stat.v[paste0("logistic",1:5,".wgt.mean")]),NA,NA,sqrt(stat.v[paste0("GBM",1:5,".wgt.mean")]),NA,NA,sqrt(stat.v[paste0("SL",1:5,".wgt.mean")]),rep(NA,16))
wgt.sd.sd <- c(rep(NA,7),sqrt(stat.v[paste0("logistic",1:5,".wgt.sd")]),NA,NA,sqrt(stat.v[paste0("GBM",1:5,".wgt.sd")]),NA,NA,sqrt(stat.v[paste0("SL",1:5,".wgt.sd")]),rep(NA,16))
wgt.min.sd <- c(rep(NA,7),sqrt(stat.v[paste0("logistic",1:5,".wgt.min")]),NA,NA,sqrt(stat.v[paste0("GBM",1:5,".wgt.min")]),NA,NA,sqrt(stat.v[paste0("SL",1:5,".wgt.min")]),rep(NA,16))
wgt.max.sd <- c(rep(NA,7),sqrt(stat.v[paste0("logistic",1:5,".wgt.max")]),NA,NA,sqrt(stat.v[paste0("GBM",1:5,".wgt.max")]),NA,NA,sqrt(stat.v[paste0("SL",1:5,".wgt.max")]),rep(NA,16))
wgt.low.sd <- c(rep(NA,7),sqrt(stat.v[paste0("logistic",1:5,".wgt.low")]),NA,NA,sqrt(stat.v[paste0("GBM",1:5,".wgt.low")]),NA,NA,sqrt(stat.v[paste0("SL",1:5,".wgt.low")]),rep(NA,16))
wgt.up.sd <- c(rep(NA,7),sqrt(stat.v[paste0("logistic",1:5,".wgt.up")]),NA,NA,sqrt(stat.v[paste0("GBM",1:5,".wgt.up")]),NA,NA,sqrt(stat.v[paste0("SL",1:5,".wgt.up")]),rep(NA,16))
wgt.q1.sd <- c(rep(NA,7),sqrt(stat.v[paste0("logistic",1:5,".wgt.q1")]),NA,NA,sqrt(stat.v[paste0("GBM",1:5,".wgt.q1")]),NA,NA,sqrt(stat.v[paste0("SL",1:5,".wgt.q1")]),rep(NA,16))
wgt.q2.sd <- c(rep(NA,7),sqrt(stat.v[paste0("logistic",1:5,".wgt.q2")]),NA,NA,sqrt(stat.v[paste0("GBM",1:5,".wgt.q2")]),NA,NA,sqrt(stat.v[paste0("SL",1:5,".wgt.q2")]),rep(NA,16))
wgt.q3.sd <- c(rep(NA,7),sqrt(stat.v[paste0("logistic",1:5,".wgt.q3")]),NA,NA,sqrt(stat.v[paste0("GBM",1:5,".wgt.q3")]),NA,NA,sqrt(stat.v[paste0("SL",1:5,".wgt.q3")]),rep(NA,16))
STAT <- rbind(bias, mse, cp, MCE.est,
              smd.x1.m, smd.x1.sd, smd.x2.m, smd.x2.sd, 
              wgt.mean.m, wgt.sd.m, wgt.min.m, wgt.max.m, wgt.low.m, wgt.up.m, wgt.q1.m, wgt.q2.m, wgt.q3.m,
              wgt.mean.sd, wgt.sd.sd, wgt.min.sd, wgt.max.sd, wgt.low.sd, wgt.up.sd, wgt.q1.sd, wgt.q2.sd, wgt.q3.sd)
write.table(STAT, paste0("results//",simname,"_stat.txt"), sep="\t")

BASIC.m <- stat.m[ c("n.pop", "mean.pop.x1", "sd.pop.x1", "min.pop.x1", "max.pop.x1", "low.pop.x1", "up.pop.x1", "q1.pop.x1", "q2.pop.x1","q3.pop.x1", "mean.pop.x2") ]
BASIC.m <- c(BASIC.m, stat.m[ paste( c("n","mean","sd","min","max","low","up","q1","q2","q3","mean"),".rct",rep(1:5, each="11"),c("",rep(".x1",9),".x2"),sep="")] )
BASIC.sd <- sqrt(stat.v[ c("n.pop", "mean.pop.x1", "sd.pop.x1", "min.pop.x1", "max.pop.x1", "low.pop.x1", "up.pop.x1", "q1.pop.x1", "q2.pop.x1","q3.pop.x1", "mean.pop.x2") ])
BASIC.sd <- c(BASIC.sd, sqrt(stat.v[ paste( c("n","mean","sd","min","max","low","up","q1","q2","q3","mean"),".rct",rep(1:5, each="11"),c("",rep(".x1",9),".x2"),sep="")]))
BASIC_STAT <- rbind(BASIC.m, BASIC.sd)
write.table(BASIC_STAT, paste0("results//",simname,"_basicstat.txt"), sep="\t")
