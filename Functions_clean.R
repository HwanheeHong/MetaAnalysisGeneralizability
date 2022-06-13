######
###### Basic functions
######
expit <- function(p) { exp(p)/(1+exp(p)) }
logit <- function(p) { log(p/(1-p)) }
q.025<-function(x){quantile(x, probs = 0.025)}
q.975<-function(x){quantile(x, probs = 0.975)}
q.25 <- function(x){quantile(x, probs = 0.25)}
q.75 <- function(x){quantile(x, probs = 0.75)}

## data.frame to list by group
split_tibble <- function(tibble, col = 'col') tibble %>% split(., .[,col])

######
###### Generate population data
######
## Generate true target population
fun.gen.true.pop <- function(n, x1_mean, x1_sd, x2_p, b0, b1, b2, tau, phi1, phi2, y_sd) {
  X1 <- rnorm(n, x1_mean, x1_sd)
  X2 <- rbinom(n,1,x2_p)
  Y0 <- rnorm(n, b0 + b1*X1 + b2*X2, y_sd)
  Y1 <- rnorm(n, b0 + b1*X1 + b2*X2 + tau + phi1*X1 + phi2*X2, y_sd)
  POP_true <- data.frame(X1=X1, X2=X2, Y0=Y0, Y1=Y1)
  return(POP_true)
}

## Sample population data from the true population
fun.sam.pop <- function(data, n_pop) {
  n <- nrow(data)
  Index <- 1:n
  temp <- sample(Index, n_pop, replace=FALSE)
  Ind_pop <- rep(0,n)
  Ind_pop[temp] <- 1 
  
  data_pop <- data[Ind_pop==1,]
  data_sam <- data[Ind_pop==0,]
  
  data_list <- list(data_pop, data_sam)
  return(data_list)
}

######
###### Sample RCTs
######
fun.sam.rct.single <- function(data, gam_single, study_ID) {
  prob <- expit(cbind(rep(1,nrow(data)),data$X1,data$X2) %*% gam_single)
  ind <- rbinom(nrow(data),1,prob)
  RCT <- data[ind==1, ]
  RCT$trt <- rbinom(sum(ind),1,0.5) 
  RCT$study <- study_ID
  
  return(RCT)
}

fun.sam.rct.all <- function(data, gam0, gam1, gam2) {
  gam <- rbind(gam0, gam1, gam2)
  RCT <- rbind(fun.sam.rct.single(data, gam[,1], 1),
               fun.sam.rct.single(data, gam[,2], 2),
               fun.sam.rct.single(data, gam[,3], 3),
               fun.sam.rct.single(data, gam[,4], 4),
               fun.sam.rct.single(data, gam[,5], 5))
  return(RCT)
}

######
###### Estimation
######
## Estimate unweighted SATE
fun.naive.ate <- function(data, TRUTH, model_name) {
  mod_unw <- glm(Y~trt, data=data, family = gaussian)
  ci_unw <- confint(mod_unw)[2,]
  res_unw <- c(summary(mod_unw)$coef[2,][1:2], ci_unw)
  res_unw <- c(res_unw, ifelse( ci_unw[1] <TRUTH & TRUTH < ci_unw[2], 1, 0))
  
  names(res_unw) <- paste0(model_name,".",c("est","se","low","up","cpind"))   
  return(res_unw)
}

## Estimate weighted SATE
fun.weighted.ate <- function(data, TRUTH, model_name) {
  design.ps <- svydesign(ids=~1, weights=~wgt, data=data)
  mod_w <- svyglm(Y ~ trt, design=design.ps)
  ci_w <- confint(mod_w)[2,]
  res_w <- c(summary(mod_w)$coef[2,][1:2], ci_w)
  res_w <- c(res_w, ifelse( ci_w[1]<TRUTH & TRUTH<ci_w[2],1,0 ))
  
  names(res_w) <- paste0(model_name,".",c("est","se","low","up","cpind"))
  return(res_w)
}

fun.weighted.trial.ate <- function(DATA,trial) {
  data <- DATA[DATA$study==trial,]
  design.ps <- svydesign(ids=~1, weights=~wgt, data=data)
  mod_w <- svyglm(Y ~ trt, design=design.ps)
  ci_w <- confint(mod_w)[2,]
  res_w <- c(as.numeric(trial),summary(mod_w)$coef[2,], ci_w)
  names(res_w) <- c('study','SATE','SE','t','pval','low','up')
  return(res_w)
}

fun.weighted.ate2 <- function(data) {
  out <- rbind(fun.weighted.trial.ate(data,'1'),
               fun.weighted.trial.ate(data,'2'),
               fun.weighted.trial.ate(data,'3'),
               fun.weighted.trial.ate(data,'4'),
               fun.weighted.trial.ate(data,'5'))
  return(out)
}

## Calculate weights, logistic regression
fun.logistic.wgt <- function(data) {
  fit <- glm(index ~ X1 + X2, family="binomial", data=data)
  ps <- predict(fit, type="response")
  wgt <- ps/(1-ps)
  
  return(wgt)
}

## Calculate weights, GBM
#### Use twang package.
#### Weights for ATT are 1 for the treatment cases (population) and p/(1-p) for the control cases (sample)
#### So, I use estimand="ATT"
fun.GBM.wgt <- function(data) {
  fit <- ps(index ~ X1 + X2, data=data, ntrees=5000, interaction.depth=3, shrinkage=0.01,
            perm.test.iters=0, stop.method=c("es.mean", "ks.max"), estimand="ATT", verbose=F)
  wgt <- get.weights(fit, stop.method="es.mean")
  
  return(wgt)
}

## BART outcome model
fun.BART.outcome <- function(data, TRUTH, model_name) {
  
  getMDs <- function(i){
    pred <- bart.out$yhat.test[i,]
    placebo <- pred[which(xtest$trt==0)]
    trt <- pred[which(xtest$trt==1)]
    
    MD <- mean(trt, na.rm=T) - mean(placebo, na.rm=T)
  }
  
  vars <- c("X1","X2")
  data$X2 <- as.factor(data$X2)
  # Training X
  xtrain <- data[data$index==0,c(vars, "trt")]
  xtrain$trt <- as.factor(xtrain$trt)
  ytrain <- data$Y[data$index==0]
  # Test X
  xtest <- data[data$index==1, vars]
  xtest <- rbind(xtest,xtest)
  xtest$trt <- as.factor( rep(c(0,1), each=nrow(xtest)/2) )
  
  bart.out <- bart(x.train=xtrain, y.train=ytrain, keeptrainfits=FALSE, ntree=200,
                   x.test=xtest, verbose=T, ndpost=1000, nskip=100)
  MDs <- sapply(1:1000, getMDs)  # 1000 = ndpost; calculate MD for each iteration
  res <- c(mean(MDs), sd(MDs), q.025(MDs), q.975(MDs))
  res <- c(res, ifelse( res[3]<TRUTH & TRUTH<res[4],1,0 ))
  names(res) <- paste0(model_name,".",c("est","se","low","up","cpind")) 
  
  return(res)
}

## Calculate weights, SuperLearner
fun.SL.select <- function(data) {
  
  vars <- c("X1","X2")
  data$X2 <- as.factor(data$X2)

  ## Luedtke and van der Laan (2016)
  ## SL.rpart causes an issue here; remove it
  ## SL.nnet has too much randomness so I remove it
  #SL.library <- c("SL.gam", "SL.glm", "SL.glm.interaction", "SL.nnet", "SL.rpart")
  SL.library <- c("SL.gam", "SL.glm", "SL.glm.interaction")
  predsSL <- SuperLearner::SuperLearner(
    Y=data$index,
    X=data[,vars],
    SL.library=SL.library,
    family="binomial",
    method="method.NNLS",
    cvControl=list(V=5, stratifyCV=TRUE, shuffle=TRUE) )$SL.predict
  
  wgt <- predsSL/(1-predsSL)
  return(wgt)
}

fun.TMLE.pop <- function(data, TRUTH, model_name){
  # Need index2: 1 for trial 0 for population
  data$index2[data$index==0] <- 1  # trial
  data$index2[data$index==1] <- 0  # population
  
  Y <- data$Y
  A <- ifelse(!is.na(data$trt),data$trt,sample(c(0,1),1))
  X1 <- data$X1
  X2 <- data$X2
  index2 <- data$index2
  n.dat <- nrow(data)
  
  ## step 1 with super learner
  SL.library <- c("SL.glm", "SL.step", "SL.glm.interaction")
  Yall <- Y
  Xall <- data.frame(A,X1,X2)
  Xall$X2 <- as.factor(Xall$X2)
  
  X1new <- X0new <- Xall
  X1new$A <- 1
  X0new$A <- 0
  newdata <- rbind(Xall,X1new,X0new)
  
  X <- Xall[index2==1,]
  y <- Yall[index2==1]
  
  Qinit <- SuperLearner(Y=y, X=X, newX=newdata, SL.library=SL.library, family="gaussian")
  Q0 <- Qinit$SL.predict[1:nrow(Xall)]
  Q0_1 <- Qinit$SL.predict[(nrow(Xall)+1):(2*nrow(Xall))]
  Q0_0 <- Qinit$SL.predict[(2*nrow(Xall)+1):(3*nrow(Xall))]
  
  q <- data.frame("q0"=Q0, "q0_0"=Q0_0, "q0_1"=Q0_1)
  
  ## step 2: H function
  # 1. membership model 
  cps <- predict(glm(index2~X1+X2, data=data, family="binomial"), type="response")
  # 2. treatment assignment model
  cpz <- rep(0.5, n.dat)  ## Because our trial is RCT
  ps0 <- mean(I(index2==0))
  # calculate clever covariate
  g0w <- (1-cpz)*cps/(1-cps)
  g1w <- (cpz*cps)/(1-cps)
  h0w <- ((1-A)*I(index2==1))/g0w
  h1w <- (A*I(index2==1))/g1w
  
  ## step 3: estimate epsilon
  data_q <- cbind(q,Y,index2,g0w,g1w,h0w,h1w)
  epsilon <- coef( glm(Y ~ -1 + offset(Q0) + h0w + h1w, data=data_q, family="gaussian", subset=index2==1)  )
  #epsilon <- coef( betareg(Yuni ~ -1 + offset(logitQ0) + h0w + h1w, data=data_q, subset=index2==1)  )
  
  ## step 4 : update initial prediction
  q1 <- data.frame("q1"=q$q0 + epsilon[1]*data_q$h0w + epsilon[2]*data_q$h1w,
                   "q1_0" = q$q0_0 + epsilon[1]/data_q$g0w,
                   "q1_1" = q$q0_1 + epsilon[2]/data_q$g1w,
                   "index2" = index2)
  tmleest <- mean(q1$q1_1[q1$index2==0]) - mean(q1$q1_0[q1$index2==0])
  
  ## step 5: inference (estimate variance of tmleest)
  # Use efficient influence curve
  aa <- ((A*h1w/ps0) - ((1-A)*h0w/ps0))
  bb <- (Y - q[,1])
  bb <- ifelse(is.na(bb),0,bb)
  eic <- (aa*bb) + (I(index2==0)/ps0*q1[,3]) - (I(index2==0)/ps0*q1[,2]) - (tmleest/ps0)
  se.eic <- sqrt(var(eic)/n.dat)
  
  out <- c(tmleest, se.eic, tmleest-(1.96*se.eic),  tmleest+(1.96*se.eic))
  out <- c(out, ifelse(out[3]<TRUTH & TRUTH<out[4],1,0))
  names(out) <- paste0(model_name,".",c("est", "se", "low", "up", "cpind"))
  return(out)
}

## Calculate SMD with weights
fun.SMD.w <- function(data, model_name){
  pop <- data[data$index==1,]
  rct <- data[data$index==0,]
  
  x.pop <- mean(pop$X1)
  s2.pop <- var(pop$X1)
  p.pop <- mean(pop$X2)
  
  x.rct <- sum(rct$wgt*rct$X1)/sum(rct$wgt)
  s2.rct <- ((sum(rct$wgt)) / (sum(rct$wgt)^2 - sum(rct$wgt^2)))*sum((rct$wgt*(rct$X1 - sum(rct$wgt*rct$X1)/sum(rct$wgt)))^2)
  p.rct <- sum(rct$wgt*rct$X2)/sum(rct$wgt)
  
  SMD <- c( abs(x.rct - x.pop)/sqrt(s2.pop),
            abs(p.rct - p.pop)/sqrt(p.pop*(1-p.pop)))
  names(SMD) <- paste0(model_name,".SMD.",c("x1","x2"))
  return(SMD)
}

## Meta-analysis with unweighted results
fun.meta.unw <- function(data, TRUTH, model_name){
  
  Summary <- data %>%
    group_by(study, trt) %>%
    summarize(n=sum(!is.na(Y)), mean=mean(Y), var=var(Y))
  
  meta.data <- data.frame(
    study <- Summary$study[Summary$trt==0],
    n.placebo <- Summary$n[Summary$trt==0],
    mean.placebo <- Summary$mean[Summary$trt==0],
    var.placebo <- Summary$var[Summary$trt==0],
    n.trt <- Summary$n[Summary$trt==1],
    mean.trt <- Summary$mean[Summary$trt==1],
    var.trt <- Summary$var[Summary$trt==1] )
  colnames(meta.data) <- c("study","n.placebo","mean.placebo","var.placebo","n.trt","mean.trt","var.trt")
  
  ma <- metacont(n.trt, mean.trt, sqrt(var.trt), n.placebo, mean.placebo, sqrt(var.placebo),
                 data=meta.data, studlab=study)
  
  FEres <- data.frame(summary(ma)$fixed)[c("TE", "seTE","lower", "upper")]
  REres <- c(data.frame(summary(ma)$random)[c("TE", "seTE","lower", "upper")], summary(ma)$tau$TE)
  
  MA <- cbind( FEres, ifelse( FEres[3]<TRUTH & TRUTH<FEres[4],1,0 ),
               REres, ifelse( REres[3]<TRUTH & TRUTH<REres[4],1,0 ))  
  names(MA) <- paste0(model_name,".",c("FMA.est","FMA.se","FMA.low","FMA.up","FMA.cpind",
                                       "RMA.est","RMA.se","RMA.low","RMA.up","RMA.tau","RMA.cpind")) 
  return(MA)
}

## Meta-analysis with weighted results
fun.meta.w <- function(data, TRUTH, model_name){
  
  meta.data <- data.frame(fun.weighted.ate2(data))
  ma <- metagen(SATE,SE,studlab=study,data=meta.data)
  
  FEres <- data.frame(summary(ma)$fixed)[c("TE", "seTE","lower", "upper")]
  REres <- c(data.frame(summary(ma)$random)[c("TE", "seTE","lower", "upper")], summary(ma)$tau$TE)
  
  MA <- cbind( FEres, ifelse( FEres[3]<TRUTH & TRUTH<FEres[4],1,0 ),
               REres, ifelse( REres[3]<TRUTH & TRUTH<REres[4],1,0 ))  
  names(MA) <- paste0(model_name,".",c("FMA.est","FMA.se","FMA.low","FMA.up","FMA.cpind",
                                       "RMA.est","RMA.se","RMA.low","RMA.up","RMA.tau","RMA.cpind")) 
  return(MA)
}

#### Function for pairwise meta-analysis with mean difference
## "data" should include
## ES: effect size (either mean difference or SMD)
## se_ES: standard error of ES
pairwise_meta <- function(ES, se_ES) {
  
  df <- length(ES)-1  ## the total number of studies minus 1
  
  var_ES <- se_ES*se_ES
  w <- 1/var_ES
  w_ES <- ES*w
  w_ES2 <- ES*ES*w
  w2 <- w*w
  
  sum.w <- sum(w)
  sum.w.ES <- sum(w_ES)
  sum.w.ES2 <- sum(w_ES2)
  sum.w2 <- sum(w2)
  
  Q <- sum.w.ES2 - (sum.w.ES^2)/sum.w
  I2 <- ((Q-df)/Q)*100
  
  ## Fixed effect
  FE_es <- sum.w.ES/sum.w
  FE_se <- sqrt(1/sum.w)
  FE_low <- FE_es - (1.96*FE_se)
  FE_up <- FE_es + (1.96*FE_se)
  
  ## Random effect
  if(Q > df) {RE_v <- (Q-df)/(sum.w - (sum.w2/sum.w))}
  if (Q <= df) {RE_v <- 0}
  
  wv <- 1/(var_ES + RE_v)
  wv_ES <- ES*wv
  wv_ES2 <- ES*ES*wv
  wv2 <- wv*wv  
  
  sum.wv <- sum(wv)
  sum.wv.ES <- sum(wv_ES)
  sum.wv.ES2 <- sum(wv_ES2)
  sum.wv2 <- sum(wv2) 
  
  Qv <- sum.wv.ES2 - (sum.wv.ES^2)/sum.wv
  I2v <- ((Qv-df)/Qv)*100
  
  RE_es <- sum.wv.ES/sum.wv
  RE_se <- sqrt(1/sum.wv)
  RE_low <- RE_es - (1.96*RE_se)
  RE_up <- RE_es + (1.96*RE_se)
  tau2 <- RE_v
  
  output <- c(Q, I2, FE_es, FE_se, FE_low, FE_up,
              Qv, I2v, RE_es, RE_se, RE_low, RE_up, tau2)
  names(output) <- c("Q", "I2", "FE_es", "FE_se", "FE_low", "FE_up",
                     "Qv", "I2v", "RE_es", "RE_se", "RE_low", "RE_up", "tau2")
  
  return(output)
}

## Obtain results
## data1 = data list for RCT
## data2 = data list for RCT+population
fun.results <- function(data1, data2, TRUTH, model_name){
  #### a. Estimate SATE
  SATE <- unlist(mapply(fun.weighted.ate, data1, TRUTH, paste0(model_name,1:5), SIMPLIFY=F, USE.NAMES=F))
  #### b. Calculate SMD 
  SMD <- unlist(mapply(fun.SMD.w, data2, paste0(model_name,1:5), SIMPLIFY=F, USE.NAMES=F))
  #### c. Distribution of weights
  WGT <- unlist( sapply(1:5, function(X) c(mean(data1[[X]]$wgt), sd(data1[[X]]$wgt), min(data1[[X]]$wgt), max(data1[[X]]$wgt), 
                                           q.025(data1[[X]]$wgt), q.975(data1[[X]]$wgt), quantile(data1[[X]]$wgt, probs=c(0.25, 0.5, 0.75))),
                        simplify=F, USE.NAMES=F))
  names(WGT) <- paste0(model_name,rep(1:5, each=9),".wgt.",c("mean","sd","min","max","low","up","q1","q2","q3"))
  #### d. Meta-analysis
  data1_sub <- data1
  data1_sub <- bind_rows(data1_sub)
  MA <- fun.meta.w(data1_sub, TRUTH, paste0(model_name,"1"))
  #### e. Final output
  res <- unlist(c(SATE, SMD, WGT, MA))
  return(res)
}

#########
#### Run Simulation
#########

fun.single.simulation <- function(...) {
  
  #### True target population & observed target population
  POP_true <- fun.gen.true.pop(n=N_pop_true, x1_mean=x1_mean, x1_sd=x1_sd, x2_p=x2_p, b0=b0, b1=b1, b2=b2, tau=tau, phi1=phi1, phi2=phi2, y_sd=y_sd)
  data_list <- fun.sam.pop(POP_true, N_pop_obs)
  data_pop <- as.data.frame(data_list[1])
  data_sam <- as.data.frame(data_list[2])
  
  #### RCTs
  RCT <- fun.sam.rct.all(data_sam, gam0, gam1, gam2)
  RCT$Y <- RCT$Y0*(1-RCT$trt) + RCT$Y1*RCT$trt
  
  #### Combine population and RCTs
  data_pop$trt <- NA; data_pop$study <- 0; data_pop$Y <- NA; data_pop$index <- 1
  RCT$index <- 0
  
  RCT_list <- split_tibble(RCT, 'study')
  names(RCT_list) <- paste(1:5)
  
  data_list <- lapply(names(RCT_list), function(X) rbind(RCT_list[[X]], data_pop))   
  names(data_list) <- paste(1:5)
  
  #### 0. Basic information of each dataset
  pop_info <- c( "n.pop"=nrow(data_pop), 
                 "mean.pop.x1"=mean(data_pop$X1), "sd.pop.x1"=sd(data_pop$X1), "min.pop.x1"=min(data_pop$X1), "max.pop.x1"=max(data_pop$X1),
                 "low.pop.x1"=q.025(data_pop$X1), "up.pop.x1"=q.975(data_pop$X1), 
                 "q1.pop.x1"=q.25(data_pop$X1), "q2.pop.x1"=median(data_pop$X1), "q3.pop.x1"=q.75(data_pop$X1),
                 "mean.pop.x2"=mean(data_pop$X2))
  names(pop_info) <- c("n.pop", "mean.pop.x1", "sd.pop.x1", "min.pop.x1", "max.pop.x1", "low.pop.x1", "up.pop.x1", "q1.pop.x1", "q2.pop.x1","q3.pop.x1", "mean.pop.x2")
  rct_info <- unlist( sapply(1:5, function(X) c( "n.rct"=nrow(RCT_list[[X]]), 
                                                 "mean.rct.x1"=mean(RCT_list[[X]]$X1), "sd.rct.x1"=sd(RCT_list[[X]]$X1), "min.rct.x1"=min(RCT_list[[X]]$X1), "max.rct.x1"=max(RCT_list[[X]]$X1),
                                                 "low.rct.x1"=q.025(RCT_list[[X]]$X1), "up.rct.x1"=q.975(RCT_list[[X]]$X1), 
                                                 "q1.rct.x1"=q.25(RCT_list[[X]]$X1), "q2.rct.x1"=median(RCT_list[[X]]$X1), "q3.rct.x1"=q.75(RCT_list[[X]]$X1),
                                                 "mean.rct.x2"=mean(RCT_list[[X]]$X2)),
                             simplify=F, USE.NAMES=F))
  names(rct_info) <- paste0( c("n","mean","sd","min","max","low","up","q1","q2","q3","mean"),".rct",rep(1:5, each="11"),c("",rep(".x1",9),".x2"))
  basic_res <- c("TRUTH"=TRUTH, pop_info, rct_info)  
  
  #### 1. Naive approach
  ######
  #### a. Estimate SATE
  naive_SATE <- unlist(mapply(fun.naive.ate, RCT_list, TRUTH, paste0("naive",1:5), SIMPLIFY=F, USE.NAMES=F))
  #### b. Calculate SMD
  sum_pop <- data.frame("x.pop"=mean(data_pop$X1), "s2.pop"=var(data_pop$X1), "p.pop"=mean(data_pop$X2))
  sum_RCT <- lapply( names(RCT_list), function(X) c( "x.rct"=mean(RCT_list[[X]]$X1), "s2.rct"=var(RCT_list[[X]]$X1), "p.rct"=mean(RCT_list[[X]]$X2) ))
  naive_SMD <- unlist(sapply(1:5, function(X) c( abs(sum_RCT[[X]]["x.rct"] - sum_pop["x.pop"])/sqrt(sum_pop["s2.pop"]),
                                                 abs(sum_RCT[[X]]["p.rct"] - sum_pop["p.pop"])/sqrt(sum_pop["p.pop"]*(1-sum_pop["p.pop"]))),
                             simplify=F, USE.NAMES=F))
  names(naive_SMD) <- paste0("naive",rep(1:5,each=2),".SMD.",rep(c("x1","x2"),5))
  #### c. Meta-analysis
  naive_MA <- fun.meta.unw(RCT, TRUTH, "naive")
  #### d. Final output
  naive_res <- unlist(c(naive_SATE, naive_SMD, naive_MA))
  
  ######
  #### 2. logistic regression
  data_wgt <- sapply(data_list, fun.logistic.wgt)
  data_full <- lapply(1:5, function(X) cbind(data_list[[X]], "wgt"=data_wgt[[X]]))
  data_rct <- lapply(1:5, function(X) cbind(data_list[[X]][which(data_list[[X]]$index==0),], "wgt"=data_wgt[[X]][which(data_list[[X]]$index==0)]))
  logistic_res <- fun.results(data_rct, data_full, TRUTH, "logistic")
  
  ######
  #### 3. GBM
  #### Use twang package.
  #### Weights for ATT are 1 for the treatment cases (population) and p/(1-p) for the control cases (sample)
  #### So, I use estimand="ATT"
  data_wgt <- sapply(data_list, fun.GBM.wgt)
  data_full <- lapply(1:5, function(X) cbind(data_list[[X]], "wgt"=data_wgt[[X]]))
  data_rct <- lapply(1:5, function(X) cbind(data_list[[X]][which(data_list[[X]]$index==0),], "wgt"=data_wgt[[X]][which(data_list[[X]]$index==0)]))
  GBM_res <- fun.results(data_rct, data_full, TRUTH, "GBM")
  
  ######
  #### 4. BART predicting outcomes
  SATE_BARTo <- unlist(mapply(fun.BART.outcome, data_list, TRUTH, paste0("BARTo",1:5), SIMPLIFY=F, USE.NAMES=F))
  ES_BARTo <- SATE_BARTo[paste0("BARTo",1:5,".est")]
  se_ES_BARTo <- SATE_BARTo[paste0("BARTo",1:5,".se")]
  MA_BARTo <-  pairwise_meta(ES_BARTo, se_ES_BARTo)[c("FE_es","FE_se", "FE_low", "FE_up","RE_es", "RE_se", "RE_low", "RE_up", "tau2")]
  MA_BARTo["tau2"] <- sqrt(MA_BARTo["tau2"])
  MA_BARTo <- c(MA_BARTo, ifelse(MA_BARTo[3]<TRUTH & TRUTH<MA_BARTo[4],1,0), ifelse(MA_BARTo[7]<TRUTH & TRUTH<MA_BARTo[8],1,0) )
  names(MA_BARTo) <- paste0("BARTo.", c("FMA.est","FMA.se","FMA.low","FMA.up",
                                        "RMA.est","RMA.se","RMA.low","RMA.up","RMA.tau", 
                                        "FMA.cpind", "RMA.cpind"))
  BARTo_res <- c(SATE_BARTo, MA_BARTo) 
  
  ######
  #### 5. Superlearner (## Luedtke and van der Laan (2016))
  #### Remove SL.rpart because SL.rpart causes an issue here. 
  #### SL.nnet has too much randomness so I remove it
  data_wgt <- sapply(data_list, fun.SL.select)
  data_full <- lapply(1:5, function(X) cbind(data_list[[X]], "wgt"=data_wgt[[X]]))
  data_rct <- lapply(1:5, function(X) cbind(data_list[[X]][which(data_list[[X]]$index==0),], "wgt"=data_wgt[[X]][which(data_list[[X]]$index==0)]))
  SL_res <- fun.results(data_rct, data_full, TRUTH, "SL")
  
  ######
  #### 6. TMLE predict outcomes from population subjects
  #### See Rudolph et al. 2017 
  #### Linear fluctuation (Gruber et al. 2010)
  #### Linear fluctuation performs good when Q0 (initial outcome model) is correctly estimated.
  #### I am using SuperLearner for Q0 estimation.
  SATE_TMLEpop <- unlist(mapply(fun.TMLE.pop, data_list, TRUTH, paste0("TMLEpop",1:5), SIMPLIFY=F, USE.NAMES=F))   
  ES_TMLEpop <- SATE_TMLEpop[paste0("TMLEpop",1:5,".est")]
  se_ES_TMLEpop <- SATE_TMLEpop[paste0("TMLEpop",1:5,".se")]
  MA_TMLEpop <-  pairwise_meta(ES_TMLEpop, se_ES_TMLEpop)[c("FE_es","FE_se", "FE_low", "FE_up","RE_es", "RE_se", "RE_low", "RE_up", "tau2")]
  MA_TMLEpop["tau2"] <- sqrt(MA_TMLEpop["tau2"])
  MA_TMLEpop <- c(MA_TMLEpop, ifelse(MA_TMLEpop[3]<TRUTH & TRUTH<MA_TMLEpop[4],1,0), ifelse(MA_TMLEpop[7]<TRUTH & TRUTH<MA_TMLEpop[8],1,0) )
  names(MA_TMLEpop) <- paste0("TMLEpop.", c("FMA.est","FMA.se","FMA.low","FMA.up",
                                            "RMA.est","RMA.se","RMA.low","RMA.up","RMA.tau",
                                            "FMA.cpind", "RMA.cpind"))
  TMLEpop_res <- c(SATE_TMLEpop, MA_TMLEpop)  
  
  ## Combine results
  res <- c(basic_res, 
           naive_res, logistic_res, GBM_res, SL_res,
           BARTo_res, TMLEpop_res)
  return(res)
}