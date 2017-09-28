#! /usr/bin/env Rscript

# get environment variables
MYSCRATCH <- Sys.getenv('MYSCRATCH')
RESULTDIR <- Sys.getenv('RESULTDIR')
STEPSIZE <- as.numeric(Sys.getenv('STEPSIZE'))
TASKID <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))

# set defaults if nothing comes from environment variables
MYSCRATCH[is.na(MYSCRATCH)] <- '.'
RESULTDIR[is.na(RESULTDIR)] <- '.'
STEPSIZE[is.na(STEPSIZE)] <- 1
TASKID[is.na(TASKID)] <- 0

# get command lines arguments
args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 1){
  stop("Not enough arguments. Please use args 'listsize', 'prepare', 'run <itemsize>' or 'merge'")
}

ns <- c(1000, 2500, 5000)
bigB <- 1000
parm <- expand.grid(seed=1:bigB,
                    n=ns)

# source in simulation Functions
source("~/hamming/makeData.R")
# load survtmle
library(survtmle, lib.loc = "/home/dbenkese/R/x86_64-unknown-linux-gnu-library/3.2/")

# get the list size #########
if (args[1] == 'listsize') {
  cat(nrow(parm))
}

# execute prepare job ##################
if (args[1] == 'prepare') {
  for(i in 1:nrow(parm)){
     set.seed(parm$seed[i])
     dat <- makeData(n=parm$n[i])
     save(dat, file=paste0("~/hamming/scratch/dataList_n=",parm$n[i],
                           "_seed=",parm$seed[i],".RData"))
   }
   print(paste0('initial datasets saved to: ~/hamming/scratch/dataList ... .RData'))
}

# execute parallel job #################################################
if (args[1] == 'run') {
  if (length(args) < 2) {
    stop("Not enough arguments. 'run' needs a second argument 'id'")
  }
  id <- as.numeric(args[2])
  print(paste(Sys.time(), "arrid:" , id, "TASKID:",
              TASKID, "STEPSIZE:", STEPSIZE))
  for (i in (id+TASKID):(id+TASKID+STEPSIZE-1)) {
    print(paste(Sys.time(), "i:" , i))
    print(parm[i,])
    
    # load data
    load(paste0("~/hamming/scratch/dataList_n=",parm$n[i],
                "_seed=",parm$seed[i], ".RData"))
    
    # set seed
    set.seed(parm$seed[i])
    
    # get glm formula for censoring and ftime
    glm.ctime <- get.ctimeForm(trt = dat$trt, site = dat$adjustVars$site, 
                               ftime = dat$ftime, ftype = dat$ftype)

    # faster to call mean.tmle
    object <- survtmle(ftime = dat$ftime,
                      ftype = dat$ftype,
                      adjustVars = dat$adjustVars,
                      trt = dat$trt,
                      glm.trt = "1",
                      glm.ftime = "trt*factor(site)",
                      glm.ctime = glm.ctime,
                      method = "mean",
                      t0=6)

    # get trend estimate
    trend <- trend_test(object)

    # true value
    nabla_g <- grad_g(object$est)
    Upsilon_n <- nabla_g %*% cov(Reduce(cbind, object$ic)) %*% t(nabla_g) 
    # beta_0 <- mean(replicate(20, getTruth(Upsilon = Upsilon_n, n = 1e6)[2]))
    beta_0 <- replicate(20, getTruth(Upsilon = Upsilon_n, n = 1e6))
    rowMeans(beta_0)

    # output should look like 
    # seed, n, truth
    # beta_n, ci, cov
    out <- c(parm$seed[i], parm$n[i], beta_0,
             trend$beta, trend$ci, 
             as.numeric(trend$ci[1] < beta_0 & trend$ci[2] > beta_0))

    # save output 
    save(out, file = paste0("~/hamming/scratch/out_n=",
                            parm$n[i],"_seed=",parm$seed[i],".RData.tmp"))
    file.rename(paste0("~/hamming/scratch/out_n=",
                       parm$n[i],"_seed=",parm$seed[i],".RData.tmp"),
                paste0("~/hamming/scratch/out_n=",
                       parm$n[i],"_seed=",parm$seed[i],".RData"))
  }
}

# merge job ###########################
if (args[1] == 'merge') {   
    ns <- c(1000, 2500, 5000)
    bigB <- 1000
    parm <- expand.grid(seed=1:bigB,
                        n=ns)

    rslt <- NULL
    for(i in 1:nrow(parm)){
        tmp <- tryCatch({
            load(paste0("~/hamming/scratch/out_n=",
                        parm$n[i],"_seed=",parm$seed[i],".RData"))
            out
        }, error=function(e){
          c(parm$seed[i], parm$n[i], rep(NA,8))
        })
        rslt <- rbind(rslt, tmp)
    }
    # format
    out <- data.frame(rslt)
    colnames(out) <- c("seed","n","truth","beta","ci_l","ci_u","cover")
    save(out, file=paste0('~/hamming/out/allOut.RData'))
    print("results saved")
    tmp1 <- by(out, out$n, function(x){ 100*mean(x$beta - x$truth, na.rm = TRUE) })
    tmp2 <- by(out, out$n, function(x){ 100*var(x$beta, na.rm = TRUE) })
    tmp3 <- by(out, out$n, function(x){ 100*mean((x$beta - x$truth)^2, na.rm = TRUE) })
    tmp4 <- by(out, out$n, function(x){ 100*mean(x$cover, na.rm = TRUE) })
    tmp <- Reduce(cbind,list(tmp1,tmp2,tmp3,tmp4))
    library(xtable)
    colnames(tmp) <- c("Bias x 1e2","Variance x 1e2","Mean squared-error x 1e2", "Coverage (%)")
    xtable(tmp)
}