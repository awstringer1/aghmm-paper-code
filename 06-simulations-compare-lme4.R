# AGHQ using the C++ L-BFGS-B implementation
# through the aghqmm package
# Comparison to lme4 for a random-intercepts model

## Set paths ##
# CHANGE the base path to whatever you want on your machine
basepath <- '~/work/projects/mixedmodel-computation/replication'
# CHANGE the name of the simulation to control how saved results are named
simname <- "sims-lme4-20230817-v1"
stopifnot(dir.exists(basepath))
resultspath <- file.path(basepath,'results')
if (!dir.exists(resultspath)) dir.create(resultspath)
simresultsname <- paste0(simname,".RData")
simsprocessedname <- paste0(simname,".csv")

# CHANGE the path to where you downloaded the aghqmm package repo
# from https://github.com/awstringer1/aghqmm
aghqmmpath <- "~/work/projects/mixedmodel-computation/aghqmm" # CHANGE this

# Everything else in the script should run without changes

## Load Packages ##

install.packages(aghqmmpath,repos=NULL,type="source")

pkgs <- c(
  'tidyverse',
  'lme4',
  'Rcpp',
  'RcppEigen',
  'parallel'
)
for (pkg in pkgs) {
  if (!require(pkg,character.only = TRUE,quietly = TRUE)) {
    cat(paste0("Could not find package ",pkg,", installing from CRAN.\n"))
    install.packages(pkg)
    require(pkg,character.only = TRUE,quietly = TRUE)
  }
}


## Set Parameters ##
numruns <- 2  # Number of times to execute the simulations
numsims <- 10 # Number of simulations in each category PER RUN
m <- c(100,200,500,1000) # m,n chosen from main manuscript
n <- c(3,5,7,9)
k <- seq(1,25,by=2)
beta <- c(-2.5,-.15)
S <- 2 # random intercept variance

# bfgsdelta <- c(1e-06)
bfgsdelta <- c(1e-01)
inner_tol <- c(1e-06)
inner_maxitr <- c(10)
bfgshist <- c(6)
past <- c(3)
wolfe <- .9
ftol <- 1e-04
max_linesearch <- 100
tol <- 1e-04
maxitr <- 100

# Create simulation objects
simstodoframe <- expand.grid(
  sim = 1:numsims,
  n = n,
  m = m,
  k = k,
  bfgsdelta = bfgsdelta,
  inner_tol = inner_tol,
  inner_maxitr = inner_maxitr,
  bfgshist = bfgshist,
  past = past,
  ftol = ftol,
  wolfe = wolfe,
  max_linesearch = max_linesearch,
  tol = tol,
  maxitr = maxitr
) %>%
  as_tibble() %>%
  as.data.frame()

simlist <- list()
idx <- 1
for (i in 1:nrow(simstodoframe)) {
  lst <- list(
    idx = idx,
    sim = simstodoframe[i,"sim"],
    n = simstodoframe[i,"n"],
    m = simstodoframe[i,"m"],
    k = simstodoframe[i,"k"],
    bfgsdelta = simstodoframe[i,"bfgsdelta"],
    inner_tol = simstodoframe[i,"inner_tol"],
    inner_maxitr = simstodoframe[i,"inner_maxitr"],
    bfgshist = simstodoframe[i,"bfgshist"],
    past = simstodoframe[i,"past"],
    wolfe = simstodoframe[i,"wolfe"],
    ftol = simstodoframe[i,"ftol"],
    max_linesearch = simstodoframe[i,"max_linesearch"],
    tol = simstodoframe[i,"tol"],
    maxitr = simstodoframe[i,"maxitr"],
    beta = beta,
    S = S
  )
  
  simlist <- c(simlist,list(lst))
  idx <- idx+1
}

options(mc.cores = parallel::detectCores())
RNGkind("L'Ecuyer-CMRG") # For reproducibility with parallel

### Function to execute simulation ###

dosim <- function(lst) {
  # lst: simulation parameters
  idx <-          lst$idx           # Keep track of which simulation
  m <-            lst$m             # Number of groups
  n <-            lst$n             # Number of observations per group
  beta <-         lst$beta          # Regression parameters
  S <-            lst$S
  k <-            lst$k             # Number of quadrature points
  bfgsdelta <-    lst$bfgsdelta     # Tolerance for function change
  inner_tol <-    lst$inner_tol     # Tolerance for inner optimization
  inner_maxitr <- lst$inner_maxitr  # Maximum iterations for inner optimization
  bfgshist <-     lst$bfgshist      # Number of iterations of gradient information to retain for Hessian approximation
  past <-         lst$past          # Number of iterations to look at for convergence based on function change
  ftol <-         lst$ftol        
  wolfe <-        lst$wolfe          
  max_linesearch<-lst$max_linesearch          
  tol <-          lst$tol
  maxitr <-       lst$maxitr
  
  # Fixed parameters
  # maxitr <-         100   # Maximum iterations for outer optimization
  # tol <-            1e-03 # Tolerance for gradient norm of outer optimization
  # max_linesearch <- 100   # Maximum number of iterations for the line search at each step 
  
  cat("Sim: ",lst$idx," of ",length(simlist),"|m=",m,"|n=",n,"|k=",k,"... ",sep="")
  
  # Create control list
  bfgscontrol <- aghqmm::aghqmm_control(
    bfgsdelta = bfgsdelta,
    bfgshist = bfgshist,
    past = past,
    wolfe = wolfe,
    ftol = ftol,
    max_linesearch = max_linesearch,
    tol = tol,
    maxitr = maxitr
  )
  
  # Simulate data
  simdata <- aghqmm::simulate_data(m,n,beta,S)
  # Fit model using both methods
  cat("doing aghq... ")
  opt_aghqmm <- tryCatch(aghqmm::aghqmm(y ~ x + (1|id),simdata,k=k,method = "both",control = bfgscontrol),
                         warning = function(w) w,
                         error = function(e) e)
  
  cat("doing lme4... ")
  opt_lme4 <- tryCatch(aghqmm::aghqmm(y ~ x + (1|id),simdata,k=k,method = "lme4",control = bfgscontrol),
                       warning = function(w) w,
                       error = function(e) e)
  
  totalsamplesize <- nrow(simdata)
  
  if (inherits(opt_aghqmm,'condition')) {
    results_aghqmm <- opt_aghqmm
  } else {
    simresults_aghqmm <- data.frame(
      totaltime = opt_aghqmm$comptime,
      nll = opt_aghqmm$nll,
      nllavg = opt_aghqmm$nll/totalsamplesize,
      nllavgbase10 = (opt_aghqmm$nll/log(10))/totalsamplesize,
      normgrad = opt_aghqmm$normgrad_2,
      normgradlogbase10avg = log(opt_aghqmm$normgrad_2,base=10) - log(totalsamplesize,base=10)
    )
    results_aghqmm <- list(simresults=simresults_aghqmm)
  }
  
  if (inherits(opt_lme4,'condition')) {
    results_lme4 <- opt_lme4
  } else {
    estsigmas_lme4 <- opt_lme4$sigmaints[ ,2]
    simresults_lme4 <- data.frame(
      totaltime = opt_lme4$comptime,
      nll = opt_lme4$nll,
      nllavg = opt_lme4$nll/totalsamplesize,
      nllavgbase10 = (opt_lme4$nll/log(10))/totalsamplesize,
      normgrad = opt_lme4$normgrad_2,
      normgradlogbase10avg = log(opt_lme4$normgrad_2,base=10) - log(totalsamplesize,base=10)
    )
    results_lme4 <- list(simresults=simresults_lme4)
  }
  cat("completed.\n",sep="")
  list(paramvalues = list(idx=idx,m=m,n=n,k=k,inner_tol=inner_tol,inner_maxitr=inner_maxitr,bfgshist=bfgshist,past=past),
       aghqmm = results_aghqmm,
       lme4 = results_lme4,
       data = simdata,
       opt = list(aghqmm=opt_aghqmm,lme4=opt_lme4))
}

## TESTING: TO DELETE
# lst <- simlist[[sample.int(length(simlist),1)]]
# lst$m <- 1000
# lst$n <- 5
# lst$k <- 1
# 
# # try mucking about with the BFGS settings
# lst$bfgsdelta<-1e-03
# 
# sim <- dosim(lst)
# sim$aghqmm
# sim$lme4
# 
# optfun <- function(tt)
#   aghqmm::logmarglik(tt,y~x+(1|id),sim$data,k=1)$nll
# gradfun <- function(tt)
#   aghqmm::logmarglik(tt,y~x+(1|id),sim$data,k=1)$grad
# 
# tm <- Sys.time()
# outsideopt <- optim(rep(0,3),optfun,gradfun,method="BFGS")
# optimtime <- as.numeric(difftime(Sys.time(),tm,units='secs'))
# cbind(outsideopt$par,sim$opt$lme4$theta)
# c(optimtime,sim$lme4$simresults$totaltime)
# 
# 
# c(sim$aghqmm$simresults$nll,sim$lme4$simresults$nll)
# c(sim$aghqmm$simresults$normgrad,sim$lme4$simresults$normgrad)
# c(sim$aghqmm$simresults$totaltime,sim$lme4$simresults$totaltime)
# 
# library(aghqmm)
# formula <- y ~ x + (1|id)
# data <- sim$data
# k=lst$k
# family=stats::binomial()
# method = "both"
# control = aghqmm_control()






processsimulation <- function(sim) {
  
  out_aghqmm <- with(sim$paramvalues,data.frame(
    idx = idx,
    method = "aghqmm",
    m = m,
    n = n,
    k = k
  ))
  out_lme4 <- with(sim$paramvalues,data.frame(
    idx = idx,
    method = "lme4",
    m = m,
    n = n,
    k = k
  ))
  
  if (inherits(sim$aghqmm,"condition")) {
    out_aghqmm['successful'] <- 0
    out_aghqmm['comptime'] <- 0
    out_aghqmm['nllavgbase10'] <- 0
    out_aghqmm['normgradavg2base10log'] <- 0
  } else {
    out_aghqmm['successful'] <- 1
    out_aghqmm['comptime'] <- sim$aghqmm$simresults$totaltime
    out_aghqmm['nllavgbase10'] <- sim$aghqmm$simresults$nllavgbase10
    out_aghqmm['normgradavg2base10log'] <- sim$aghqmm$simresults$normgradlogbase10avg
  }
  if (inherits(sim$lme4,"condition")) {
    out_lme4['successful'] <- 0
    out_lme4['comptime'] <- 0
    out_lme4['nllavgbase10'] <- 0
    out_lme4['normgradavg2base10log'] <- 0
  } else {
    out_lme4['successful'] <- 1
    out_lme4['comptime'] <- sim$lme4$simresults$totaltime
    out_lme4['nllavgbase10'] <- sim$lme4$simresults$nllavgbase10
    out_lme4['normgradavg2base10log'] <- sim$lme4$simresults$normgradlogbase10avg
  }
  
  out <- rbind(out_aghqmm,out_lme4)
  
  out
}

### Do Simulations ###
set.seed(92633)
mc.reset.stream() # Reproducbility in parallel
# Do the simulations
cat("Doing",length(simlist),"simulations...\n")
tm <- Sys.time()
# execute the simulations numruns times
simruns <- list()
length(simruns) <- numruns
for (b in 1:numruns) {
  simruns[[b]] <- mclapply(simlist,dosim)
}
sims <- Reduce(c,simruns)
simtime <- as.numeric(difftime(Sys.time(),tm,units='secs'))
cat("Finished simulations, they took",simtime,"seconds.\n")
cat("Saving simulations...\n")
save(sims,file=file.path(resultspath,simresultsname))
#load(file=file.path(resultspath,simresultsname))
cat("Saved simulations to file:",file.path(resultspath,simresultsname),"\n")
cat("Processing simulations...\n")
simsprocessed <- dplyr::bind_rows(lapply(sims,processsimulation))
readr::write_csv(simsprocessed,file=file.path(resultspath,simsprocessedname))
cat("Finished processing simulations.\n")
cat("Wrote results to:",file.path(resultspath,simsprocessedname),"\n")


