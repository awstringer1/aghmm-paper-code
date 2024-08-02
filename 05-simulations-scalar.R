### Simulations to replicate Section 3.3 of the supplementary materials.
### Alex Stringer
### 2023/08

## Set paths ##
# CHANGE the base path to whatever you want on your machine
basepath <- getwd()
stopifnot(dir.exists(basepath))
# CHANGE the name of the simulation to control how saved results are named
simname <- "sims-scalar-20230817-v1"
resultspath <- file.path(basepath,'results')
if (!dir.exists(resultspath)) dir.create(resultspath)
simresultsname <- paste0(simname,".RData")
simsprocessedname <- paste0(simname,".csv")


## Load Packages ##

pkgs <- c(
  'ggplot2',
  'dplyr',
  'tidyr',
  'fastmatrix',
  'readr',
  'lme4',
  'Rcpp',
  'RcppEigen'
  # 'parallel'
)
for (pkg in pkgs) {
  if (!require(pkg,character.only = TRUE,quietly = TRUE)) {
    cat(paste0("Could not find package ",pkg,", installing from CRAN.\n"))
    install.packages(pkg)
    require(pkg,character.only = TRUE,quietly = TRUE)
  }
}

## Install aghqmm

# Local
# If installing remotely from Github doesn't work, download the package repository to your basepath
# and uncomment the next two lines:
# aghqmmpath <- file.path(basepath, "aghqmm")
# install.packages(aghqmmpath, repos=NULL, type="source")
# Remote: if you have a Github PAT set up in your R session, this should work:
remotes::install_github("awstringer1/aghqmm", force = TRUE)
# If you want to set up remotes, this tutorial is helpful:
# https://carpentries.github.io/sandpaper-docs/github-pat.html
library(aghqmm)


## Set Parameters ##
# These settings reproduce what's in the paper:
# numruns <- 10  # Number of times to execute the simulations
# numsims <- 100 # Number of simulations in each category PER RUN
# m <- c(100,200,500,1000)
# n <- c(3,5,7,9)
# k <- c(1,3,5,7,9,11,13,15,17,19,21,23,25)
# These settings are used for testing/continuous integration purposes:
numruns <- 2  # Number of times to execute the simulations
numsims <- 2 # Number of simulations in each category PER RUN
m <- c(100, 1000)
n <- c(3,5)
k <- c(1,3,5)

beta <- c(-2.5,-.15)
S <- 2
bfgsdelta <- c(1e-01)        
inner_tol <- c(1e-06)
inner_maxitr <- c(10)
bfgshist <- c(6)
past <- c(3)

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
  past = past
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
    doprofile = simstodoframe[i,"doprofile"],
    beta = beta,
    S = S
  )
  
  simlist <- c(simlist,list(lst))
  idx <- idx+1
}

# options(mc.cores = parallel::detectCores())
# RNGkind("L'Ecuyer-CMRG") # For reproducibility with parallel

### Function to execute simulation ###

dosim <- function(lst) {
  # lst: simulation parameters
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
  
  # Fixed parameters
  maxitr <-         100   # Maximum iterations for outer optimization
  tol <-            1e-05 # Tolerance for gradient norm of outer optimization
  max_linesearch <- 100   # Maximum number of iterations for the line search at each step 
  
  cat("Sim: ",lst$idx," of ",length(simlist),"|m=",m,"|n=",n,"|k=",k,"...",sep="")
  
  # Create control list
  bfgscontrol <- aghqmm::aghqmm_control(
    bfgsdelta = bfgsdelta,
    inner_tol = inner_tol,
    inner_maxitr = inner_maxitr,
    bfgshist = bfgshist,
    past = past
  )
  
  # Simulate data
  simdata <- aghqmm::simulate_data(m,n,beta,S)
  # Fit model
  opt <- tryCatch(aghqmm::aghqmm(y ~ x + (1|id),simdata,k=k,method = "both",control = bfgscontrol),error = function(e) e)
  
  if (inherits(opt,'condition')) return(opt)
  
  # Return results
  simresults <- data.frame(
    m=m,n=n,k=k,inner_tol=inner_tol,inner_maxitr=inner_maxitr,bfgshist=bfgshist,past=past,
    totaltime = opt$comptime
  )
  # get the true sigmas from delta/phi
  truesigmas <- S
  estsigmas <- opt$sigmaints[2]
  
  paramresults <- data.frame(
    true = c(beta,truesigmas),
    est = c(opt$theta[1:2],estsigmas),
    lowerWald = c(opt$betaints[ ,1],opt$sigmaints[1]),
    upperWald = c(opt$betaints[ ,3],opt$sigmaints[3])
  ) %>%
    dplyr::mutate(
      covrWald = lowerWald <= true & upperWald >= true,
      lengthWald = upperWald - lowerWald
    )
  cat(" completed.\n",sep="")
  list(simresults=simresults,paramresults=paramresults)
}

processsimulation <- function(sim) {
  if (inherits(sim,'condition')) return(NULL)
  # Change this if the parameters change
  out <- sim$simresults
  out$beta1true <- sim$paramresults[1,'true']
  out$beta1est <- sim$paramresults[1,'est']
  out$beta1bias <- out$beta1est-out$beta1true
  out$beta1covrWald <- sim$paramresults[1,'covrWald']
  out$beta1lengthWald <- sim$paramresults[1,'lengthWald']
  out$beta2true <- sim$paramresults[2,'true']
  out$beta2est <- sim$paramresults[2,'est']
  out$beta2bias <- out$beta2est-out$beta2true
  out$beta2covrWald <- sim$paramresults[2,'covrWald']
  out$beta2lengthWald <- sim$paramresults[2,'lengthWald']
  out$sigmasqtrue <- sim$paramresults[3,'true']
  out$sigmasqest <- sim$paramresults[3,'est']
  out$sigmasqbias <- out$sigmasqest-out$sigmasqtrue
  out$sigmasqcovrWald <- sim$paramresults[3,'covrWald']
  out$sigmasqlengthWald <- sim$paramresults[3,'lengthWald']
  out
}

### Do Simulations ###
set.seed(92633)
# mc.reset.stream() # Reproducbility in parallel
# Do the simulations
cat("Doing",length(simlist),"simulations...\n")
tm <- Sys.time()
# execute the simulations numruns times
simruns <- list()
length(simruns) <- numruns
for (b in 1:numruns) {
  # simruns[[b]] <- mclapply(simlist,dosim)
  simruns[[b]] <- lapply(simlist,dosim)
}
sims <- Reduce(c,simruns)
simtime <- as.numeric(difftime(Sys.time(),tm,units='secs'))
cat("Finished simulations, they took",simtime,"seconds.\n")
cat("Saving simulations...\n")
save(sims,file=file.path(resultspath,simresultsname))
cat("Saved simulations to file:",file.path(resultspath,simresultsname),"\n")
cat("Processing simulations...\n")
simsprocessed <- dplyr::bind_rows(lapply(sims,processsimulation))
readr::write_csv(simsprocessed,file=file.path(resultspath,simsprocessedname))
cat("Finished processing simulations.\n")
cat("Wrote results to:",file.path(resultspath,simsprocessedname),"\n")


