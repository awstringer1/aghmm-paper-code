### Simulations to replicate Section 4.3, simulation 1
### in the paper.
### Alex Stringer
### 2023/08

## Set paths ##
# CHANGE the base path to whatever you want on your machine
basepath <- '~/work/projects/mixedmodel-computation/replication'
# CHANGE the name of the simulation to control how saved results are named
simname <- "sims-20230817-v1"
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
  'GLMMadaptive',
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
# Will run numsims simulations in parallel, repeated numruns times in series.
numruns <- 10  # Number of times to execute the simulations
numsims <- 100 # Number of simulations in each category PER RUN
m <- c(100,200,500,1000)
n <- c(3,5,7,9)
k <- c(1,3,5,7,9,11,13,15,17,19,21,23,25)
beta <- c(-2.5,-.15,.1,.2)
S <- matrix(c(2,1,1,1),ncol=2)
Q <- solve(S)
stopifnot(all(eigen(Q,only.values=TRUE)$values > 0))
Qchol <- fastmatrix::ldl(Q)
stopifnot(abs(with(Qchol, lower %*% diag(d) %*% t(lower)) - Q) < 1e-06)
delta <- log(Qchol$d)
phi <- with(Qchol,lower[lower.tri(lower)])

bfgsdelta <- c(1e-06)        
inner_tol <- c(1e-06)
inner_maxitr <- c(10)
bfgshist <- c(4)
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
    delta = delta,
    phi = phi,
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
  m <-            lst$m             # Number of groups
  n <-            lst$n             # Number of observations per group
  beta <-         lst$beta          # Regression parameters
  delta <-        lst$delta         # Variance parameters: diagonal of Cholesky of precision matrix
  phi <-          lst$phi           # Variance parameters: off-diagonal of Cholesky of precision matrix
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
  bfgscontrol <- aghqmm::aghqmm_control()
  
  # Simulate data
  simdata <- aghqmm::simulate_data(m,n,beta,S)
  # Fit model
  opt <- tryCatch(aghqmm::aghqmm(y ~ x*t + (t|id),simdata,k=k,method = "both",control = bfgscontrol),error = function(e) e)
  
  if (inherits(opt,'condition')) return(opt)
  
  # Return results
  simresults <- data.frame(
    m=m,n=n,k=k,inner_tol=inner_tol,inner_maxitr=inner_maxitr,bfgshist=bfgshist,past=past,
    totaltime = opt$comptime
  )
  # get the true sigmas from delta/phi
  truesigmas <- c(diag(S),S[2,1])
  estsigmas <- opt$sigmaints[ ,2]
  
  paramresults <- data.frame(
    true = c(beta,truesigmas),
    est = c(opt$theta[1:4],estsigmas),
    lowerWald = c(opt$betaints[ ,1],opt$sigmaints[ ,1]),
    upperWald = c(opt$betaints[ ,3],opt$sigmaints[ ,3])
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
  out$beta3true <- sim$paramresults[3,'true']
  out$beta3est <- sim$paramresults[3,'est']
  out$beta3bias <- out$beta3est-out$beta3true
  out$beta3covrWald <- sim$paramresults[3,'covrWald']
  out$beta3lengthWald <- sim$paramresults[3,'lengthWald']
  out$beta4true <- sim$paramresults[4,'true']
  out$beta4est <- sim$paramresults[4,'est']
  out$beta4bias <- out$beta4est-out$beta4true
  out$beta4covrWald <- sim$paramresults[4,'covrWald']
  out$beta4lengthWald <- sim$paramresults[4,'lengthWald']
  out$sigmasq1true <- sim$paramresults[5,'true']
  out$sigmasq1est <- sim$paramresults[5,'est']
  out$sigmasq1bias <- out$sigmasq1est-out$sigmasq1true
  out$sigmasq1covrWald <- sim$paramresults[5,'covrWald']
  out$sigmasq1lengthWald <- sim$paramresults[5,'lengthWald']
  out$sigmasq2true <- sim$paramresults[6,'true']
  out$sigmasq2est <- sim$paramresults[6,'est']
  out$sigmasq2bias <- out$sigmasq2est-out$sigmasq2true
  out$sigmasq2covrWald <- sim$paramresults[6,'covrWald']
  out$sigmasq2lengthWald <- sim$paramresults[6,'lengthWald']
  out$sigmacov1true <- sim$paramresults[7,'true']
  out$sigmacov1est <- sim$paramresults[7,'est']
  out$sigmacov1bias <- out$sigmacov1est-out$sigmacov1true
  out$sigmacov1covrWald <- sim$paramresults[7,'covrWald']
  out$sigmacov1lengthWald <- sim$paramresults[7,'lengthWald']
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
cat("Saved simulations to file:",file.path(resultspath,simresultsname),"\n")
cat("Processing simulations...\n")
simsprocessed <- dplyr::bind_rows(lapply(sims,processsimulation))
readr::write_csv(simsprocessed,file=file.path(resultspath,simsprocessedname))
cat("Finished processing simulations.\n")
cat("Wrote results to:",file.path(resultspath,simsprocessedname),"\n")
cat("Next, execute the script '01-summarize-simulations-absolute.R'\n")

