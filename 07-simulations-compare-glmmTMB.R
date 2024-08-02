### Simulations to replicate Section 4.3, simulation 2,
### comparison to GLMMadaptive, in the paper.
### Alex Stringer
### 2023/08

## Set paths ##
# CHANGE the base path to whatever you want on your machine
basepath <- getwd()
# CHANGE the name of the simulation to control how saved results are named
simname <- "sims-glmmTMB-20240724-v1"
stopifnot(dir.exists(basepath))
resultspath <- file.path(basepath,'results')
if (!dir.exists(resultspath)) dir.create(resultspath)
figurepath <- file.path(basepath,'figures')
if (!dir.exists(figurepath)) dir.create(figurepath)
simresultsname <- paste0(simname,".RData")
simsprocessedname <- paste0(simname,".csv")

# Everything else in the script should run without changes

## Load Packages ##

pkgs <- c(
  'ggplot2',
  'dplyr',
  'tidyr',
  'fastmatrix',
  'readr',
  'lme4',
  'glmmTMB',
  'Rcpp',
  'RcppEigen',
  'remotes'
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
# Will run numsims simulations in parallel, repeated numruns times in series.
# These settings reproduce what's in the paper:
# numruns <- 1  # Number of times to execute the simulations
# numsims <- 1000 # Number of simulations in each category PER RUN
# m <- c(100,200,500,1000) # m,n chosen from main manuscript
# n <- c(3,5,7,9)
# k <- 1 # glmmTMB only has Laplace
# These settings are used for testing/continuous integration purposes:
numruns <- 1  # Number of times to execute the simulations
numsims <- 2 # Number of simulations in each category PER RUN
m <- c(100, 200) # m,n chosen from main manuscript
n <- c(3, 5)
k <- 1 # glmmTMB only has Laplace



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

# options(mc.cores = parallel::detectCores())
# RNGkind("L'Ecuyer-CMRG") # For reproducibility with parallel

### Function to execute simulation ###

dosim <- function(lst) {
  # lst: simulation parameters
  idx <-          lst$idx           # Keep track of which simulation
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
  
  cat("Sim: ",lst$idx," of ",length(simlist),"|m=",m,"|n=",n,"|k=",k,"... ",sep="")
  
  # Create control list
  bfgscontrol <- aghqmm::aghqmm_control()
  
  # Simulate data
  simdata <- aghqmm::simulate_data(m,n,beta,S)
  # Fit model using both methods
  # cat("doing aghq... ")
  opt_aghqmm <- tryCatch(aghqmm::aghqmm(y ~ x*t + (t|id),simdata,k=k,method = "both",control = bfgscontrol),
                         warning = function(w) w,
                         error = function(e) e)

  # formula <- y ~ x*t + (t|id)
  # data <- simdata
  # k <- 1
  # method <- "glmmTMB"
  # family <- stats::binomial()
  # control <- bfgscontrol


  # cat("doing glmmTMB... ")
  opt_glmmTMB <- tryCatch(aghqmm::aghqmm(y ~ x*t + (t|id),simdata,k=k,method = "glmmTMB",control = bfgscontrol),
                               warning = function(w) w,
                               error = function(e) e)

  cat("completed.\n",sep="")
  list(paramvalues = list(idx=idx,m=m,n=n,k=k,inner_tol=inner_tol,inner_maxitr=inner_maxitr,bfgshist=bfgshist,past=past),
       aghqmm = opt_aghqmm$comptime,
       glmmTMB = opt_glmmTMB$comptime)
}

processsimulation <- function(sim) {
  
  out_aghqmm <- with(sim$paramvalues,data.frame(
    idx = idx,
    method = "aghqmm",
    m = m,
    n = n,
    k = k
  ))
  out_glmmTMB <- with(sim$paramvalues,data.frame(
    idx = idx,
    method = "glmmTMB",
    m = m,
    n = n,
    k = k
  ))
  
  if (inherits(sim$aghqmm,"condition") | is.null(sim$aghqmm)) {
    out_aghqmm['successful'] <- 0
    out_aghqmm['comptime'] <- 0
  } else {
    out_aghqmm['successful'] <- 1
    out_aghqmm['comptime'] <- sim$aghqmm
  }
  if (inherits(sim$glmmTMB,"condition") | is.null(sim$glmmTMB)) {
    out_glmmTMB['successful'] <- 0
    out_glmmTMB['comptime'] <- 0
  } else {
    out_glmmTMB['successful'] <- 1
    out_glmmTMB['comptime'] <- sim$glmmTMB
  }
  
  out <- rbind(out_aghqmm,out_glmmTMB)
  
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
#load(file=file.path(resultspath,simresultsname))
cat("Saved simulations to file:",file.path(resultspath,simresultsname),"\n")
cat("Processing simulations...\n")
simsprocessed <- dplyr::bind_rows(lapply(sims,processsimulation))
readr::write_csv(simsprocessed,file=file.path(resultspath,simsprocessedname))
cat("Finished processing simulations.\n")
cat("Wrote results to:",file.path(resultspath,simsprocessedname),"\n")
cat("Summarizing simulations... ")

# check successes, not for paper, just to make sure it worked
simsprocessed %>%
  group_by(method, m, n) %>%
  summarize(numtotal = n(), numsucces = sum(successful)) %>%
  print(n = Inf)
# Almost all successful.

## Relative Computation Times ##

aghq_times <- simsprocessed %>%
  filter(method == "aghqmm") %>%
  dplyr::select(idx,aghqtime = comptime)

reltimes <- simsprocessed %>%
  left_join(aghq_times,by='idx') %>%
  mutate(reltime = comptime / aghqtime)


m_labeller <- function(vl) paste0("m = ",vl)
n_labeller <- function(vl) paste0("n = ",vl)

relative_timeplot <- reltimes %>%
  filter(successful == 1, method == "glmmTMB") %>%
  ggplot(aes(y = reltime)) +
  theme_bw() +
  facet_grid(n~m,labeller = labeller(m = m_labeller,n = n_labeller)) +
  geom_boxplot(outlier.size = .7,outlier.alpha = .05) +
  labs(
    title = "Relative Computation Times, glmmTMB / New Approach",
    x = "Laplace Approximation",
    y = "Relative computation time, glmmTMB / New Approach"
  ) +
  scale_y_continuous(breaks = seq(0,20,by=2)) +
  coord_cartesian(ylim = c(0,20)) +
  scale_x_discrete(breaks = seq(1,25,by=4)) +
  geom_hline(yintercept=1,lty='dashed')

ggsave(filename = file.path(figurepath,"glmmTMB-relativetimeplot.pdf"),plot = relative_timeplot,width=7,height=7)  

