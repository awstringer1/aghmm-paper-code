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
figurepath <- file.path(basepath,'figures')
if (!dir.exists(figurepath)) dir.create(figurepath)
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


simsummary <- simsprocessed %>%
  mutate(sigmasqrelbias = (sigmasqest-sigmasqtrue)/sigmasqtrue)

successes <- simsummary %>%
  group_by(m,n,k) %>%
  summarize(successful = n())

## Summaries ##

m_labeller <- function(vl) paste0("m = ",vl)
n_labeller <- function(vl) paste0("n = ",vl)

## Bias
biasboxplot <- function(var,data = simsummary,ylim=NULL) {
  vr <- sym(paste0(var,'bias'))
  plt <- ggplot(data,aes(x=as.factor(k),y=!!vr)) +
    facet_grid(n~m,labeller = labeller(m = m_labeller,n = n_labeller)) +
    theme_bw() +
    geom_boxplot(outlier.size = .7,outlier.alpha = .05) +
    labs(x = "Number of Quadrature Points") +
    scale_x_discrete(breaks = seq(1,25,by=4))
  # theme(axis.text.x = element_text(angle=45,size=7))
  if (!is.null(plt)) plt <- plt + coord_cartesian(ylim=ylim)
  plt
}
# Save nicely annotated plots

beta0biasplot <- biasboxplot("beta1",ylim = c(-8,3)) +
  labs(title = bquote(Bias*", "*widehat(beta)[0] - beta[0]),y = bquote(Bias(widehat(beta)[0])))

beta1biasplot <- biasboxplot("beta2",ylim = c(-8,3)) +
  labs(title = bquote(Bias*", "*widehat(beta)[1] - beta[1]),y = bquote(Bias(widehat(beta)[1])))

sigmasqrelbiasplot <- biasboxplot("sigmasqrel",ylim = c(-1,20)) +
  labs(title = bquote(Relative~Bias*", "*widehat(sigma)^2/sigma^2),y = bquote(Relative~Bias(widehat(sigma)^2)))
sigmasqrelbiasplotzoom <- biasboxplot("sigmasqrel",ylim = c(-1,3)) +
  labs(title = bquote(Relative~Bias*", "*widehat(sigma)^2/sigma^2),y = bquote(Relative~Bias(widehat(sigma)^2)))

ggsave(filename = file.path(figurepath,"beta0-bias-scalar.pdf"),plot=beta0biasplot,width=7,height=7)
ggsave(filename = file.path(figurepath,"beta1-bias-scalar.pdf"),plot=beta1biasplot,width=7,height=7)
ggsave(filename = file.path(figurepath,"sigmasq-relbias-scalar.pdf"),plot=sigmasqrelbiasplot,width=7,height=7)
ggsave(filename = file.path(figurepath,"sigmasq-relbias-scalar-zoom.pdf"),plot=sigmasqrelbiasplotzoom,width=7,height=7)




## Coverage
coverageplot <- function(var,data = simsummary,ylim=NULL) {
  vW <- sym(paste0(var,"covrWald"))
  plt <- data %>%
    group_by(m,n,k) %>%
    summarize(covrWald = mean(!!vW,na.rm=TRUE),
              covrWaldsd = sqrt(covrWald*(1-covrWald)/n())
    ) %>%
    ggplot(aes(x=k)) +
    theme_bw() +
    facet_grid(n~m,labeller = labeller(m = m_labeller,n = n_labeller)) +
    geom_line(aes(y=covrWald)) +
    geom_line(aes(y=covrWald-2*covrWaldsd),linetype='dashed') +
    geom_line(aes(y=covrWald+2*covrWaldsd),linetype='dashed') +
    geom_point(aes(y=covrWald),size = 1) +
    scale_x_continuous(breaks = seq(1,25,by=4)) +
    geom_hline(yintercept = .95,linetype = 'dotted') +
    scale_y_continuous(limits = c(0,1),breaks = seq(0,1,by=.1)) +
    labs(x = "Number of Quadrature Points",y="Empirical Coverage Proportion")
  if (!is.null(plt)) {
    plt <- plt + coord_cartesian(ylim=ylim)
  } else {
    plt <- plt + coord_cartesian(ylim=c(0,1))
  }
  plt
}
# nice plots
beta0covrplot <- coverageplot("beta1") + 
  labs(title = bquote(Empirical~Coverage*","~beta[0]))

beta1covrplot <- coverageplot("beta2") + 
  labs(title = bquote(Empirical~Coverage*","~beta[1]))

sigmasqcovrplot <- coverageplot("sigmasq") + 
  labs(title = bquote(Empirical~Coverage*","~sigma^2))

ggsave(filename = file.path(figurepath,"beta0-covr-scalar.pdf"),plot=beta0covrplot,width=7,height=7)
ggsave(filename = file.path(figurepath,"beta1-covr-scalar.pdf"),plot=beta1covrplot,width=7,height=7)
ggsave(filename = file.path(figurepath,"sigmasq-covr-scalar.pdf"),plot=sigmasqcovrplot,width=7,height=7)

## Lengths
lengthboxplot <- function(var,data = simsummary,ylim=NULL) {
  vr <- sym(paste0(var,'lengthWald'))
  plt <- ggplot(data,aes(x=as.factor(k),y=!!vr)) +
    facet_grid(n~m,labeller = labeller(m = m_labeller,n = n_labeller)) +
    theme_bw() +
    geom_boxplot(outlier.size = .7,outlier.alpha = .05) +
    labs(x = "Number of Quadrature Points",y = "Upper Limit - Lower Limit") +
    scale_x_discrete(breaks = seq(1,25,by=4))
  if (!is.null(plt)) plt <- plt + coord_cartesian(ylim=ylim)
  plt
}

beta0lengthplot <- lengthboxplot("beta1",ylim = c(0,10)) + 
  labs(title = bquote(Length~of~Confidence~Intervals*","~beta[0]))

beta1lengthplot <- lengthboxplot("beta2",ylim = c(0,10)) + 
  labs(title = bquote(Length~of~Confidence~Intervals*","~beta[1]))

sigmasqlengthplot <- lengthboxplot("sigmasq",ylim = c(0,100)) + 
  labs(title = bquote(Length~of~Confidence~Intervals*","~sigma^2))
sigmasqlengthplotzoom <- lengthboxplot("sigmasq",ylim = c(0,10)) + 
  labs(title = bquote(Length~of~Confidence~Intervals*","~sigma^2))

ggsave(filename = file.path(figurepath,"beta0-length-scalar.pdf"),plot=beta0lengthplot,width=7,height=7)
ggsave(filename = file.path(figurepath,"beta1-length-scalar.pdf"),plot=beta1lengthplot,width=7,height=7)
ggsave(filename = file.path(figurepath,"sigmasq-length-scalar.pdf"),plot=sigmasqlengthplot,width=7,height=7)
ggsave(filename = file.path(figurepath,"sigmasq-length-scalar-zoom.pdf"),plot=sigmasqlengthplotzoom,width=7,height=7)

