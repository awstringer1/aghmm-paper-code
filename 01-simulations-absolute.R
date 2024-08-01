### Simulations to replicate Section 4.3, simulation 1
### in the paper.
### Alex Stringer
### 2023/08

## Set paths ##
# CHANGE the base path to whatever you want on your machine
basepath <- getwd()
# CHANGE the name of the simulation to control how saved results are named
simname <- "sims-20240801-v1"
stopifnot(dir.exists(basepath))
resultspath <- file.path(basepath,'results')
if (!dir.exists(resultspath)) dir.create(resultspath)
figurepath <- file.path(basepath,'figures')
if (!dir.exists(figurepath)) dir.create(figurepath)
simresultsname <- paste0(simname,".RData")
simsprocessedname <- paste0(simname,".csv")

# CHANGE the path to where you downloaded the aghqmm package repo
# from https://github.com/awstringer1/aghqmm
aghqmmpath <- "~/work/projects/mixedmodel-computation/aghqmm" # CHANGE this

# Everything else in the script should run without changes

## Load Packages ##

pkgs <- c(
  'ggplot2',
  'dplyr',
  'tidyr',
  'readr',
  'lme4',
  'GLMMadaptive',
  'Rcpp',
  'RcppEigen',
  # 'parallel', # Disable parallel, for compatibility with Windows
  'fastmatrix',
  'remotes'
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
# numruns <- 10  # Number of times to execute the simulations
# numsims <- 100 # Number of simulations in each category PER RUN
# m <- c(100,200,500,1000)
# n <- c(3,5,7,9)
# k <- c(1,3,5,7,9,11,13,15,17,19,21,23,25)
# These settings are used for testing/continuous integration purposes:
numruns <- 2  # Number of times to execute the simulations
numsims <- 10 # Number of simulations in each category PER RUN
m <- c(100, 1000)
n <- c(3,5)
k <- c(1,3,5)

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
# mc.reset.stream() # Reproducbility in parallel- diabled for testing
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

cat("Summarizing simulations...\n")
simsummary <- simsprocessed |>
  mutate(sigmasq1relbias = (sigmasq1est-sigmasq1true)/sigmasq1true,
         sigmasq2relbias = (sigmasq2est-sigmasq2true)/sigmasq2true,
         sigmacovrelbias = (sigmacov1true-sigmacov1est)/sigmacov1true
  )


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
    scale_x_discrete(breaks = seq(min(k),max(k),by=4))
  # theme(axis.text.x = element_text(angle=45,size=7))
  if (!is.null(plt)) plt <- plt + coord_cartesian(ylim=ylim)
  plt
}


beta0biasplot <- biasboxplot("beta1",ylim = c(-10,10)) +
  labs(title = bquote(Bias*", "*widehat(beta)[0] - beta[0]),y = bquote(Bias(widehat(beta)[0])))

beta1biasplot <- biasboxplot("beta2",ylim = c(-10,10)) +
  labs(title = bquote(Bias*", "*widehat(beta)[1] - beta[1]),y = bquote(Bias(widehat(beta)[1])))

beta2biasplot <- biasboxplot("beta3",ylim = c(-10,10)) +
  labs(title = bquote(Bias*", "*widehat(beta)[2] - beta[2]),y = bquote(Bias(widehat(beta)[2])))

beta3biasplot <- biasboxplot("beta4",ylim = c(-10,10)) +
  labs(title = bquote(Bias*", "*widehat(beta)[3] - beta[3]),y = bquote(Bias(widehat(beta)[3])))

sigmasq1relbiasplot <- biasboxplot("sigmasq1rel",ylim = c(-1,250)) +
  labs(title = bquote(Relative~Bias*", "*widehat(sigma)[1]^2/sigma[1]^2),y = bquote(Relative~Bias(widehat(sigma)[1]^2)))
sigmasq1relbiasplotzoom <- biasboxplot("sigmasq1rel",ylim = c(-1,10)) +
  labs(title = bquote(Relative~Bias*", "*widehat(sigma)[1]^2/sigma[1]^2),y = bquote(Relative~Bias(widehat(sigma)[1]^2)))

sigmasq2relbiasplot <- biasboxplot("sigmasq2rel",ylim = c(-1,250)) +
  labs(title = bquote(Relative~Bias*", "*widehat(sigma)[2]^2/sigma[2]^2),y = bquote(Relative~Bias(widehat(sigma)[2]^2)))
sigmasq2relbiasplotzoom <- biasboxplot("sigmasq2rel",ylim = c(-1,10)) +
  labs(title = bquote(Relative~Bias*", "*widehat(sigma)[2]^2/sigma[2]^2),y = bquote(Relative~Bias(widehat(sigma)[2]^2)))

sigmacovbiasplot <- biasboxplot("sigmacov1",ylim = c(-1,250)) +
  labs(title = bquote(Bias*", "*widehat(sigma)[12]-sigma[12]),y = bquote(Bias(widehat(sigma)[12])))
sigmacovbiasplotzoom <- biasboxplot("sigmacov1",ylim = c(-1,10)) +
  labs(title = bquote(Bias*", "*widehat(sigma)[12]-sigma[12]),y = bquote(Bias(widehat(sigma)[12])))

# Save
ggsave(filename = file.path(figurepath,"beta0-bias.pdf"),plot=beta0biasplot,width=7,height=7)
ggsave(filename = file.path(figurepath,"beta1-bias.pdf"),plot=beta1biasplot,width=7,height=7)
ggsave(filename = file.path(figurepath,"beta2-bias.pdf"),plot=beta2biasplot,width=7,height=7)
ggsave(filename = file.path(figurepath,"beta3-bias.pdf"),plot=beta3biasplot,width=7,height=7)
ggsave(filename = file.path(figurepath,"sigmasq1-relbias.pdf"),plot=sigmasq1relbiasplot,width=7,height=7)
ggsave(filename = file.path(figurepath,"sigmasq1-relbias-zoom.pdf"),plot=sigmasq1relbiasplotzoom,width=7,height=7)
ggsave(filename = file.path(figurepath,"sigmasq2-relbias.pdf"),plot=sigmasq2relbiasplot,width=7,height=7)
ggsave(filename = file.path(figurepath,"sigmasq2-relbias-zoom.pdf"),plot=sigmasq2relbiasplotzoom,width=7,height=7)
ggsave(filename = file.path(figurepath,"sigmacov1-bias.pdf"),plot=sigmacovbiasplot,width=7,height=7)
ggsave(filename = file.path(figurepath,"sigmacov1-bias-zoom.pdf"),plot=sigmacovbiasplotzoom,width=7,height=7)


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
    scale_x_continuous(breaks = seq(min(k), max(k), by=4)) +
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

beta2covrplot <- coverageplot("beta3") + 
  labs(title = bquote(Empirical~Coverage*","~beta[2]))

beta3covrplot <- coverageplot("beta4") + 
  labs(title = bquote(Empirical~Coverage*","~beta[3]))

sigmasq1covrplot <- coverageplot("sigmasq1") + 
  labs(title = bquote(Empirical~Coverage*","~sigma[1]^2))

sigmasq2covrplot <- coverageplot("sigmasq2") + 
  labs(title = bquote(Empirical~Coverage*","~sigma[2]^2))

sigmacovcovrplot <- coverageplot("sigmacov1") + 
  labs(title = bquote(Empirical~Coverage*","~sigma[12]))

ggsave(filename = file.path(figurepath,"beta0-covr.pdf"),plot=beta0covrplot,width=7,height=7)
ggsave(filename = file.path(figurepath,"beta1-covr.pdf"),plot=beta1covrplot,width=7,height=7)
ggsave(filename = file.path(figurepath,"beta2-covr.pdf"),plot=beta2covrplot,width=7,height=7)
ggsave(filename = file.path(figurepath,"beta3-covr.pdf"),plot=beta3covrplot,width=7,height=7)
ggsave(filename = file.path(figurepath,"sigmasq1-covr.pdf"),plot=sigmasq1covrplot,width=7,height=7)
ggsave(filename = file.path(figurepath,"sigmasq2-covr.pdf"),plot=sigmasq2covrplot,width=7,height=7)
ggsave(filename = file.path(figurepath,"sigmacov1-covr.pdf"),plot=sigmacovcovrplot,width=7,height=7)


## Lengths
lengthboxplot <- function(var,data = simsummary,ylim=NULL) {
  vr <- sym(paste0(var,'lengthWald'))
  plt <- ggplot(data,aes(x=as.factor(k),y=!!vr)) +
    facet_grid(n~m,labeller = labeller(m = m_labeller,n = n_labeller)) +
    theme_bw() +
    geom_boxplot(outlier.size = .7,outlier.alpha = .05) +
    labs(x = "Number of Quadrature Points",y = "Upper Limit - Lower Limit") +
    scale_x_discrete(breaks = seq(min(k), max(k), by=4))
  if (!is.null(plt)) plt <- plt + coord_cartesian(ylim=ylim)
  plt
}

beta0lengthplot <- lengthboxplot("beta1",ylim = c(0,10)) + 
  labs(title = bquote(Length~of~Confidence~Intervals*","~beta[0]))

beta1lengthplot <- lengthboxplot("beta2",ylim = c(0,10)) + 
  labs(title = bquote(Length~of~Confidence~Intervals*","~beta[1]))

beta2lengthplot <- lengthboxplot("beta3",ylim = c(0,10)) + 
  labs(title = bquote(Length~of~Confidence~Intervals*","~beta[2]))

beta3lengthplot <- lengthboxplot("beta4",ylim = c(0,10)) + 
  labs(title = bquote(Length~of~Confidence~Intervals*","~beta[3]))

sigmasq1lengthplot <- lengthboxplot("sigmasq1",ylim = c(0,1000)) + 
  labs(title = bquote(Length~of~Confidence~Intervals*","~sigma[1]^2))
sigmasq1lengthplotzoom <- lengthboxplot("sigmasq1",ylim = c(0,10)) + 
  labs(title = bquote(Length~of~Confidence~Intervals*","~sigma[1]^2))

sigmasq2lengthplot <- lengthboxplot("sigmasq2",ylim = c(0,1000)) + 
  labs(title = bquote(Length~of~Confidence~Intervals*","~sigma[2]^2))
sigmasq2lengthplotzoom <- lengthboxplot("sigmasq2",ylim = c(0,10)) + 
  labs(title = bquote(Length~of~Confidence~Intervals*","~sigma[2]^2))

sigmacovlengthplot <- lengthboxplot("sigmacov1",ylim = c(0,1000)) + 
  labs(title = bquote(Length~of~Confidence~Intervals*","~sigma[12]))
sigmacovlengthplotzoom <- lengthboxplot("sigmacov1",ylim = c(0,10)) + 
  labs(title = bquote(Length~of~Confidence~Intervals*","~sigma[12]))

ggsave(filename = file.path(figurepath,"beta0-length.pdf"),plot=beta0lengthplot,width=7,height=7)
ggsave(filename = file.path(figurepath,"beta1-length.pdf"),plot=beta1lengthplot,width=7,height=7)
ggsave(filename = file.path(figurepath,"beta2-length.pdf"),plot=beta2lengthplot,width=7,height=7)
ggsave(filename = file.path(figurepath,"beta3-length.pdf"),plot=beta3lengthplot,width=7,height=7)
ggsave(filename = file.path(figurepath,"sigmasq1-length.pdf"),plot=sigmasq1lengthplot,width=7,height=7)
ggsave(filename = file.path(figurepath,"sigmasq1-length-zoom.pdf"),plot=sigmasq1lengthplotzoom,width=7,height=7)
ggsave(filename = file.path(figurepath,"sigmasq2-length.pdf"),plot=sigmasq2lengthplot,width=7,height=7)
ggsave(filename = file.path(figurepath,"sigmasq2-length-zoom.pdf"),plot=sigmasq2lengthplotzoom,width=7,height=7)
ggsave(filename = file.path(figurepath,"sigmacov1-length.pdf"),plot=sigmacovlengthplot,width=7,height=7)
ggsave(filename = file.path(figurepath,"sigmacov1-length-zoom.pdf"),plot=sigmacovlengthplotzoom,width=7,height=7)



### Main manuscript
# Save the results to be shown in the main manuscript
# This is Figure 1.

MAINTEXTSIZE <- 23
if (1000 %in% m & 5 %in% n) {
  mainsummary <- simsummary %>%
    filter(m==1000,n==5) %>%
    dplyr::select(
      k,m,n,
      contains('beta1'),
      contains('sigmasq1'),
      contains('sigmacov1')
    )


  main_biasplot <- function(var,data = mainsummary,ylim=NULL) {
    vr <- sym(paste0(var,'bias'))
    plt <- ggplot(data,aes(x=as.factor(k),y=!!vr)) +
      theme_bw() +
      geom_boxplot(outlier.size = .7,outlier.alpha = .05) +
      labs(x = "Number of Quadrature Points") +
      scale_x_discrete(breaks = seq(min(k), max(k), by=4)) +
      theme(text = element_text(size=MAINTEXTSIZE))
    if (!is.null(plt)) plt <- plt + coord_cartesian(ylim=ylim)
    plt
  }

  main_coverageplot <- function(var,data = mainsummary,ylim=NULL) {
    vW <- sym(paste0(var,"covrWald"))
    plt <- data %>%
      group_by(m,n,k) %>%
      summarize(covrWald = mean(!!vW,na.rm=TRUE),
                covrWaldsd = sqrt(covrWald*(1-covrWald)/n())
      ) %>%
      ggplot(aes(x=k)) +
      theme_bw() +
      geom_line(aes(y=covrWald)) +
      geom_line(aes(y=covrWald-2*covrWaldsd),linetype='dashed') +
      geom_line(aes(y=covrWald+2*covrWaldsd),linetype='dashed') +
      geom_point(aes(y=covrWald)) +
      scale_x_continuous(breaks = seq(min(k), max(k), by=4)) +
      geom_hline(yintercept = .95,linetype = 'dotted') +
      scale_y_continuous(limits = c(0,1),breaks = seq(0,1,by=.1)) +
      labs(x = "Number of Quadrature Points",y="Empirical Coverage Proportion") +
      theme(text = element_text(size=MAINTEXTSIZE))
    if (!is.null(plt)) {
      plt <- plt + coord_cartesian(ylim=ylim)
    } else {
      plt <- plt + coord_cartesian(ylim=c(0,1))
    }
    plt
  }

  main_lengthboxplot <- function(var,data = mainsummary,ylim=NULL) {
    vr <- sym(paste0(var,'lengthWald'))
    plt <- ggplot(data,aes(x=as.factor(k),y=!!vr)) +
      theme_bw() +
      geom_boxplot(outlier.size = .7,outlier.alpha = .05) +
      labs(x = "Number of Quadrature Points",y = "Upper Limit - Lower Limit") +
      scale_x_discrete(breaks = seq(min(k), max(k) ,by=4)) +
      theme(text = element_text(size=MAINTEXTSIZE))
    if (!is.null(plt)) plt <- plt + coord_cartesian(ylim=ylim)
    plt
  }

  main_beta0_bias_plot <- main_biasplot("beta1",data = mainsummary,ylim = c(-5,5)) +
    labs(title = bquote(Bias*", "*widehat(beta)[0] - beta[0]),y = bquote(Bias(widehat(beta)[0])))
  main_sigmasq1_bias_plot <- main_biasplot("sigmasq1rel",data = mainsummary,ylim = c(-5,5)) +
    labs(title = bquote(Relative~Bias*", "*widehat(sigma)[1]^2/sigma[1]^2),y = bquote(Relative~Bias(widehat(sigma)[1]^2)))
  main_sigmacovbiasplot <- main_biasplot("sigmacov1",ylim = c(-5,5)) +
    labs(title = bquote(Bias*", "*widehat(sigma)[12]-sigma[12]),y = bquote(Bias(widehat(sigma)[12])))


  main_beta0_covrplot <- main_coverageplot("beta1") + 
    labs(title = bquote(Empirical~Coverage*","~beta[0]))
  main_sigmasq1covrplot <- main_coverageplot("sigmasq1") + 
    labs(title = bquote(Empirical~Coverage*","~sigma[1]^2))
  main_sigmacovcovrplot <- main_coverageplot("sigmacov1") + 
    labs(title = bquote(Empirical~Coverage*","~sigma[12]))

  main_beta0lengthplot <- main_lengthboxplot("beta1",ylim = c(0,10)) + 
    labs(title = bquote(Length~of~Confidence~Intervals*","~beta[0]))
  main_sigmasq1lengthplot <- main_lengthboxplot("sigmasq1",ylim = c(0,10)) + 
    labs(title = bquote(Length~of~Confidence~Intervals*","~sigma[1]^2))
  main_sigmacovlengthplot <- main_lengthboxplot("sigmacov1",ylim = c(0,10)) + 
    labs(title = bquote(Length~of~Confidence~Intervals*","~sigma[12]))

  ggsave(filename = file.path(figurepath,"beta0-bias-main.pdf"),plot=main_beta0_bias_plot,width=7,height=7)
  ggsave(filename = file.path(figurepath,"sigmasq1-bias-main.pdf"),plot=main_sigmasq1_bias_plot,width=7,height=7)
  ggsave(filename = file.path(figurepath,"sigmacov1-bias-main.pdf"),plot=main_sigmacovbiasplot,width=7,height=7)

  ggsave(filename = file.path(figurepath,"beta0-covr-main.pdf"),plot=main_beta0_covrplot,width=7,height=7)
  ggsave(filename = file.path(figurepath,"sigmasq1-covr-main.pdf"),plot=main_sigmasq1covrplot,width=7,height=7)
  ggsave(filename = file.path(figurepath,"sigmacov1-covr-main.pdf"),plot=main_sigmacovcovrplot,width=7,height=7)

  ggsave(filename = file.path(figurepath,"beta0-length-main.pdf"),plot=main_beta0lengthplot,width=7,height=7)
  ggsave(filename = file.path(figurepath,"sigmasq1-length-main.pdf"),plot=main_sigmasq1lengthplot,width=7,height=7)
  ggsave(filename = file.path(figurepath,"sigmacov1-length-main.pdf"),plot=main_sigmacovlengthplot,width=7,height=7)
}

