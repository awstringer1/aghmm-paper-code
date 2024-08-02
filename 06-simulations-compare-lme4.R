# AGHQ using the C++ L-BFGS-B implementation
# through the aghqmm package
# Comparison to lme4 for a random-intercepts model

## Set paths ##
# CHANGE the base path to whatever you want on your machine
basepath <- getwd()
# CHANGE the name of the simulation to control how saved results are named
simname <- "sims-lme4-20240802-v1"
stopifnot(dir.exists(basepath))
resultspath <- file.path(basepath,'results')
if (!dir.exists(resultspath)) dir.create(resultspath)
figurepath <- file.path(basepath,'figures')
if (!dir.exists(figurepath)) dir.create(figurepath)
simresultsname <- paste0(simname,".RData")
simsprocessedname <- paste0(simname,".csv")

# Everything else in the script should run without changes

## Load Packages ##

install.packages(aghqmmpath,repos=NULL,type="source")

pkgs <- c(
  'ggplot2',
  'dplyr',
  'tidyr',
  'readr',
  'fastmatrix',
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
# These settings were used for the paper:
# numruns <- 2  # Number of times to execute the simulations
# numsims <- 10 # Number of simulations in each category PER RUN
# m <- c(100,200,500,1000) # m,n chosen from main manuscript
# n <- c(3,5,7,9)
# k <- seq(1,25,by=2)
# These are used for continuous integration/testing:
numruns <- 1  # Number of times to execute the simulations
numsims <- 2 # Number of simulations in each category PER RUN
m <- c(100, 1000) # m,n chosen from main manuscript
n <- c(3,5)
k <- seq(1,5,by=2)

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

# options(mc.cores = parallel::detectCores())
# RNGkind("L'Ecuyer-CMRG") # For reproducibility with parallel

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

simsummary <- simsprocessed

m_labeller <- function(vl) paste0("m = ",vl)
n_labeller <- function(vl) paste0("n = ",vl)

## Barplot of successes ##
successplot <- simsummary %>%
  group_by(method,m,n,k) %>%
  summarize(num=n(),successful=mean(successful)) %>%
  mutate(se = sqrt(successful*(1-successful)/num),
         lower = successful - 2*se,
         upper = successful + 2*se
  ) %>%
  # ggplot(aes(x = factor(k),y = successful,fill = method)) +
  ggplot(aes(x = factor(k),y = successful,group = method,colour = method)) +
  theme_bw() +
  facet_grid(n~m,labeller = labeller(m = m_labeller,n = n_labeller)) +
  geom_line() +
  geom_line(aes(y = lower),linetype='dashed') +
  geom_line(aes(y = upper),linetype='dashed') +
  # geom_bar(stat = "identity",position = position_dodge(),colour="black") +
  # geom_errorbar(aes(ymin=lower,ymax=upper),position=position_dodge(width=.9),size=1,width=.5) +
  # scale_fill_grey(start=.5,end=.8,labels = c("aghqmm" = "New Approach","lme4" = "lme4")) +
  scale_colour_grey(start=.5,end=.8,labels = c("aghqmm" = "New Approach","lme4" = "lme4")) +
  labs(
    title = "Successful Simulations",
    x = "Number of quadrature points",
    y = "Proportion of simulations that were successful",
    fill = "Method"
  ) +
  geom_hline(yintercept = 1,linetype='dotted') +
  scale_x_discrete(breaks = seq(1,25,by=4)) +
  scale_y_continuous(breaks = seq(0,1,by=.1),labels = scales::percent_format()) +
  coord_cartesian(ylim=c(0,1))

ggsave(filename = file.path(figurepath,"lme4-successplot.pdf"),plot = successplot,width=7,height=7)  

## Boxplots of Absolute Computation Time ##
# Show the absolute computation times; less relevant than relative, but necessary to know
absolute_timeplot <- simsummary %>%
  filter(successful == 1) %>%
  ggplot(aes(x = factor(k),y = comptime,fill = method)) +
  theme_bw() +
  facet_grid(n~m,labeller = labeller(m = m_labeller,n = n_labeller)) +
  geom_boxplot(outlier.size = .7,outlier.alpha = .05) +
  scale_fill_grey(start=.5,end=.8,labels = c("aghqmm" = "New Approach","lme4" = "lme4")) +
  labs(
    title = "Absolute Computation Times (seconds)",
    x = "Number of quadrature points",
    y = "Computation time (seconds)",
    fill = "Method"
  ) +
  scale_y_continuous(limits = c(0,2),breaks = seq(0,5,by=.25)) +
  scale_x_discrete(breaks = seq(1,25,by=4)) +
  coord_cartesian(ylim = c(0,2))

ggsave(filename = file.path(figurepath,"lme4-absolutetimeplot.pdf"),plot = absolute_timeplot,width=7,height=7)  

## Boxplots of Relative Computation Time ##
# More relevant to direct comparison of a method; less affected by hardware

aghq_times <- simsummary %>%
  filter(method == "aghqmm") %>%
  dplyr::select(idx,aghqtime = comptime)

reltimes <- simsummary %>%
  left_join(aghq_times,by='idx') %>%
  mutate(reltime = comptime / aghqtime)


relative_timeplot <- reltimes %>%
  filter(successful == 1,method == "lme4") %>%
  ggplot(aes(x = factor(k),y = reltime)) +
  theme_bw() +
  facet_grid(n~m,labeller = labeller(m = m_labeller,n = n_labeller)) +
  geom_boxplot(outlier.size = .7,outlier.alpha = .05) +
  labs(
    title = "Relative Computation Times, lme4 / New Approach",
    x = "Number of quadrature points",
    y = "Relative computation time, lme4 / New Approach"
  ) +
  scale_y_continuous(breaks = seq(-1,10,by=1)) +
  coord_cartesian(ylim = c(-1,10)) +
  scale_x_discrete(breaks = seq(1,25,by=4)) +
  geom_hline(yintercept=1,lty='dashed')

ggsave(filename = file.path(figurepath,"lme4-relativetimeplot.pdf"),plot = relative_timeplot,width=7,height=7)  

## Boxplots of log-likelihood values ##
# Do the actual nll, by method

raw_nll_plot <- simsummary %>%
  filter(successful == 1) %>%
  ggplot(aes(x = factor(k),y = nllavgbase10,fill = method)) +
  theme_bw() +
  facet_grid(n~m,labeller = labeller(m = m_labeller,n = n_labeller)) +
  geom_boxplot(outlier.size = .7,outlier.alpha = .05) +
  scale_fill_grey(start=.5,end=.8,labels = c("aghqmm" = "New Approach","lme4" = "lme4")) +
  labs(
    title = "Average negative log-likelihood values by method",
    x = "Number of quadrature points",
    y = bquote("-log"[10]~"likelihood / (mn)"),
    fill = "Method"
  ) +
  scale_y_continuous(breaks = seq(.1,.2,by=.05)) +
  scale_x_discrete(breaks = seq(1,25,by=4)) +
  coord_cartesian(ylim = c(.1,.2))

ggsave(filename = file.path(figurepath,"lme4-rawnllplot.pdf"),plot = raw_nll_plot,width=7,height=7)  

## Boxplots of differences in log-likelihood values ##
# Difference of nll for the same datasets by method

aghq_nll <- simsummary %>%
  filter(method == "aghqmm") %>%
  dplyr::select(idx,aghqnll = nllavgbase10)

nlldiffdata <- simsummary %>%
  left_join(aghq_nll,by='idx') %>%
  mutate(nlldiff = nllavgbase10 - aghqnll)

nlldiffplot <- nlldiffdata %>%
  filter(successful == 1,method == "lme4") %>%
  ggplot(aes(x = factor(k),y = nlldiff)) +
  theme_bw() +
  facet_grid(n~m,labeller = labeller(m = m_labeller,n = n_labeller)) +
  geom_boxplot(outlier.size = .7,outlier.alpha = .05) +
  labs(
    title = "Difference in negative log-likelihood, lme4 - New Approach",
    x = "Number of quadrature points",
    y = "Difference in negative log-likelihood, lme4 - New Approach"
  ) +
  # coord_cartesian(ylim = c(0,16)) +
  geom_hline(yintercept=0,lty='dashed') +
  scale_y_continuous(breaks = seq(-2,.2,by=.05)) +
  scale_x_discrete(breaks = seq(1,25,by=4)) +
  coord_cartesian(ylim = c(-.2,.2))


ggsave(filename = file.path(figurepath,"lme4-nlldiffplot.pdf"),plot = nlldiffplot,width=7,height=7)  



## Boxplots of gradient norm values ##
# Actual values, by method

raw_norm_plot <- simsummary %>%
  filter(successful == 1) %>%
  ggplot(aes(x = factor(k),y = normgradavg2base10log,fill = method)) +
  theme_bw() +
  facet_grid(n~m,labeller = labeller(m = m_labeller,n = n_labeller)) +
  geom_boxplot(outlier.size = .7,outlier.alpha = .05) +
  scale_fill_grey(start=.5,end=.8,labels = c("aghqmm" = "New Approach","lme4" = "lme4")) +
  labs(
    title = "Average log-norm of gradient by method",
    x = "Number of quadrature points",
    y = bquote("log"[10]~"||gradient||"[2]~"/ (mn)"),
    fill = "Method"
  ) +
  scale_y_continuous(breaks = seq(-20,-2,by=2)) +
  scale_x_discrete(breaks = seq(1,25,by=4)) +
  coord_cartesian(ylim = c(-20,-2))

ggsave(filename = file.path(figurepath,"lme4-rawnormplot.pdf"),plot = raw_norm_plot,width=7,height=7)  

## Boxplots of differences in log-norm of gradient values ##
# Difference of log gradient norms for the same datasets by method

aghq_norms <- simsummary %>%
  filter(method == "aghqmm") %>%
  dplyr::select(idx,aghqnorm = normgradavg2base10log)

normdiffdata <- simsummary %>%
  left_join(aghq_norms,by='idx') %>%
  mutate(normdiff = normgradavg2base10log - aghqnorm)

normdiffplot <- normdiffdata %>%
  filter(successful == 1,method == "lme4") %>%
  ggplot(aes(x = factor(k),y = normdiff)) +
  theme_bw() +
  facet_grid(n~m,labeller = labeller(m = m_labeller,n = n_labeller)) +
  geom_boxplot(outlier.size = .7,outlier.alpha = .05) +
  labs(
    title = "Difference in log-norm of gradient, lme4 - New Approach",
    x = "Number of quadrature points",
    y = "Difference in log-norm of gradient, lme4 - New Approach"
  ) +
  # coord_cartesian(ylim = c(0,16)) +
  geom_hline(yintercept=0,lty='dashed') +
  scale_y_continuous(breaks = seq(-8,12,by=2)) +
  scale_x_discrete(breaks = seq(1,25,by=4)) +
  coord_cartesian(ylim = c(-8,12))

ggsave(filename = file.path(figurepath,"lme4-normdiffplot.pdf"),plot = normdiffplot,width=7,height=7)  

### Main manuscript
# Save the results to be shown in the main manuscript

MAINTEXTSIZE <- 23

if (1000 %in% m & 5 %in% n) {
  main_reltimes <- reltimes %>%
    filter(m==1000,n==5,successful == 1,method == "lme4")

  main_nlldiffdata <- nlldiffdata %>%
    filter(m==1000,n==5,successful == 1,method == "lme4")

  main_normdiffdata <- normdiffdata %>%
    filter(m==1000,n==5,successful == 1,method == "lme4")


  relative_timeplot_manuscript <- main_reltimes %>%
    ggplot(aes(x = factor(k),y = reltime)) +
    theme_bw() +
    geom_boxplot(outlier.size = .7,outlier.alpha = .05) +
    labs(
      title = "Rel. Comp. Times, lme4",
      x = "Number of quadrature points",
      y = "Comp. Time, lme4 / New Approach"
    ) +
    scale_y_continuous(breaks = seq(0,6,by=1)) +
    coord_cartesian(ylim = c(0,6)) +
    scale_x_discrete(breaks = seq(1,25,by=4)) +
    geom_hline(yintercept=1,lty='dashed') +
    theme(text = element_text(size=MAINTEXTSIZE))

  ggsave(filename = file.path(figurepath,"lme4-relativetimeplot_manuscript.pdf"),plot = relative_timeplot_manuscript,width=7,height=7)


  nlldiffplot_manuscript <- main_nlldiffdata %>%
    ggplot(aes(x = factor(k),y = nlldiff)) +
    theme_bw() +
    geom_boxplot(outlier.size = .7,outlier.alpha = .05) +
    labs(
      title = "Diff. in neg. log-lik, lme4",
      x = "Number of quadrature points",
      y = "neg. log-lik, lme4 - New Approach"
    ) +
    geom_hline(yintercept=0,lty='dashed') +
    scale_y_continuous(breaks = seq(-.1,.1,by=.02)) +
    scale_x_discrete(breaks = seq(1,25,by=4)) +
    coord_cartesian(ylim = c(-.1,.1)) +
    theme(text = element_text(size=MAINTEXTSIZE))


  ggsave(filename = file.path(figurepath,"lme4-nlldiffplot_manuscript.pdf"),plot = nlldiffplot_manuscript,width=7,height=7)

  normdiffplot_manuscript <- main_normdiffdata %>%
    ggplot(aes(x = factor(k),y = normdiff)) +
    theme_bw() +
    geom_boxplot(outlier.size = .7,outlier.alpha = .05) +
    labs(
      title = "Diff. in log-norm of grad, lme4",
      x = "Number of quadrature points",
      y = "log-norm(grad), lme4 - New Approach"
    ) +
    geom_hline(yintercept=0,lty='dashed') +
    scale_y_continuous(breaks = seq(-4,12,by=2)) +
    scale_x_discrete(breaks = seq(1,25,by=4)) +
    coord_cartesian(ylim = c(-4,12)) +
    theme(text = element_text(size=MAINTEXTSIZE))

  ggsave(filename = file.path(figurepath,"lme4-normdiffplot_manuscript.pdf"),plot = normdiffplot_manuscript,width=7,height=7)
}
