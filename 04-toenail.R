### Toenail data ###

## BEGIN SETUP ##

## Set paths
basepath <- getwd()
stopifnot(dir.exists(basepath))
resultspath <- file.path(basepath,'results')
if (!dir.exists(resultspath)) dir.create(resultspath)
figurepath <- file.path(basepath,'figures')
if (!dir.exists(figurepath)) dir.create(figurepath)

## Set parameters 
# These are the values used to reproduce the results in the manuscript:
# NUMBOOT <- 200
# ktodo <- seq(1,25,by=2)
# NUMRUNS <- 5

# These are the values used for continuous integration using Github Actions:
NUMBOOT <- 5
ktodo <- seq(1,5,by=2)
NUMRUNS <- 1 # For computation times

## Libraries ----

# Test for, and install if not found, the packages required to run this script.
pkgs <- c(
  # 'tidyverse',
  'dplyr',
  'ggplot2',
  'tidyr',
  'Matrix',
  'lme4',
  'GLMMadaptive',
  'mice',
  'remotes',
  'fastmatrix'
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

## Load data
data(toenail) # From package 'mice', loaded above

## END SETUP ##

## Fit the model ----

LMElist <- list()
length(LMElist) <- length(ktodo)
AQlist <- GAlist <- LMElist
aghq_control <- aghqmm_control(
  bfgsdelta = c(1e-01),
  inner_tol = c(1e-06),
  inner_maxitr = c(10),
  bfgshist = c(6),
  past = c(3),
  wolfe = .9,
  ftol = 1e-04,
  max_linesearch = 100,
  tol = 1e-04,
  maxitr = 100
)
for (j in 1:length(ktodo)) {
  k <- ktodo[j]
  cat("Fitting for k = ",k,"\n",sep="")
  cat("lme4...")
  LMElist[[j]] <- tryCatch(
    aghqmm::aghqmm(outcome ~ treatment*month + (1|ID),data=toenail,k=k,method="lme4"),
    error = function(e) e
  )
  cat(" done.\n")
  cat("AGHQ...")
  AQlist[[j]] <- tryCatch(
    aghqmm::aghqmm(outcome ~ treatment*month + (1|ID),data=toenail,k=k,method="both",control = aghq_control),
    error = function(e) e
  )
  cat(" done.\n")
  cat("GLMMa...")
  GAlist[[j]] <- tryCatch(
    aghqmm::aghqmm(outcome ~ treatment*month + (1|ID),data=toenail,k=k,method="GLMMadaptive",control = aghq_control),
    error = function(e) e
  )
  cat(" done.\n")
  LMElist[[j]]$k <- k
  AQlist[[j]]$k <- k
  GAlist[[j]]$k <- k
}

LMEerrs <- Reduce(c,Map(inherits,LMElist,what="condition"))
AQerrs <- Reduce(c,Map(inherits,AQlist,what="condition"))
GAerrs <- Reduce(c,Map(inherits,GAlist,what="condition"))

LMElist <- LMElist[!LMEerrs]
AQlist <- AQlist[!AQerrs]
GAlist <- GAlist[!GAerrs]

paramsummary_beta <- tibble(
  method = c(Reduce(c,Map("[[",LMElist,"method")),Reduce(c,Map("[[",AQlist,"method")),Reduce(c,Map("[[",GAlist,"method"))),
  k = c(Reduce(c,Map("[[",LMElist,"k")),Reduce(c,Map("[[",AQlist,"k")),Reduce(c,Map("[[",GAlist,"k"))),
  beta0_lower = c(Reduce(c,Map("[[",Map("[[",LMElist,"betaints"),1)),Reduce(c,Map("[[",Map("[[",AQlist,"betaints"),1)),Reduce(c,Map("[[",Map("[[",GAlist,"betaints"),1))),
  beta1_lower = c(Reduce(c,Map("[[",Map("[[",LMElist,"betaints"),2)),Reduce(c,Map("[[",Map("[[",AQlist,"betaints"),2)),Reduce(c,Map("[[",Map("[[",GAlist,"betaints"),2))),
  beta2_lower = c(Reduce(c,Map("[[",Map("[[",LMElist,"betaints"),3)),Reduce(c,Map("[[",Map("[[",AQlist,"betaints"),3)),Reduce(c,Map("[[",Map("[[",GAlist,"betaints"),3))),
  beta3_lower = c(Reduce(c,Map("[[",Map("[[",LMElist,"betaints"),4)),Reduce(c,Map("[[",Map("[[",AQlist,"betaints"),4)),Reduce(c,Map("[[",Map("[[",GAlist,"betaints"),4))),
  beta0_point = c(Reduce(c,Map("[[",Map("[[",LMElist,"betaints"),5)),Reduce(c,Map("[[",Map("[[",AQlist,"betaints"),5)),Reduce(c,Map("[[",Map("[[",GAlist,"betaints"),5))),
  beta1_point = c(Reduce(c,Map("[[",Map("[[",LMElist,"betaints"),6)),Reduce(c,Map("[[",Map("[[",AQlist,"betaints"),6)),Reduce(c,Map("[[",Map("[[",GAlist,"betaints"),6))),
  beta2_point = c(Reduce(c,Map("[[",Map("[[",LMElist,"betaints"),7)),Reduce(c,Map("[[",Map("[[",AQlist,"betaints"),7)),Reduce(c,Map("[[",Map("[[",GAlist,"betaints"),7))),
  beta3_point = c(Reduce(c,Map("[[",Map("[[",LMElist,"betaints"),8)),Reduce(c,Map("[[",Map("[[",AQlist,"betaints"),8)),Reduce(c,Map("[[",Map("[[",GAlist,"betaints"),8))),
  beta0_upper = c(Reduce(c,Map("[[",Map("[[",LMElist,"betaints"),9)),Reduce(c,Map("[[",Map("[[",AQlist,"betaints"),9)),Reduce(c,Map("[[",Map("[[",GAlist,"betaints"),9))),
  beta1_upper = c(Reduce(c,Map("[[",Map("[[",LMElist,"betaints"),10)),Reduce(c,Map("[[",Map("[[",AQlist,"betaints"),10)),Reduce(c,Map("[[",Map("[[",GAlist,"betaints"),10))),
  beta2_upper = c(Reduce(c,Map("[[",Map("[[",LMElist,"betaints"),11)),Reduce(c,Map("[[",Map("[[",AQlist,"betaints"),11)),Reduce(c,Map("[[",Map("[[",GAlist,"betaints"),11))),
  beta3_upper = c(Reduce(c,Map("[[",Map("[[",LMElist,"betaints"),12)),Reduce(c,Map("[[",Map("[[",AQlist,"betaints"),12)),Reduce(c,Map("[[",Map("[[",GAlist,"betaints"),12)))
) |>
  pivot_longer(
    beta0_lower:beta3_upper,
    names_to = c("param","type"),
    names_sep = "_",
    values_to = "value"
  ) |>
  pivot_wider(names_from = type,values_from=value)

# lme4 does not provide Wald intervals for sigma
# prepare the frame then add the profile and bootstrap intervals after, with their time comparisons
paramsummary_sigma <- data.frame(
  method = c(Reduce(c,Map("[[",GAlist,"method")),Reduce(c,Map("[[",AQlist,"method")),Reduce(c,Map("[[",LMElist,"method"))),
  k = c(Reduce(c,Map("[[",GAlist,"k")),Reduce(c,Map("[[",AQlist,"k")),Reduce(c,Map("[[",LMElist,"k"))),
  sigmasq_lower = c(Reduce(c,Map("[[",Map("[[",GAlist,"sigmaints"),1)),Reduce(c,Map("[[",Map("[[",AQlist,"sigmaints"),1)),rep(-1,length(LMElist))),
  sigmasq_point = c(Reduce(c,Map("[[",Map("[[",GAlist,"sigmaints"),2)),Reduce(c,Map("[[",Map("[[",AQlist,"sigmaints"),2)),rep(-1,length(LMElist))),
  sigmasq_upper = c(Reduce(c,Map("[[",Map("[[",GAlist,"sigmaints"),3)),Reduce(c,Map("[[",Map("[[",AQlist,"sigmaints"),3)),rep(-1,length(LMElist)))
)
lmethetas <- as.numeric(Reduce(c,Map("[",Map("[[",LMElist,"theta"),5)))
paramsummary_sigma[paramsummary_sigma$method == "lme4", ]$sigmasq_point <- exp(-lmethetas/2)
# now replicate the lme4 results and change the method to profile/boot
lme4results <- paramsummary_sigma |> filter(method == "lme4")
bootresults <- profileresults <- as.data.frame(lme4results)
bootresults[["method"]] <- rep("boot", nrow(bootresults))
profileresults[["method"]] <- rep("profile", nrow(bootresults))

# Now, do the bootstrapping and the profile
profilelist <- bootlist <- list()
length(profilelist) <- length(bootlist) <- length(ktodo)
profiletimes <- boottimes <- numeric(length(ktodo))
names(profiletimes) <- names(boottimes) <- ktodo

for (j in 1:length(ktodo)) {
  k <- ktodo[j]
  cat("Fitting for k = ",k,"\n",sep="")
  cat("model...")
  mod <- tryCatch(lme4::glmer(outcome ~ treatment*month + (1|ID),data=toenail,family=binomial(),nAGQ = k),error=function(e) e)
  if (inherits(mod,'condition')) {
    profilelist[[j]] <- mod
    bootlist[[j]] <- mod
    profiletimes[j] <- -1
    boottimes[j] <- -1
    next
  }
  cat("\nprofile...")
  tm <- Sys.time()
  tmp <- tryCatch(confint(mod,parm=".sig01",method="profile"),error = function(e) e)
  profilelist[[j]] <- tmp
  if (inherits(tmp,'condition')) {
    profiletimes[j] <- -1
  } else {
    profiletimes[j] <- as.numeric(difftime(Sys.time(),tm,units='secs'))
  }
  cat("\nboot...")
  tm <- Sys.time()
  tmp <- tryCatch(confint(mod,parm=".sig01",method="boot",nsim=NUMBOOT),error = function(e) e)
  bootlist[[j]] <- tmp
  if (inherits(tmp,'condition')) {
    boottimes[j] <- -1
  } else {
    boottimes[j] <- as.numeric(difftime(Sys.time(),tm,units='secs'))
  }
  cat("\n")
}
save(bootlist,profilelist,boottimes,profiletimes,file = file.path(resultspath,"smoking-profilebootresults.RData"))
# load(file.path(resultspath,"smoking-profilebootresults.RData"))

# append to results
get_bounds <- function(lst) {
  if (inherits(lst,'condition'))
    return(c(-1,-1))
  return(unname(as.numeric(lst)))
}

profilebounds <- as.matrix(Reduce(rbind, Map(get_bounds,profilelist)))
bootbounds <- as.matrix(Reduce(rbind, Map(get_bounds,bootlist)))

profileresults[ ,'sigmasq_lower'] <- as.numeric(profilebounds[ ,1])
profileresults[ ,'sigmasq_upper'] <- as.numeric(profilebounds[ ,2])
bootresults[ ,'sigmasq_lower'] <- as.numeric(bootbounds[ ,1])
bootresults[ ,'sigmasq_upper'] <- as.numeric(bootbounds[ ,2])

# convert the aghq/glmma results to sigma
paramsummary_sigma_all <- paramsummary_sigma |>
  filter(method != "lme4") |>
  mutate_at(vars(contains("sigma")),sqrt) |>
  bind_rows(profileresults) |>
  bind_rows(bootresults) |>
  rename(sigma_lower = sigmasq_lower,sigma_point = sigmasq_point,sigma_upper = sigmasq_upper)

## plots


# beta
DODGEWIDTH <- 1
YLIM <- c(-5,3)
YBREAKS <- seq(min(YLIM),max(YLIM),by=1)
XBREAKS <- seq(1,25,by=4)


paramsummary_beta$paramf <- factor(paramsummary_beta$param,labels = c("beta[0]","beta[1]","beta[2]","beta[3]"))

param_plot_beta <- paramsummary_beta |>
  ggplot(aes(y = point,x=factor(k),ymin = lower,ymax = upper,linetype=method)) +
  theme_bw() +
  facet_wrap(~paramf,labeller = label_parsed) +
  geom_linerange(position = position_dodge(width=DODGEWIDTH)) +
  geom_point(size = 2,position = position_dodge(width=DODGEWIDTH)) +
  scale_y_continuous(breaks = YBREAKS) +
  scale_x_discrete(breaks = XBREAKS) +
  scale_linetype_manual(
    values = c("both" = "solid","GLMMadaptive" = "dashed","lme4" = "dotted"),
    labels = c("both" = "New method","GLMMadaptive" = "GLMMadaptive","lme4" = "lme4")
  ) +
  coord_cartesian(ylim = YLIM) +
  # coord_flip() +
  labs(title = bquote("Point and interval estimates for"~beta~"parameters"),
       y = "Estimated value",
       x = "Number of quadrature points",
       linetype = "Method")

ggsave(filename = file.path(figurepath,"toenail-beta-supplement.pdf"),plot=param_plot_beta,width=7,height=7)

# for manuscript
MAINTEXTSIZE <- 23
POINTSIZE <- 1.5
pointest_beta0 <- paramsummary_beta |> filter(param=="beta0",k==25,method=="both") |> pull(point)
beta0_plot <- paramsummary_beta |>
  filter(param == "beta0") |>
  ggplot(aes(y = point,x=factor(k),ymin = lower,ymax = upper,linetype=method)) +
  theme_bw() +
  geom_linerange(position = position_dodge(width=DODGEWIDTH)) +
  geom_point(size = POINTSIZE,position = position_dodge(width=DODGEWIDTH)) +
  geom_hline(yintercept = pointest_beta0,linetype='dotted') +
  scale_y_continuous(breaks = seq(-6,0,by=0.5)) +
  scale_linetype_manual(
    values = c("both" = "solid","GLMMadaptive" = "dashed","lme4" = "dotted"),
    labels = c("both" = "New method","GLMMadaptive" = "GLMMadaptive","lme4" = "lme4")
  ) +  
  coord_cartesian(ylim = c(-6,0)) +
  # coord_flip() +
  labs(title = bquote("Point and interval estimates for"~beta[0]),
       y = "Estimated value",
       x = "Number of quadrature points",
       linetype = "Method") +
  guides(linetype="none") +
  theme(text = element_text(size = MAINTEXTSIZE))

ggsave(filename = file.path(figurepath,"toenail-beta0-manuscript.pdf"),plot=beta0_plot,width=7,height=7)

pointest_beta1 <- paramsummary_beta |> filter(param=="beta1",k==25,method=="both") |> pull(point)
beta1_plot <- paramsummary_beta |>
  filter(param == "beta1") |>
  ggplot(aes(y = point,x=factor(k),ymin = lower,ymax = upper,linetype=method)) +
  theme_bw() +
  geom_linerange(position = position_dodge(width=DODGEWIDTH)) +
  geom_point(size = POINTSIZE,position = position_dodge(width=DODGEWIDTH)) +
  geom_hline(yintercept = pointest_beta1,linetype='dotted') +
  scale_y_continuous(breaks = seq(-2,2,by=0.5)) +
  scale_linetype_manual(
    values = c("both" = "solid","GLMMadaptive" = "dashed","lme4" = "dotted"),
    labels = c("both" = "New method","GLMMadaptive" = "GLMMadaptive","lme4" = "lme4")
  ) +  
  coord_cartesian(ylim = c(-2,2)) +
  # coord_flip() +
  labs(title = bquote("Point and interval estimates for"~beta[1]),
       y = "Estimated value",
       x = "Number of quadrature points",
       linetype = "Method") +
  guides(linetype="none") +
  theme(text = element_text(size = MAINTEXTSIZE))

ggsave(filename = file.path(figurepath,"toenail-beta1-manuscript.pdf"),plot=beta1_plot,width=7,height=7)

pointest_sigma <- paramsummary_sigma_all |> filter(k==25,method=="both") |> pull(sigma_point)
sigma_plot <- paramsummary_sigma_all |>
  ggplot(aes(y = sigma_point,x=factor(k),ymin = sigma_lower,ymax = sigma_upper,linetype=method)) +
  theme_bw() +
  geom_linerange(position = position_dodge(width=DODGEWIDTH)) +
  geom_point(size = POINTSIZE,position = position_dodge(width=DODGEWIDTH)) +
  geom_hline(yintercept = pointest_sigma,linetype='dotted') +
  scale_y_continuous(breaks = seq(0,10,by=0.5)) +
  scale_linetype_manual(
    values = c("both" = "solid","GLMMadaptive" = "dashed","profile" = "dotted","boot" = "dotdash"),
    labels = c("both" = "New method","GLMMadaptive" = "GLMMadaptive","profile" = "Profile","boot" = "Bootstrap")
  ) +  
  coord_cartesian(ylim = c(0,10)) +
  # coord_flip() +
  labs(title = bquote("Point and interval estimates for"~sigma),
       y = "Estimated value",
       x = "Number of quadrature points",
       linetype = "Method") +
  guides(linetype="none") +
  theme(text = element_text(size = MAINTEXTSIZE))

ggsave(filename = file.path(figurepath,"toenail-sigma-manuscript.pdf"),plot=sigma_plot,width=7,height=7)

# sigma
DODGEWIDTH <- 1
YLIM <- c(0,12)
YBREAKS <- seq(min(YLIM),max(YLIM),by=2)
XBREAKS <- seq(1,25,by=4)


param_plot_sigma <- paramsummary_sigma_all |>
  ggplot(aes(y = sigma_point,x=factor(k),ymin = sigma_lower,ymax = sigma_upper,linetype=method)) +
  theme_bw() +
  geom_linerange(position = position_dodge(width=DODGEWIDTH)) +
  geom_point(size = 2,position = position_dodge(width=DODGEWIDTH)) +
  scale_y_continuous(breaks = YBREAKS) +
  scale_x_discrete(breaks = XBREAKS) +
  # scale_linetype_discrete(labels = c("both" = "New method","GLMMadaptive" = "GLMMadaptive")) +
  scale_linetype_manual(
    values = c("both" = "solid","GLMMadaptive" = "dashed","profile" = "dotted","boot" = "dotdash"),
    labels = c("both" = "New method","GLMMadaptive" = "GLMMadaptive","profile" = "Profile","boot" = "Bootstrap")
  ) +
  coord_cartesian(ylim = YLIM) +
  # coord_flip() +
  labs(title = bquote("Point and interval estimates for"~sigma^2),
       y = "Estimated value",
       x = "Number of quadrature points",
       linetype = "Method")

ggsave(filename = file.path(figurepath,"toenail-sigma-supplement.pdf"),plot=param_plot_sigma,width=7,height=7)

## computation time plots

compsummary <- tibble(
  method = c(Reduce(c,Map("[[",LMElist,"method")),Reduce(c,Map("[[",AQlist,"method")),Reduce(c,Map("[[",GAlist,"method"))),
  k = c(Reduce(c,Map("[[",LMElist,"k")),Reduce(c,Map("[[",AQlist,"k")),Reduce(c,Map("[[",GAlist,"k"))),
  time = c(Reduce(c,Map("[[",LMElist,"comptime")),Reduce(c,Map("[[",AQlist,"comptime")),Reduce(c,Map("[[",GAlist,"comptime"))),
  nll = c(Reduce(c,Map("[[",LMElist,"nll")),Reduce(c,Map("[[",AQlist,"nll")),Reduce(c,Map("[[",GAlist,"nll"))),
  normgrad = c(Reduce(c,Map("[[",LMElist,"normgrad_2")),Reduce(c,Map("[[",AQlist,"normgrad_2")),Reduce(c,Map("[[",GAlist,"normgrad_2")))
) |>
  mutate(
    nll = nll / (log(10)*nrow(toenail)), # convert to base 10
    lognormgrad = log(normgrad / nrow(toenail),base=10)
  )

# relative times
aghqtimes <- compsummary |>
  filter(method == "both") |>
  dplyr::select(k,aghqtime = time,aghqnll = nll,aghqlognormgrad = lognormgrad)

reltimes_lme4 <- compsummary |>
  filter(method == "lme4") |>
  left_join(aghqtimes,by='k') |>
  mutate(reltime = time / aghqtime,nlldiff = nll - aghqnll,lognormdiff = lognormgrad - aghqlognormgrad) |>
  dplyr::select(k,reltime,nlldiff,lognormdiff)
reltimes_GA <- compsummary |>
  filter(method == "GLMMadaptive") |>
  left_join(aghqtimes,by='k') |>
  mutate(reltime = time / aghqtime,nlldiff = nll - aghqnll,lognormdiff = lognormgrad - aghqlognormgrad) |>
  dplyr::select(k,reltime,nlldiff,lognormdiff)


# computation time
# re-run for a long time
LMEcomptimes <- GAcomptimes <- AQcomptimes <- list()
length(LMEcomptimes) <- length(GAcomptimes) <- length(AQcomptimes) <- length(ktodo)
tm <- Sys.time()
for (j in 1:length(ktodo)) {
  k <- ktodo[j]
  GAcomptimes[[j]]$k <- k
  GAcomptimes[[j]]$comptimes <- numeric(NUMRUNS)
  for (i in 1:NUMRUNS) {
    cat("k = ",k,", run = ",i," of ",NUMRUNS," for GLMMadaptive...","\n",sep="")
    tmp <- tryCatch(
      aghqmm::aghqmm(outcome ~ treatment*month + (1|ID),data=toenail,k=k,method="GLMMadaptive",control = aghq_control),
      error = function(e) e
    )
    if (inherits(tmp,'condition')) {
      GAcomptimes[[j]]$comptimes[i] <- -1
    } else {
      GAcomptimes[[j]]$comptimes[i] <- tmp$comptime
    }
  }
  cat("Done GLMMadaptive.\n")
  AQcomptimes[[j]]$k <- k
  AQcomptimes[[j]]$comptimes <- numeric(NUMRUNS)
  for (i in 1:NUMRUNS) {
    cat("k = ",k,", run = ",i," of ",NUMRUNS," for aghqmm...","\n",sep="")
    tmp <- tryCatch(
      aghqmm::aghqmm(outcome ~ treatment*month + (1|ID),data=toenail,k=k,method="both",control = aghq_control),
      error = function(e) e
    )
    if (inherits(tmp,'condition')) {
      AQcomptimes[[j]]$comptimes[i] <- -1
    } else {
      AQcomptimes[[j]]$comptimes[i] <- tmp$comptime
    }
  }
  cat("Done aghqmm.\n")
  LMEcomptimes[[j]]$k <- k
  LMEcomptimes[[j]]$comptimes <- numeric(NUMRUNS)
  for (i in 1:NUMRUNS) {
    cat("k = ",k,", run = ",i," of ",NUMRUNS," for lme4...","\n",sep="")
    tmp <- tryCatch(
      aghqmm::aghqmm(outcome ~ treatment*month + (1|ID),data=toenail,k=k,method="lme4"),
      error = function(e) e
    )
    if (inherits(tmp,'condition')) {
      LMEcomptimes[[j]]$comptimes[i] <- -1
    } else {
      LMEcomptimes[[j]]$comptimes[i] <- tmp$comptime
    }
  }
  cat("Done lme4.\n")
}
totalcomptime <- as.numeric(difftime(Sys.time(),tm,units='secs'))
cat("Computation times took",totalcomptime,"seconds.\n")

save(GAcomptimes,AQcomptimes,LMEcomptimes,file = file.path(resultspath,"toenail-comptimes.RData"))
load(file.path(resultspath,"toenail-comptimes.RData"))


for (j in 1:length(GAcomptimes)) {
  GAcomptimes[[j]]$rep <- 1:length(GAcomptimes[[j]]$comptimes)
  GAcomptimes[[j]]$method <- "GLMMadaptive"
}
for (j in 1:length(AQcomptimes)) {
  AQcomptimes[[j]]$rep <- 1:length(AQcomptimes[[j]]$comptimes)
  AQcomptimes[[j]]$method <- "AQ"
}
for (j in 1:length(AQcomptimes)) {
  LMEcomptimes[[j]]$rep <- 1:length(LMEcomptimes[[j]]$comptimes)
  LMEcomptimes[[j]]$method <- "lme4"
}

GAreltimeslist <- LMEreltimeslist <- AQcomptimes
for (j in 1:length(AQcomptimes)) {
  GAreltimeslist[[j]]$comptime <- GAcomptimes[[j]]$comptimes/AQcomptimes[[j]]$comptimes
  GAreltimeslist[[j]]$method <- "GArelative"
  LMEreltimeslist[[j]]$comptime <- LMEcomptimes[[j]]$comptimes/AQcomptimes[[j]]$comptimes
  LMEreltimeslist[[j]]$method <- "LMErelative"
}

process_row <- function(lst) {
  data.frame(
    method = lst$method,
    k = lst$k,
    comptime = lst$comptime,
    rep = lst$rep
  )
}

GAreltimesframe <- bind_rows(Map(process_row,GAreltimeslist))
LMEreltimesframe <- bind_rows(Map(process_row,LMEreltimeslist))
reltimesframe <- bind_rows(GAreltimesframe,LMEreltimesframe) |>
  filter(comptime > 0)


relcomptimesplot <- reltimesframe |>
  ggplot(aes(x = factor(k),y = comptime,fill = method)) +
  theme_bw() +
  geom_boxplot(outlier.size = .7,outlier.alpha = .05) +
  labs(
    title = "Relative computation times",
    x = "Number of quadrature points",
    y = "Rel. comp. time, Method / New Approach",
    fill = "Method"
  ) +
  scale_fill_grey(start=.5,end=.8,labels = c("LMErelative" = "lme4","GArelative" = "GLMMadaptive")) +
  scale_y_continuous(breaks = seq(0,8,by=1)) +
  coord_cartesian(ylim = c(0,8)) +
  scale_x_discrete(breaks = seq(1,25,by=4)) +
  geom_hline(yintercept=1,lty='dashed')

relcomptimesplotmain <- relcomptimesplot +
  theme(text = element_text(size = MAINTEXTSIZE))

ggsave(filename = file.path(figurepath,"toenail-relativecomputation-supplement.pdf"),plot=relcomptimesplot,width=7,height=7)
ggsave(filename = file.path(figurepath,"toenail-relativecomputation-manuscript.pdf"),plot=relcomptimesplotmain,width=7,height=7)

# relative time of the single boot/profile fit
profiletimes / aghqtimes$aghqtime
boottimes / aghqtimes$aghqtime


# absolute computation times

GAcomptimesframe <- bind_rows(Map(process_row,GAcomptimes))
AQcomptimesframe <- bind_rows(Map(process_row,AQcomptimes))
LMEcomptimesframe <- bind_rows(Map(process_row,LMEcomptimes))
comptimesframe <- bind_rows(GAcomptimesframe,AQcomptimesframe,LMEcomptimesframe)

abscomptimesplot <- comptimesframe |>
  filter(comptime > 0) |>
  ggplot(aes(x = factor(k),y = comptime,fill = method)) +
  theme_bw() +
  geom_boxplot(outlier.size = .7,outlier.alpha = .05) +
  labs(
    title = "Absolute computation times (seconds), toenail fungus treatment data",
    x = "Number of quadrature points",
    y = "Computation time (seconds)",
    fill = "Method"
  ) +
  scale_y_continuous(breaks = seq(0,1,by=.1)) +
  scale_fill_grey(start=.4,end=.9,labels = c("AQ" = "New method","GLMMadaptive" = "GLMMadaptive","lme4" = "lme4")) +
  coord_cartesian(ylim = c(0,1)) +
  scale_x_discrete(breaks = seq(1,25,by=4))

ggsave(filename = file.path(figurepath,"toenail-abscomptimes-supplement.pdf"),plot=abscomptimesplot,width=7,height=7)
