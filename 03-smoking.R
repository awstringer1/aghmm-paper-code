### Data Analysis 1: Smoking Cessation
### Alex Stringer
### 2023/07

## Set paths
# CHANGE the base path to whatever you want on your machine
basepath <- getwd()
stopifnot(dir.exists(basepath))
resultspath <- file.path(basepath,'results')
if (!dir.exists(resultspath)) dir.create(resultspath)
# CHANGE the data path to wherever you saved the data on your machine
datapath <- file.path(basepath,"data")
dataname <- "SmkStudy.dat"
stopifnot(dir.exists(datapath))
stopifnot(file.exists(file.path(datapath,dataname))) # Data
figurepath <- file.path(basepath,'figures')
if (!dir.exists(figurepath)) dir.create(figurepath)

## Load external packages
pkgs <- c(
  'ggplot2',
  'dplyr',
  'tidyr',
  'readr',
  'fastmatrix',
  'lme4',
  'GLMMadaptive',
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


## Load data

smoking <- read_table(file.path(datapath,dataname),col_names = FALSE)
colnames(smoking) <- c("id","quit","time","group","groupxtime")

## Fit models
# Number of times to re-reun whole procedure for assessing computation time
# This is what's used in the paper:
# NUMRUNS <- 500
# This is what I use for testing the code:
NUMRUNS <- 5


# Fit using GLMMadaptive and aghqmm, for multiple k

ktodo <- seq(1,25,by=2)
GAlist <- list()
length(GAlist) <- length(ktodo)
AQlist <- GAlist
for (j in 1:length(ktodo)) {
  k <- ktodo[j]
  GAlist[[j]] <- tryCatch(
    aghqmm::aghqmm(quit ~ time*group + (time|id),data=smoking,k=k,method="GLMMadaptive"),
    error = function(e) e
  )
  AQlist[[j]] <- tryCatch(
    aghqmm::aghqmm(quit ~ time*group + (time|id),data=smoking,k=k,method="both"),
    error = function(e) e
  )
  GAlist[[j]]$k <- k
  AQlist[[j]]$k <- k
}
GAerrs <- Reduce(c,Map(inherits,GAlist,what="condition"))
AQerrs <- Reduce(c,Map(inherits,AQlist,what="condition"))

GAlist <- GAlist[!GAerrs]
AQlist <- AQlist[!AQerrs]


## summaries
# computation time, norm(gradient), nll,
# beta and sigma intervals
compsummary <- tibble(
  method = c(Reduce(c,Map("[[",GAlist,"method")),Reduce(c,Map("[[",AQlist,"method"))),
  k = c(Reduce(c,Map("[[",GAlist,"k")),Reduce(c,Map("[[",AQlist,"k"))),
  time = c(Reduce(c,Map("[[",GAlist,"comptime")),Reduce(c,Map("[[",AQlist,"comptime"))),
  nll = c(Reduce(c,Map("[[",GAlist,"nll")),Reduce(c,Map("[[",AQlist,"nll"))),
  normgrad = c(Reduce(c,Map("[[",GAlist,"normgrad_2")),Reduce(c,Map("[[",AQlist,"normgrad_2")))
) %>%
  mutate(
    nll = nll / (log(10)*nrow(smoking)), # convert to base 10
    lognormgrad = log(normgrad / nrow(smoking),base=10)
  )

# relative times
aghqtimes <- compsummary %>%
  filter(method == "both") %>%
  dplyr::select(k,aghqtime = time,aghqnll = nll,aghqlognormgrad = lognormgrad)

reltimes <- compsummary %>%
  filter(method == "GLMMadaptive") %>%
  left_join(aghqtimes,by='k') %>%
  mutate(reltime = time / aghqtime,nlldiff = nll - aghqnll,lognormdiff = lognormgrad - aghqlognormgrad) %>%
  dplyr::select(k,reltime,nlldiff,lognormdiff)

# point estimates and confidence intervals
# data frame with rows for parameter/method/k/type (lower/point/higher)
paramsummary <- tibble(
  method = c(Reduce(c,Map("[[",GAlist,"method")),Reduce(c,Map("[[",AQlist,"method"))),
  k = c(Reduce(c,Map("[[",GAlist,"k")),Reduce(c,Map("[[",AQlist,"k"))),
  beta0_lower = c(Reduce(c,Map("[[",Map("[[",GAlist,"betaints"),1)),Reduce(c,Map("[[",Map("[[",AQlist,"betaints"),1))),
  beta1_lower = c(Reduce(c,Map("[[",Map("[[",GAlist,"betaints"),2)),Reduce(c,Map("[[",Map("[[",AQlist,"betaints"),2))),
  beta2_lower = c(Reduce(c,Map("[[",Map("[[",GAlist,"betaints"),3)),Reduce(c,Map("[[",Map("[[",AQlist,"betaints"),3))),
  beta3_lower = c(Reduce(c,Map("[[",Map("[[",GAlist,"betaints"),4)),Reduce(c,Map("[[",Map("[[",AQlist,"betaints"),4))),
  beta0_point = c(Reduce(c,Map("[[",Map("[[",GAlist,"betaints"),5)),Reduce(c,Map("[[",Map("[[",AQlist,"betaints"),5))),
  beta1_point = c(Reduce(c,Map("[[",Map("[[",GAlist,"betaints"),6)),Reduce(c,Map("[[",Map("[[",AQlist,"betaints"),6))),
  beta2_point = c(Reduce(c,Map("[[",Map("[[",GAlist,"betaints"),7)),Reduce(c,Map("[[",Map("[[",AQlist,"betaints"),7))),
  beta3_point = c(Reduce(c,Map("[[",Map("[[",GAlist,"betaints"),8)),Reduce(c,Map("[[",Map("[[",AQlist,"betaints"),8))),
  beta0_upper = c(Reduce(c,Map("[[",Map("[[",GAlist,"betaints"),9)),Reduce(c,Map("[[",Map("[[",AQlist,"betaints"),9))),
  beta1_upper = c(Reduce(c,Map("[[",Map("[[",GAlist,"betaints"),10)),Reduce(c,Map("[[",Map("[[",AQlist,"betaints"),10))),
  beta2_upper = c(Reduce(c,Map("[[",Map("[[",GAlist,"betaints"),11)),Reduce(c,Map("[[",Map("[[",AQlist,"betaints"),11))),
  beta3_upper = c(Reduce(c,Map("[[",Map("[[",GAlist,"betaints"),12)),Reduce(c,Map("[[",Map("[[",AQlist,"betaints"),12))),
  sigma1sq_lower = c(Reduce(c,Map("[[",Map("[[",GAlist,"sigmaints"),1)),Reduce(c,Map("[[",Map("[[",AQlist,"sigmaints"),1))),
  sigma2sq_lower = c(Reduce(c,Map("[[",Map("[[",GAlist,"sigmaints"),2)),Reduce(c,Map("[[",Map("[[",AQlist,"sigmaints"),2))),
  sigma12_lower = c(Reduce(c,Map("[[",Map("[[",GAlist,"sigmaints"),3)),Reduce(c,Map("[[",Map("[[",AQlist,"sigmaints"),3))),
  sigma1sq_point = c(Reduce(c,Map("[[",Map("[[",GAlist,"sigmaints"),4)),Reduce(c,Map("[[",Map("[[",AQlist,"sigmaints"),4))),
  sigma2sq_point = c(Reduce(c,Map("[[",Map("[[",GAlist,"sigmaints"),5)),Reduce(c,Map("[[",Map("[[",AQlist,"sigmaints"),5))),
  sigma12_point = c(Reduce(c,Map("[[",Map("[[",GAlist,"sigmaints"),6)),Reduce(c,Map("[[",Map("[[",AQlist,"sigmaints"),6))),
  sigma1sq_upper = c(Reduce(c,Map("[[",Map("[[",GAlist,"sigmaints"),7)),Reduce(c,Map("[[",Map("[[",AQlist,"sigmaints"),7))),
  sigma2sq_upper = c(Reduce(c,Map("[[",Map("[[",GAlist,"sigmaints"),8)),Reduce(c,Map("[[",Map("[[",AQlist,"sigmaints"),8))),
  sigma12_upper = c(Reduce(c,Map("[[",Map("[[",GAlist,"sigmaints"),9)),Reduce(c,Map("[[",Map("[[",AQlist,"sigmaints"),9)))
) %>%
  pivot_longer(
    beta0_lower:sigma12_upper,
    names_to = c("param","type"),
    names_sep = "_",
    values_to = "value"
  ) %>%
  pivot_wider(names_from = type,values_from=value)


# beta
DODGEWIDTH <- 1
YLIM <- c(-5,3)
YBREAKS <- seq(min(YLIM),max(YLIM),by=1)
XBREAKS <- seq(1,25,by=4)

parambetadata <- paramsummary %>%
  filter(grepl('beta',param))
parambetadata$paramf <- factor(parambetadata$param,labels = c("beta[0]","beta[1]","beta[2]","beta[3]"))

param_plot_beta <- parambetadata %>%
  ggplot(aes(y = point,x=factor(k),ymin = lower,ymax = upper,linetype=method)) +
  theme_bw() +
  facet_wrap(~paramf,labeller = label_parsed) +
  geom_linerange(position = position_dodge(width=DODGEWIDTH)) +
  geom_point(size = 2,position = position_dodge(width=DODGEWIDTH)) +
  scale_y_continuous(breaks = YBREAKS) +
  scale_x_discrete(breaks = XBREAKS) +
  scale_linetype_discrete(labels = c("both" = "New method","GLMMadaptive" = "GLMMadaptive")) +
  coord_cartesian(ylim = YLIM) +
  # coord_flip() +
  labs(title = bquote("Point and interval estimates for"~beta~"parameters"),
       y = "Estimated value",
       x = "Number of quadrature points",
       linetype = "Method")

ggsave(filename = file.path(figurepath,"smoking-beta-supplement.pdf"),plot=param_plot_beta,width=7,height=7)

# sigma
DODGEWIDTH <- 1
YLIM <- c(-1,12)
YBREAKS <- seq(min(YLIM),max(YLIM),by=2)
XBREAKS <- seq(1,25,by=4)

paramsigmadata <- paramsummary %>%
  filter(grepl('sigma',param))
paramsigmadata$paramf <- factor(paramsigmadata$param,labels = c("sigma[12]","sigma[1]^2","sigma[2]^2"))

param_plot_sigma <- paramsigmadata %>%
  ggplot(aes(y = point,x=factor(k),ymin = lower,ymax = upper,linetype=method)) +
  theme_bw() +
  facet_wrap(~paramf,labeller = label_parsed,ncol=2,nrow=2) +
  geom_linerange(position = position_dodge(width=DODGEWIDTH)) +
  geom_point(size = 2,position = position_dodge(width=DODGEWIDTH)) +
  scale_y_continuous(breaks = YBREAKS) +
  scale_x_discrete(breaks = XBREAKS) +
  scale_linetype_discrete(labels = c("both" = "New method","GLMMadaptive" = "GLMMadaptive")) +
  coord_cartesian(ylim = YLIM) +
  # coord_flip() +
  labs(title = bquote("Point and interval estimates for"~sigma~"parameters"),
       y = "Estimated value",
       x = "Number of quadrature points",
       linetype = "Method")

ggsave(filename = file.path(figurepath,"smoking-sigma-supplement.pdf"),plot=param_plot_sigma,width=7,height=7)

## plots for manuscript
MAINTEXTSIZE <- 23
POINTSIZE <- 1.5

# beta0
pointest_beta0 <- paramsummary %>% filter(param=="beta0",k==25,method=="both") %>% pull(point)
beta0_plot <- paramsummary %>%
  filter(param == "beta0") %>%
  ggplot(aes(y = point,x=factor(k),ymin = lower,ymax = upper,linetype=method)) +
  theme_bw() +
  geom_linerange(position = position_dodge(width=DODGEWIDTH)) +
  geom_point(size = POINTSIZE,position = position_dodge(width=DODGEWIDTH)) +
  geom_hline(yintercept = pointest_beta0,linetype='dotted') +
  scale_y_continuous(breaks = seq(-5,-1,by=0.5)) +
  scale_linetype_discrete(labels = c("both" = "New method","GLMMadaptive" = "GLMMadaptive")) +
  coord_cartesian(ylim = c(-5,-1)) +
  # coord_flip() +
  labs(title = bquote("Point and interval estimates for"~beta[0]),
       y = "Estimated value",
       x = "Number of quadrature points",
       linetype = "Method") +
  guides(linetype="none") +
  theme(text = element_text(size = MAINTEXTSIZE))

ggsave(filename = file.path(figurepath,"smoking-beta0-manuscript.pdf"),plot=beta0_plot,width=7,height=7)



# sigmasq1
pointest_sigma1 <- paramsummary %>% filter(param=="sigma1sq",k==25,method=="both") %>% pull(point)
sigmasq1_plot <- paramsummary %>%
  filter(param == "sigma1sq") %>%
  ggplot(aes(y = point,x=factor(k),ymin = lower,ymax = upper,linetype=method)) +
  theme_bw() +
  geom_linerange(position = position_dodge(width=DODGEWIDTH)) +
  geom_point(size = POINTSIZE,position = position_dodge(width=DODGEWIDTH)) +
  geom_hline(yintercept = pointest_sigma1,linetype='dotted') +
  scale_y_continuous(breaks = seq(0,12,by=1)) +
  scale_linetype_discrete(labels = c("both" = "New method","GLMMadaptive" = "GLMMadaptive")) +
  coord_cartesian(ylim = c(0,12)) +
  # coord_flip() +
  labs(title = bquote("Point and interval estimates for"~sigma[1]^2),
       y = "Estimated value",
       x = "Number of quadrature points",
       linetype = "Method") +
  guides(linetype="none") +
  theme(text = element_text(size = MAINTEXTSIZE))


ggsave(filename = file.path(figurepath,"smoking-sigma1sq-manuscript.pdf"),plot=sigmasq1_plot,width=7,height=7)



# sigma12
pointest_sigma12 <- paramsummary %>% filter(param=="sigma12",k==25,method=="both") %>% pull(point)
sigma12_plot <- paramsummary %>%
  filter(param == "sigma12") %>%
  ggplot(aes(y = point,x=factor(k),ymin = lower,ymax = upper,linetype=method)) +
  theme_bw() +
  geom_linerange(position = position_dodge(width=DODGEWIDTH)) +
  geom_point(size = POINTSIZE,position = position_dodge(width=DODGEWIDTH)) +
  geom_hline(yintercept = pointest_sigma12,linetype='dotted') +
  scale_y_continuous(breaks = seq(-2,2,by=.5)) +
  scale_linetype_discrete(labels = c("both" = "New method","GLMMadaptive" = "GLMMadaptive")) +
  coord_cartesian(ylim = c(-2,2)) +
  # coord_flip() +
  labs(title = bquote("Point and interval estimates for"~sigma[12]),
       y = "Estimated value",
       x = "Number of quadrature points",
       linetype = "Method") +
  guides(linetype="none") +
  theme(text = element_text(size = MAINTEXTSIZE))


ggsave(filename = file.path(figurepath,"smoking-sigma12-manuscript.pdf"),plot=sigma12_plot,width=7,height=7)

# sigma2sq
pointest_sigma2sq <- paramsummary %>% filter(param=="sigma2sq",k==25,method=="both") %>% pull(point)
sigma2sq_plot <- paramsummary %>%
  filter(param == "sigma2sq") %>%
  ggplot(aes(y = point,x=factor(k),ymin = lower,ymax = upper,linetype=method)) +
  theme_bw() +
  geom_linerange(position = position_dodge(width=DODGEWIDTH)) +
  geom_point(size = POINTSIZE,position = position_dodge(width=DODGEWIDTH)) +
  geom_hline(yintercept = pointest_sigma2sq,linetype='dotted') +
  scale_y_continuous(breaks = seq(0,6.5,by=.5)) +
  scale_linetype_discrete(labels = c("both" = "New method","GLMMadaptive" = "GLMMadaptive")) +
  coord_cartesian(ylim = c(0,6.5)) +
  # coord_flip() +
  labs(title = bquote("Point and interval estimates for"~sigma[2]^2),
       y = "Estimated value",
       x = "Number of quadrature points",
       linetype = "Method") +
  guides(linetype="none") +
  theme(text = element_text(size = MAINTEXTSIZE))

ggsave(filename = file.path(figurepath,"smoking-sigma2sq-manuscript.pdf"),plot=sigma2sq_plot,width=7,height=7)


# computation time
# re-run for a long time
GAcomptimes <- AQcomptimes <- list()
length(GAcomptimes) <- length(AQcomptimes) <- length(ktodo)
tm <- Sys.time()
for (j in 1:length(ktodo)) {
  k <- ktodo[j]
  GAcomptimes[[j]]$k <- k
  GAcomptimes[[j]]$comptimes <- numeric(NUMRUNS)
  for (i in 1:NUMRUNS) {
    cat("k = ",k,", run = ",i," of ",NUMRUNS," for GLMMadaptive...","\n",sep="")
    tmp <- tryCatch(
      aghqmm::aghqmm(quit ~ time*group + (time|id),data=smoking,k=k,method="GLMMadaptive"),
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
      aghqmm::aghqmm(quit ~ time*group + (time|id),data=smoking,k=k,method="both"),
      error = function(e) e
    )
    if (inherits(tmp,'condition')) {
      AQcomptimes[[j]]$comptimes[i] <- -1
    } else {
      AQcomptimes[[j]]$comptimes[i] <- tmp$comptime
    }
  }
  cat("Done aghqmm.\n")
}
totalcomptime <- as.numeric(difftime(Sys.time(),tm,units='secs'))
cat("Computation times took",totalcomptime,"seconds.\n")

save(GAcomptimes,AQcomptimes,file = file.path(resultspath,"smoking-comptimes.RData"))
load(file.path(resultspath,"smoking-comptimes.RData"))

for (j in 1:length(GAcomptimes)) {
  GAcomptimes[[j]]$rep <- 1:length(GAcomptimes[[j]]$comptimes)
  GAcomptimes[[j]]$method <- "GLMMadaptive"
}
for (j in 1:length(AQcomptimes)) {
  AQcomptimes[[j]]$rep <- 1:length(AQcomptimes[[j]]$comptimes)
  AQcomptimes[[j]]$method <- "AQ"
}

reltimeslist <- AQcomptimes
for (j in 1:length(AQcomptimes)) {
  reltimeslist[[j]]$comptime <- GAcomptimes[[j]]$comptimes/AQcomptimes[[j]]$comptimes
  reltimeslist[[j]]$method <- "relative"
}

# post-process more
process_row <- function(lst) {
  data.frame(
    method = lst$method,
    k = lst$k,
    comptime = lst$comptime,
    rep = lst$rep
  )
}

GAcomptimesframe <- bind_rows(Map(process_row,GAcomptimes))
AQcomptimesframe <- bind_rows(Map(process_row,AQcomptimes))
relcomptimesframe <- bind_rows(Map(process_row,reltimeslist))

relcomptimesplot <- relcomptimesframe %>%
  ggplot(aes(x = factor(k),y = comptime)) +
  theme_bw() +
  geom_boxplot(outlier.size = .7,outlier.alpha = .05) +
  labs(
    title = "Relative computation times",
    x = "Number of quadrature points",
    y = "Rel. comp. time, GLMMa / New Approach"
  ) +
  scale_y_continuous(breaks = seq(0,10,by=1)) +
  coord_cartesian(ylim = c(0,10)) +
  scale_x_discrete(breaks = seq(1,25,by=4)) +
  geom_hline(yintercept=1,lty='dashed')

relcomptimesplotmain <- relcomptimesplot +
  theme(text = element_text(size = MAINTEXTSIZE))

ggsave(filename = file.path(figurepath,"smoking-relcomptimes-supplement.pdf"),plot=relcomptimesplot,width=7,height=7)
ggsave(filename = file.path(figurepath,"smoking-relcomptimes-manuscript.pdf"),plot=relcomptimesplotmain,width=7,height=7)

abscomptimesplot <- GAcomptimesframe %>%
  bind_rows(AQcomptimesframe) %>%
  filter(comptime > 0) %>%
  ggplot(aes(x = factor(k),y = comptime,fill = method)) +
  theme_bw() +
  geom_boxplot(outlier.size = .7,outlier.alpha = .05) +
  labs(
    title = "Absolute computation times (seconds), smoking cessation data",
    x = "Number of quadrature points",
    y = "Computation time (seconds)",
    fill = "Method"
  ) +
  scale_y_continuous(breaks = seq(0,35,by=5)) +
  scale_fill_grey(start=.5,end=.8,labels = c("AQ" = "New method","GLMMadaptive" = "GLMMadaptive")) +
  coord_cartesian(ylim = c(0,35)) +
  scale_x_discrete(breaks = seq(1,25,by=4))

ggsave(filename = file.path(figurepath,"smoking-abscomptimes-supplement.pdf"),plot=abscomptimesplot,width=7,height=7)





