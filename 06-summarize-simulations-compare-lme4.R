### Summarize sumulation results for supplement
### Alex Stringer
### 2023/07

library(tidyverse)

# CHANGE this to wherever you saved the simulation results
basepath <- '~/work/projects/mixedmodel-computation/replication'
resultspath <- file.path(basepath,'results')
stopifnot(dir.exists(resultspath)) # This should have been created in the previous script
figurepath <- file.path(basepath,'figures')
if (!dir.exists(figurepath)) dir.create(figurepath)
# CHANGE this to whatever you named the simulations
simname <- "sims-lme4-20230817-v1"
simresultsname <- paste0(simname,".RData")
simsummaryname <- paste0(simname,".csv")

simsummary <- readr::read_csv(file.path(resultspath,simsummaryname))

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








