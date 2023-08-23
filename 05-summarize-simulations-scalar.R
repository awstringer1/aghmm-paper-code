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
simname <- "sims-scalar-20230817-v1"
simresultsname <- paste0(simname,".RData")
simsummaryname <- paste0(simname,".csv")

simsummary <- readr::read_csv(file.path(resultspath,simsummaryname)) %>%
  mutate(sigmasqrelbias = (sigmasqest-sigmasqtrue)/sigmasqtrue)

successes <- simsummary %>%
  group_by(m,n,k) %>%
  summarize(successful = n())
arrange(successes,successful)

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

