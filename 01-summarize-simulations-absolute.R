### Summarize sumulation results for supplement
### Alex Stringer
### 2023/07
### You must run the script '01-simulations-absolute.R' BEFORE running this one.

library(tidyverse)

# CHANGE this to wherever you saved the simulation results
basepath <- '~/work/projects/mixedmodel-computation/replication'
resultspath <- file.path(basepath,'results')
stopifnot(dir.exists(resultspath)) # This should have been created in the previous script
figurepath <- file.path(basepath,'figures')
if (!dir.exists(figurepath)) dir.create(figurepath)
# CHANGE this to whatever you named the simulations
simname <- "sims-20230817-v1"
simresultsname <- paste0(simname,".RData")
simsummaryname <- paste0(simname,".csv")

simsummary <- readr::read_csv(file.path(resultspath,simsummaryname)) %>%
  mutate(sigmasq1relbias = (sigmasq1est-sigmasq1true)/sigmasq1true,
         sigmasq2relbias = (sigmasq2est-sigmasq2true)/sigmasq2true,
         sigmacovrelbias = (sigmacov1true-sigmacov1est)/sigmacov1true
  )

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
    scale_x_discrete(breaks = seq(1,25,by=4))
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
    scale_x_discrete(breaks = seq(1,25,by=4)) +
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
    scale_x_continuous(breaks = seq(1,25,by=4)) +
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
    scale_x_discrete(breaks = seq(1,25,by=4)) +
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


