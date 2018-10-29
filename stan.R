library(rstan)
library(ggplot2)
library(Cairo)

###############################################################################
#### http://mc-stan.org/users/documentation/case-studies/pool-binary-trials.html

sample_size <- c(47, 103, 316, 235, 46, (196-10+368), 249, 295, 88)
somatic <- c(3, 5, 19, 11, 4, 56, 20, 0, 2)

N = length(sample_size)

mydata <- list(K = sample_size, y = somatic, N= N)


## ------- code ----------------
no_pool_stan <- '
data {
  int<lower=0> N;           // items
  int<lower=0> K[N];        // initial trials
  int<lower=0> y[N];        // initial successes
}
parameters {
  real<lower=0, upper=1> phi;         // population chance of success
real<lower=1> kappa;                // population concentration
vector<lower=0, upper=1>[N] theta;  // chance of success 
}
model {
  kappa ~ pareto(1, 1.5);                        // hyperprior
  theta ~ beta(phi * kappa, (1 - phi) * kappa);  // prior
  y ~ binomial(K, theta);                        // likelihood
}
'
fit <- stan(model_code = no_pool_stan, data = mydata, iter = 1000, 
            chains = 3, control = list(adapt_delta = 0.99,
                                       max_treedepth = 10))
## ------------ theta --------------------
plot(fit, pars ='theta')

## ------------ kappa ---------------------
plot(fit, pars = 'phi')

## ------------ trace -----------------
traceplot(fit, pars = 'theta')

## ------------ diag ------------------
stan_diag(fit)

CairoPNG(file ='prevalence.png',
         width =1400, height = 700)
p <- stan_plot(fit, pars = 'theta')

p + scale_y_discrete(name ="Study", limits = names(fit)[3:12],
                     labels = c('Current study',
                                'Pujade, SOLO2/ENGOT-Ov21, 2017',
                                'Pennington, 2014',
                                'Coleman, ARIEL3, 2017',
                                'Chao, 2016',
                                'Hennessy, 2010', 'TCGA, 2011',
                                'McAlpine, 2011',
                                'Mafficini, 2016')) +
  theme(axis.text.y = element_text(face="bold", size=36))+
  theme(axis.text.x = element_text(face="bold", size=36))
dev.off()

p

