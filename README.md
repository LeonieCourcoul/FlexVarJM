
<!-- README.md is generated from README.Rmd. Please edit that file -->

# FlexVarJM

<!-- badges: start -->
<!-- badges: end -->

The goal of VarJM is to estimate joint model with subject-specific
time-dependent variability.

The global function is FlexVar_JM. It handles to estimate joint model
with a marker which has a subject-specific time-dependent variability
and competing events with the possibility to take into account the left
truncation.

## Installation

You can install the development version of FlexVarJM from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("LeonieCourcoul/FlexVarJM")
```

## Exemple

# Estimation

This is an exemple in a simulated dataset. We estimated the following
model with two competing events k=1 or k=2 :

$$y_i(t_{ij}) = \color{blue}\tilde{y_i}(t_{ij}) \color{black} + \epsilon_{ij} = \beta_0 + b_{0i} + (\beta_1 + b_{1i})t_{ij} +\epsilon_{ij} $$
$$y_i(t_{ij})= \color{blue}\tilde{y}_i(t_{ij}) \color{black} + \epsilon_{ij} = \beta_0 + b_{0i} + (\beta_1 + b_{1i})t_{ij}+\epsilon_{ij}$$

$$ \left\{
\begin{array}{l}
y_i(t_{i j})= \color{blue}\tilde{y}_i(t_{ij}) \color{black} + \epsilon_{ij} = \beta_0 + b_{0i} + (\beta_1 + b_{1i})t_{ij}+\epsilon_{i j}\\
\lambda_{ik}(t)=\lambda_{0k}(t) \exp(\color{blue}\alpha_{1k}\tilde{y}_i(t)\color{black} + \color{red}\alpha_{\sigma k} \sigma_i(t) \color{black})
\end{array}
\right. $$

where :

- $\epsilon_{i}(t_{ij}) \sim \mathcal{N}(0, \color{red}\sigma_i^2\color{black})$
  with
  $\color{red}\log(\sigma_i(t_{ij})) = \mu_0 + \tau_{0i} + (\mu_1 + \tau_{1i})\times t_{ij}$

- with $b_i=\left(b_{0i},b_{1i}\right)^{\top}$ and
  $\tau_i=\left(\tau_{0i},\tau_{1i}\right)^{\top}$ assuming that the two
  sets of random effects $b_i$ and $\tau_i$ are not independent:
  $$\quad\left(\begin{array}{c}
  b_{i} \\
  \tau_i
  \end{array}\right) \sim N(0, \Sigma)$$

with the following Cholesky decomposition for the covariance matrix of
the random effects: $$\Sigma = \left(\begin{array}{cccc}
s_0 & 0 & 0 & 0  \\
s_{01} & s_1 & 0 & 0 \\
s_{02} & s_{12} & s_2 & 0 \\
s_{03} & s_{13} & s_{23} & s_3 
\end{array}\right)\left(\begin{array}{cccc}
s_0 & 0 & 0 & 0  \\
s_{01} & s_1 & 0 & 0 \\
s_{02} & s_{12} & s_2 & 0 \\
s_{03} & s_{13} & s_{23} & s_3 
\end{array}\right)^\top = LL^\top$$

- $\lambda_{0k}(t) = shape_k^2 t^{shape_k^2-1}e^{\zeta_{0k}}$ : Weibull
  function

``` r
exemple <- FlexVar_JM(formFixed = y~visit,
                      formRandom = ~ visit,
                      formGroup = ~ID,
                      formSurv = Surv(time, event ==1 ) ~ 1,
                      timeVar = "visit",
                      data.long = Data_toy,
                      variability_hetero = TRUE,
                      sharedtype = "CV",
                      hazard_baseline = "Weibull",
                      competing_risk = TRUE,
                      formSurv_CR = Surv(time, event ==2 ) ~ 1,
                      hazard_baseline_CR = "Weibull",
                      sharedtype_CR = "CV",
                      formFixedVar =~visit, 
                      formRandomVar =~visit,
                      correlated_re = TRUE,
                      S1 = 1000,
                      S2 = 8000,
                      nproc = 5
                      )
                      
summary.FlexVarJM(exemple)
```

You can access to the table of estimations and standard deviation with :

``` r
exemple$table.res
```

The computing time is given by :

``` r
exemple$time.compute
```

The output of the marqLevAlg algorithm is in :

``` r
exemple$result
```

Finally, some elements of control are in :

``` r
exemple$control
```

# Goodness-of-fit

You can check the goodness-of-fit of the longitudinal submodel and of
the survival submodel by computing the predicted random effects :

``` r
goodness <- goodness_of_fit(exemple, graph = T)
```

# Predictions

You can compute the probability for a (new) individual to have event 1
or 2 between time s and time S+t years given that he did not experience
any event before time s, its trajectory of marker until time s ans the
set of estimated parameters. For exemple, for the individual 15 to
experiment the event 1 between 3 and 5 years :

``` r
newdata <- Data_toy[which(Data_toy$ID == 15),]
predictions <- pred_s.t(newdata, exemple, s = 3, window = 1, event = 1, tirage = NULL)
```

To have a confidence interval, you can compute the predictions L times
(for exemple L = 500) by drawing the parameter in the normal
distribution of parameter the estimates and the covariance matrix. Then
you can take the 2.5th and 97.5th percentiles of the predictions :

``` r
  Hess <- matrix(rep(0,length(exemple$result$grad)**2),nrow=length(exemple$result$grad),ncol=length(exemple$result$grad))
  Hess[upper.tri(Hess, diag=T)] <- exemple$result$v
  Hess2 = Hess + t(Hess)
  diag(Hess2) <- diag(Hess2) - diag(Hess)
  predictions.boot <- c()
  for(k in 1:500){
    tirage <- rmvnorm(1, mean = estimation$table.res$Estimation, sigma = Hess2)
    predictions.boot <- c(predictions.boot, pred_s.t.ponctuel.tps(newdata = newdata, estimation, s = 3, window = 2, event = 1, tirage = tirage))
  }
```
