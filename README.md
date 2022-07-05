
<!-- README.md is generated from README.Rmd. Please edit that file -->

# FlexVarJM

<!-- badges: start -->
<!-- badges: end -->

The goal of VarJM is to estimate joint model with subject-specific
variability.

The global function is FlexVar_JM. It handles to estimate joint model
with a marker which has a subject-specific variability and competing
events with the possibility to take into account the left truncation.

## Installation

You can install the development version of FlexVarJM from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("LeonieCourcoul/FlexVarJM")
```

## Exemple

This is an exemple in a simulated dataset. We estimated the following
model with two competing events k=1 or k=2 :

![\\left\\{
\\begin{array}{ll}
y_i(t\_{i j})= \\color{blue}\\tilde{y}\_i(t\_{ij}) \\color{black} + \\epsilon\_{ij} = \\beta_0 + b\_{0i} + (\\beta_1 + b\_{1i})t\_{ij}+\\epsilon\_{i j}\\\\
\\lambda\_{ik}(t)=\\lambda\_{0k}(t) \\exp(\\color{blue}\\alpha\_{1k}\\tilde{y}\_i(t)\\color{black} + \\color{red}\\alpha\_{\\sigma k} \\sigma_i \\color{black})
\\end{array}
\\right.](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cleft%5C%7B%0A%5Cbegin%7Barray%7D%7Bll%7D%0Ay_i%28t_%7Bi%20j%7D%29%3D%20%5Ccolor%7Bblue%7D%5Ctilde%7By%7D_i%28t_%7Bij%7D%29%20%5Ccolor%7Bblack%7D%20%2B%20%5Cepsilon_%7Bij%7D%20%3D%20%5Cbeta_0%20%2B%20b_%7B0i%7D%20%2B%20%28%5Cbeta_1%20%2B%20b_%7B1i%7D%29t_%7Bij%7D%2B%5Cepsilon_%7Bi%20j%7D%5C%5C%0A%5Clambda_%7Bik%7D%28t%29%3D%5Clambda_%7B0k%7D%28t%29%20%5Cexp%28%5Ccolor%7Bblue%7D%5Calpha_%7B1k%7D%5Ctilde%7By%7D_i%28t%29%5Ccolor%7Bblack%7D%20%2B%20%5Ccolor%7Bred%7D%5Calpha_%7B%5Csigma%20k%7D%20%5Csigma_i%20%5Ccolor%7Bblack%7D%29%0A%5Cend%7Barray%7D%0A%5Cright. "\left\{
\begin{array}{ll}
y_i(t_{i j})= \color{blue}\tilde{y}_i(t_{ij}) \color{black} + \epsilon_{ij} = \beta_0 + b_{0i} + (\beta_1 + b_{1i})t_{ij}+\epsilon_{i j}\\
\lambda_{ik}(t)=\lambda_{0k}(t) \exp(\color{blue}\alpha_{1k}\tilde{y}_i(t)\color{black} + \color{red}\alpha_{\sigma k} \sigma_i \color{black})
\end{array}
\right.")

where :

-   ![\\epsilon\_{ij} \\sim \\mathcal{N}(0, \\color{red}\\sigma_i^2\\color{black})](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cepsilon_%7Bij%7D%20%5Csim%20%5Cmathcal%7BN%7D%280%2C%20%5Ccolor%7Bred%7D%5Csigma_i%5E2%5Ccolor%7Bblack%7D%29 "\epsilon_{ij} \sim \mathcal{N}(0, \color{red}\sigma_i^2\color{black})")
    with
    ![\\color{red}\\log(\\sigma_i) \\sim \\mathcal{N}(\\mu\_\\sigma, \\tau\_\\sigma^2)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ccolor%7Bred%7D%5Clog%28%5Csigma_i%29%20%5Csim%20%5Cmathcal%7BN%7D%28%5Cmu_%5Csigma%2C%20%5Ctau_%5Csigma%5E2%29 "\color{red}\log(\sigma_i) \sim \mathcal{N}(\mu_\sigma, \tau_\sigma^2)")

-   ![b_i \\sim \\mathcal{N}(0, \\Sigma_b)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;b_i%20%5Csim%20%5Cmathcal%7BN%7D%280%2C%20%5CSigma_b%29 "b_i \sim \mathcal{N}(0, \Sigma_b)")
    with
    ![\\begin{pmatrix} s\_{0} & 0\\\\ s\_{01} & s_1 \\end{pmatrix}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cbegin%7Bpmatrix%7D%20s_%7B0%7D%20%26%200%5C%5C%20s_%7B01%7D%20%26%20s_1%20%5Cend%7Bpmatrix%7D "\begin{pmatrix} s_{0} & 0\\ s_{01} & s_1 \end{pmatrix}")
    the decomposition cholesky of
    ![\\Sigma_b](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5CSigma_b "\Sigma_b")

-   ![\\lambda\_{0k}(t) = \\eta_k t^{\\eta_k-1}e^{\\zeta\_{0k}}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Clambda_%7B0k%7D%28t%29%20%3D%20%5Ceta_k%20t%5E%7B%5Ceta_k-1%7De%5E%7B%5Czeta_%7B0k%7D%7D "\lambda_{0k}(t) = \eta_k t^{\eta_k-1}e^{\zeta_{0k}}")
    : Weibull function

``` r
exemple <- FlexVar_JM(formFixed = y~visit,
                      formRandom = ~ visit,
                      formGroup = ~ID,
                      formSurv = Surv(time, event ==1 ) ~ 1,
                      timeVar = "visit",
                      nb.e.a = 2,
                      data.long = Data_exemple,
                      variability_hetero = TRUE,
                      sharedtype = "CV",
                      hazard_baseline = "Weibull",
                      competing_risk = TRUE,
                      formSurv_CR = Surv(time, event ==2 ) ~ 1,
                      hazard_baseline_CR = "Weibull",
                      sharedtype_CR = "CV",
                      S1 = 1000,
                      S2 = 8000,
                      nproc = 5
                      )
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
