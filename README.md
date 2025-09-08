# sprintr

`sprintr` is an R package for **fast, scalable interaction modeling** in **generalized linear models (GLMs)**. It implements the *sprinter* (sparse reluctant interaction) algorithm: fit a regularized **main-effects** GLM once, then score each candidate interaction by a **one-dimensional MLE** conditional on that fit, and finally refit on the selected set. This realizes **interaction reluctance**—prefer main effects when an interaction adds no clear extra predictive utility—**without** imposing hierarchy constraints. The full method and theory are developed in **Lu & Yu, *Reluctant Interaction Modeling in Generalized Linear Models***.&#x20;

## Why *sprintr*?

* **No hierarchy assumptions.** *sprintr* encodes **interaction reluctance**: it prioritizes main effects over interactions when their predictive gains are comparable, so it **does not require hierarchical constraints**. This extends the idea in **Yu, Bien & Tibshirani (2019), *Reluctant Interaction Modeling*** to GLMs.
* **Broad GLM coverage.** The framework applies widely to a broader range of GLMs (e.g., **logistic**, **Poisson**, **multinomial logistic**, **ordinal logistic**) under a unified, likelihood-based conditional screening framework.
* **Statistical guarantees.** Provides **finite-sample selection consistency** for recovering “pure” interactions (those adding signal beyond main effects) and **size control** of the selected interaction set under high-dimensional regimes (see theorems in the paper).&#x20;
* **Efficiency at scale.** Scores each interaction via a **1-D MLE** and selects top-m in a **single pass** over all $q$ candidates—scaling to **millions** (even **billions**) of interactions in practice, with strong empirical speedups while maintaining competitive accuracy.&#x20;

## Installation

```r
# install devtools if needed
install.packages("devtools")

# install the latest from GitHub
devtools::install_github("hluscience/sprintr")
```

## Quick start
```r
set.seed(123)
library(sprintr)

n <- 100
p <- 100

x  <- matrix(rnorm(n * p), n, p)
mu <- x[,1] - 2*x[,2] + 3*(x[,1]*x[,3]) - 4*(x[,4]*x[,5])

## -------------------------
## Gaussian regression
## -------------------------
y_gau  <- mu + rnorm(n)
mod_gau <- cv.sprinter(x = x, y = y_gau)

## -------------------------
## Binomial (logistic) regression
## -------------------------
prob   <- plogis(mu)
y_bin  <- rbinom(n = n, size = 1, prob = prob)
mod_bin <- cv.sprinter(x = x, y = y_bin,
                       family = "binomial",
                       type.measure = "deviance")
```

## Reproducible studies

All code to reproduce the paper’s simulation studies and application is available at:
**[https://github.com/hluscience/reproducible/tree/main/sprintr](https://github.com/hluscience/reproducible/tree/main/sprintr)**.&#x20;
