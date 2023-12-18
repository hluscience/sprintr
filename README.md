# sprintr
The `sprintr` package implements a computationally efficient method for fitting large-scale interaction modeling in generalized linear models (GLMs). It honors the reluctance principle from [Yu, Bien, and Tibshirani (2021) Reluctant Interaction Modeling](https://arxiv.org/pdf/1907.08414.pdf) and extends its application beyond Gaussian linear regression interaction modeling to a broader range of GLMs, including Logistic regression, Poisson regression, Multinomial logistic regression, and Ordinal logistic regression.

To install `sprintr` from GitHub, enter the following command in the R console:

```r
devtools::install_github("hluscience/sprintr")
```

Please note that this installation process requires the R package [devtools](https://cran.r-project.org/web/packages/devtools/index.html). If you don't have it installed, you can do so using the command: `install.packages("devtools")`
