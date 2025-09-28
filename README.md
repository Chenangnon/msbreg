# msbreg
Multistage Binomial (MSB) Modeling in R

This repository provides supplementary material to the paper: ``A multistage binomial model with measurement errors: application to protein viability prediction``

MSB regression is suited for modelling binomial dependent variables resulting from multistep data generation processes: the success probability of the target event results from many individual stages, each modelled with a separate linear predictor. This assumes multiplicative absolute risk between stages, and the common assumption of multiplicative odds (between independent variables) within each stage. Bias-reduced estimation based on Jeffreys prior is provided in addition to maximum likelihood estimation. The fitting function also properly handles predictors prone to measurement errors.

The zip file ``msbreg_0.1.0.tar.gz`` is the source R package (``msbreg``) for MSB models.
This package requires R (>= 4.0.2) and depend on a few other R packages. To install all dependencies:

``` r
# Install devtools if not yet
if(!require("devtools")) install.packages("devtools")

# Install dependencies for wavefinder
install.packages(c('Rdpack', 'methods', 'boot', 'parallel', 'RcppParallel', 'grDevices', 'graphics', 'lattice'), dependencies = TRUE)
```

To install ``msbreg`` on your computer, first download the files ``msbreg_0.1.0.tar.gz`` from the msbreg repository.

Then, run the following line with ``~`` replaced by the path to the directory where ``msbreg_0.1.0.tar.gz``was saved on your computer.

``` r
# Install the library mrbglm
install.packages("~/msbreg_0.1.0.tar.gz", repos = NULL, type = "source")
```

##** Infertility data
``` r
data ("infert", package = "datasets")

## Logistic regression fit to the infert data
GLMres <- glm (case ~ spontaneous,
               data = infert,
               family = binomial())

summary (GLMres)
```

#* MSB model fit with a maximum success probability in [0, 1)
``` r
require(msbreg)
MSBres0 <- msbreg (case ~ spontaneous,
                   lambda.formula = ~ 1,
                   start = c(GLMres$coefficients, `(Intercept).lambda` = 0),
                   criterion = "ML",
                   data = infert)

summary (MSBres0)
```

#* The ML estimate of the asymptote is lambda = 0.8469 < 1.
#* But adding the additional asymptote parameter does not seem justified
#* (by the larger AIC as compared with the simple GLM).

#* This is in line with the profile likelihood plot:
profilelambda(MSBres0, plotit = TRUE, lambda = seq(0.5, 1, .01))

#* We can use a penalized likelihood estimation to make the likelihood surface less flat (default procedure)
``` r
MSBres <- msbreg (case ~ spontaneous,
                  lambda.formula = ~ 1,
                  start = c(GLMres$coefficients, `(Intercept).lambda` = 0),
                  data = infert)

summary (MSBres)
```

#* The penalized ML estimate of lambda is 1, on parameter space boundary
#* (as can be seen on the profile likelihood plot)
``` r
profilelambda(MSBres, plotit = TRUE, lambda = seq(0.6, 1, .01))
```

#* Use the score statistic to test H0: lambda = 1 versus H1: 0 < lambda < 1.
``` r
Sctest <- score.test(MSBres)
Sctest # lambda is not significantly lower than one
```

#* Plot the two model fits with the observed data
``` r
courbe(MSBres, col = 'red', xlim = c(-4, 5))
courbe(glm.to.msbm(GLMres), xlim = c(-4, 5), col = 'blue', add = TRUE)
legend(x = -4, y = 0.9, legend = c("GLM", "MSB"),
       col = c("red", "blue"), lty = c(1, 1))
```

#* The two fits are different because of two distinct estimation methods:
#* Maximum Likelihood vs Maximum A Posteriori with Jeffreys invariant prior

#** Simulated example 1: two-stage logistic model
``` r
data(test1data)
attr(test1data$y, "formula")

# Standard logistic fit
GLMres <- glm (y/Total ~ x1 + offset(off1) + x2,
               weights = Total,
               data = test1data,
               family = binomial())

summary (GLMres)
```

#* Two-stage logistic model fit with a maximum success probability set to 1
#* No multiplicative intercept
``` r
MSBres0 <- msbreg (y/Total ~ x1 + offset(off1) | x2,
                  lambda.formula = ~ 0,
                  weights = Total,
                  data = test1data)

summary (MSBres0)
```

#* Two-stage logistic model fit with unknown maximum success probability
#* A multiplicative intercept included
``` r
MSBres <- msbreg (y/Total ~ x1 + offset(off1) | x2,
                   lambda.formula = ~ 1,
                   weights = Total,
                   data = test1data)

summary (MSBres)
```

#* AIC, HQC and BIC criteria: MSB models are fitter
``` r
cbind(AIC(GLMres, MSBres0, MSBres),
      HQC = HQC(GLMres, MSBres0, MSBres)[,2],
      BIC = BIC(GLMres, MSBres0, MSBres)[,2])
```

#* Pseudo R-squared
#* Squared-correlation between response and predictions
``` r
rsquared (GLMres, MSBres0, MSBres, method = "COR")

# Deviance reduction ratio (KL)
# (not appropriate for comparing GLMres vs MSBres0 or MSBres
#  because the null deviances differ)
rsquared (MSBres0, MSBres)
```

#* True coefficients versus estimates from the best fit
``` r
cbind(Truth = attr(test1data$y, "theta"),
      Estimate = MSBres$coefficients)
```

#* Fitting a model with known alpha and lambda components:
#* lambda = 0.85 and alpha = .05/.85 (so that 0.05 is the minimum probability)
``` r
MSBresk <- msbreg (y/Total ~ x1 + offset(off1) | x2,
                   alpha.formula = ~ offset(qlogis(0.05/0.85)) - 1,
                   lambda.formula = ~ offset(qlogis(0.85)) - 1,
                   weights = Total,
                   data = test1data)
# In 'alpha.formula' and 'lambda.formula':
# 'qlogis' transforms alpha or lambda to link scale (use 'qnorm' for probit link)
# 'offset' declares the transformed values as known, to be kept fixed.
# '-1' indicates the absence of an intercept (none to be estimated).

summary (MSBresk)
```

#* Note that the null deviance is different for this model as compared to other models.

``` r
AIC(GLMres, MSBres0, MSBres, MSBresk)
```

#** Simulated example 2: fitting MSB model
#*       with error-prone continuous predictors,
#*       and a categorical predictor.

#* Simulate predictors
``` r
set.seed(167)
mdf <- data.frame(x1 = rexp(1000, 0.5),
                  # First continuous predictor
                  x2 = rnorm(1000),
                  # Second continuous predictor
                  x3 = rpois(1000, 5),
                  # Numeric, but discrete predictor
                  f1 = rep(c('a', 'b', 'c'),
                           length.out = 1000))
                  # Factor with three levels
```

#* Generate a response vector using the predictors
#*   considering two stages
``` r
mframe <- sim.msb (formula = ~ x1 + f1 | x2 + x3,
                   data = mdf,
                   link = 'probit',
                   theta = c(3, -1, 0, 0.5, 3, -1, -0.5, 1.7),
                   seed = NULL)

mframe$frame # model frame
summary(mframe$y) # binary response
```

#* Add the simulated response to the dataset
``` r
mdf$y <- mframe$y
```

#* Simulate error-prone predictors
``` r
mdf$SDx1 <- abs(rnorm(1000, sd = 0.5))
mdf$SDx2 <- abs(rnorm(1000, sd = 0.5))
mdf$x1err <- rnorm(n = 1000, mean = mdf$x1, sd = mdf$SDx1)
mdf$x2err <- rnorm(n = 1000, mean = mdf$x2, sd = mdf$SDx2)
```

#* Fit the model, considering the error-prone predictors
``` r
MSBfit <- msbreg(formula = y ~ x1err + f1| x2err + x3,
                 input.me = list(x1err = ~ SDx1,
                                 x2err = ~ SDx2),
                 data = mdf,
                 link = 'probit',
                 criterion = 'ML')

summary (MSBfit)
```
#* Compare the coefficient estimates to the true
#* coefficients: (3, -1, 0, 0.5, 3, -1, -0.5, 1.7)

