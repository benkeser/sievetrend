Installation of `sievetrend`
----------------------------

The `sievetrend` repository is available for download as an `R` package and may be downloaded directly from GitHub as follows.

``` r
# install sievetrend from GitHub
devtools::install_github("benkeser/sievetrend")
```

    ## Downloading GitHub repo benkeser/sievetrend@master
    ## from URL https://api.github.com/repos/benkeser/sievetrend/zipball/master

    ## Installing sievetrend

    ## '/Library/Frameworks/R.framework/Resources/bin/R' --no-site-file  \
    ##   --no-environ --no-save --no-restore --quiet CMD INSTALL  \
    ##   '/private/var/folders/9d/nqhs6v950v39sthl1pljm6s40000gp/T/RtmpoSwAio/devtools19bb729b17e/benkeser-sievetrend-ff4f5b2'  \
    ##   --library='/Library/Frameworks/R.framework/Versions/3.4/Resources/library'  \
    ##   --install-tests

    ## 

``` r
library(sievetrend)
# also will require package survtmle
if(!("survtmle" %in% row.names(installed.packages()))){
    install.packages("survtmle")    
}
library(survtmle)
```

    ## survtmle: Targeted Learning for Survival Analysis

    ## Version: 1.0.0

Below, we include the code that was used to perform the simulation study and to analyze the RTS,S data for the manuscript.

Simulation
----------

Here, we demonstrate how to obtain results for a simulated data set.

``` r
# sample size
n <- 1000

# set random seed
set.seed(1234)

# simulate data set
dat <- makeData(n = n)

# get the formula for the empirical estimate of the
# censoring distribution needed for input to survtmle
glm.ctime <- get.ctimeForm(trt = dat$trt, site = dat$adjustVars$site, 
                           ftime = dat$ftime, ftype = dat$ftype)

# call survtmle to estimate cumulative incidence
object <- survtmle(ftime = dat$ftime,
                  ftype = dat$ftype,
                  adjustVars = dat$adjustVars,
                  trt = dat$trt,
                  glm.trt = "1",
                  glm.ftime = "trt*factor(site)",
                  glm.ctime = glm.ctime,
                  method = "mean",
                  t0=6)

# call trend_test to estimate projection
trend <- trend_test(object)
trend
```

    ## Trend in efficacy across failure type levels: 
    ##     beta lower_95%CI upper_95%CI  pval
    ## 1 -0.254      -0.486      -0.022 0.032

A numerical approximation of the true value of *β*<sub>0, *n*</sub> may be obtained as follows.

``` r
# get covariance matrix estimate
nabla_g <- sievetrend:::grad_g(object$est)
Upsilon_n <- nabla_g %*% cov(Reduce(cbind, object$ic)) %*% t(nabla_g) 
# call getTruth function with this estimate several times 
# and average (to increase accuracy)
beta_0n <- rowMeans(replicate(20, getTruth(Upsilon = Upsilon_n, n = 1e6)))[2]

# compare estimate to truth
c(est = trend$beta, truth = beta_0n)
```

    ##        est      truth 
    ## -0.2540477 -0.2001120

The full code used to execute the simulation study is included in the simulation subdirectory. This includes the script `cent.R`, which is batched to a slurm via the script `sce.sh`. The function `makeIllustrationPlot` was used to produce Figure 1.

RTS,S Analysis
--------------

Due to existing privacy agreements, it is difficult to obtain access to the real RTS,S data. Nevertheless, a mock data set is distributed with the `survtmle` package that will serve as illustration for the outputation methodology. See `?rtss` for further description of the data set.

The `rtss` data set is a list of ten data sets each representing an outputed data set. They are formatted for analysis of a binary genetic mark, so we first replace the failure type column with a simulated genetic distance.

``` r
# replace the ftype column in each data set
rtss_mod <- lapply(rtss, function(data){
    # which were observed failures
    fail_idx <- which(data$ftype > 0)
    # number of observed failures
    fail_n <- length(fail_idx)
    # replace with a simulated version
    data$ftype[fail_idx] <- rbinom(fail_n, 4, plogis(data$vaccine)) + 1
    # collapse site variable into a single column
    data$site <- as.numeric(1*(data$site1==1) + 2*(data$site2==1) + 
                            3*(data$site3==1) + 4*(data$site4==1) + 
                            5*(data$site5==1))
    return(data)
})

# check out new ftype distribution for first
# outputed data set
table(rtss_mod[[1]]$ftype)
```

    ## 
    ##    0    1    2    3    4    5 
    ## 4806   48  272  583  746  435

Now we use `survtmle` to estimate the cumulative incidence for each data set.

``` r
rslt <- lapply(rtss_mod, function(data){
    glm.ctime <- get.ctimeForm(trt = data$trt, site = data$site, 
                           ftime = data$ftime, ftype = data$ftype)

    object <- survtmle(ftime = data$ftime,
                  ftype = data$ftype,
                  adjustVars = data[,"site",drop=FALSE],
                  trt = data$vaccine,
                  glm.trt = "1",
                  glm.ftime = "trt*factor(site)",
                  glm.ctime = glm.ctime,
                  method = "mean",
                  t0=6)
    return(object)
})
```

We can use the `getMO` function to obtain the averaged cumulative incidence results and apply the trend\_test on this scale.

``` r
# obtain averaged results
mo_rslt <- getMO(rslt)

# apply trend test
mo_trend <- trend_test(mo_rslt)
mo_trend
```

    ## Trend in efficacy across failure type levels: 
    ##     beta lower_95%CI upper_95%CI  pval
    ## 1 -0.021      -0.916       0.874 0.963
