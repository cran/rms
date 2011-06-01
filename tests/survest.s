# From John <arrayprofile@yahoo.com>
require(rms)
set.seed(1)
rstime <- rnorm(50,30,4)
rsentry <- rstime-runif(50,2,5)
#rstatus <- as.numeric(rnorm(50) > 0.1)
rstatus <- sample(0:1, 50, replace=TRUE)
#trt <- as.factor(as.numeric(rnorm(50) > 0.15))
trt <- sample(0:1, 50, replace=TRUE)
str <- as.factor(sample(c('a','b','c'), 50, replace=TRUE))

    ### This is the one works
fit <- cph(Surv(rstime,rstatus) ~ trt + strat(str), x=TRUE, y=TRUE)
survest(fit, expand.grid(trt=c(0,1), str=c('a','b','c')),
        times=c(20,25,30), conf.int=.95)

    ### This one doesn't
fit <- cph(Surv(rsentry,rstime,rstatus) ~ trt + strat(str),
           x=TRUE, y=TRUE)
survest(fit, expand.grid(trt=c(0,1), str=c('a','b','c')),
        times=c(1,2,3), conf.int=.95)
