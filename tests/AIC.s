require(rms)
set.seed(1)
x <- runif(20)
y <- runif(20)
f <- lm(y ~ pol(x,2))
logLik(f)
AIC(f)
g <- ols(y ~ pol(x,2))
logLik(g)
AIC(g)

y <- sample(0:1, 20, replace=TRUE)
f <- glm(y ~ pol(x,2), family=binomial)
logLik(f)
AIC(f)
g <- lrm(y ~ pol(x,2))
logLik(g)
AIC(g)
g$stats
AIC(g, type='chisq')