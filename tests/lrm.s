require(rms)
n <- 50000
x1 <- runif(n)
x2 <- runif(n)
x3 <- runif(n)
L <- x1 + x2 + x3 - 1.5
y <- ifelse(runif(n) <= plogis(L), 1, 0)
system.time(f <- glm(y ~ x1 + x2 + x3, family=binomial))
print(summary(f), digits=7)
system.time(g <- lrm(y ~ x1 + x2 + x3))
print(g, digits=7)
coef(f) - coef(g)
sqrt(diag(vcov(f)))/sqrt(diag(vcov(g)))

require(MASS)
n <- 300
y <- factor(sample(0:4, n, TRUE))
x1 <- runif(n)
x2 <- runif(n)
x3 <- runif(n)
system.time(f <- polr(y ~ x1 + x2 + x3))
print(summary(f, digits=7))
system.time(g <- lrm(y ~ x1 + x2 + x3))
print(g, digits=7)
c(-f$zeta, f$coefficients) - coef(g)
print( (diag(vcov(f))[c(4:7, 1:3)])/diag(vcov(g)), digits=10)

w <- function(m) {
  x <- runif(200)
  if(m > 0) x[1:m] <- NA
  x
}
set.seed(1)
y <- sample(0:1, 200, TRUE)
x1 <- w(50)
x2 <- w(1)
x3 <- w(2)
x4 <- w(0)
x5 <- w(10)
x6 <- w(11)
x7 <- w(13)
x8 <- w(8)
x9 <- w(7)
x10 <- w(6)
x11 <- w(5)
x12 <- w(4)
x13 <- w(3)
x14 <- w(7)
x15 <- w(18)
x16 <- w(19)
x17 <- w(21)
x18 <- w(23)
x19 <- w(25)
x20 <- w(27)
f <- lrm(y ~ x1+x2+x3+x4+x5+x6+x7+x8+x9+x10+x11+x12+x13+x14+x15+x16+x17+x18+x19+x20)

source('~/R/rms/R/lrm.s');source('~/R/Hmisc/R/latexDotchart.s');source('~/R/rms/R/rmsMisc.s')
sink('/tmp/t.tex')
cat('\\documentclass{report}\\usepackage{color,epic,longtable}\\begin{document}',
    sep='\n')
print(f, latex=TRUE)
cat('\\end{document}\n')
sink()

