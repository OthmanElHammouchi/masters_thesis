library(patternBreak)
suppressPackageStartupMessages(library(ChainLadder))

triangle <- UKMotor
triangle <- cum2incr(triangle)
triangle[is.na(triangle)] <- 0
ndev <- nrow(triangle)
nboot <- 1e3

n_dev <- 7
N <- 28
p <- 13

X <- matrix(rep(0, N * p), nrow = N)
y <- numeric(N)

X[, 1] <- 1

k <- 1
for (i in 1:n_dev) {
    for (j in 1:(n_dev + 1 - i)) {
        cat(sprintf("i: %i, j: %i\n", i, j))
        if (i != 1) X[k, i] <- 1
        if (j != 1) X[k, n_dev + j - 1] <- 1
        y[k] <- triangle[i, j]
        k <- k + 1
    }
}

y[1] <- triangle[1, 1]

data <- data.frame(as.data.frame(X), response = y)

# trace(glm.fit, quote(cat <- function(...) {
#   base::cat(...)
#   if (...length() >= 3 && identical(..3, " Iterations - ")) print(coefold)
# }))
# glm.out = glm(response ~ .,
#                      family=poisson(), data=data,
#                      control = glm.control(trace = TRUE))
# untrace(glm.fit)

betas<-c(log(mean(y)), rep(0, p-1))
eta<-X%*%betas

ylogy<-function(y)
{
return(ifelse(y==0,rep(0,length(y)),y*log(y)))
}

deviance<-2*sum(ylogy(y)-y*eta-(y-exp(eta)))

devianceOld<-1e30

tol<-1e-6
iteration=0
# while(((devianceOld-deviance)/devianceOld)>tol)
# {
# for (i in 1:10) {
diff <- 1e5
while(diff > tol) {
#start IRLS UPDATE STEP
iteration=iteration+1
#calculate pseudo data based on current betas
z=eta+exp(-eta)*(y-exp(eta))
#calculate new weights: diagonal elements
w<-c(exp(eta))

#update betas
# lmUpdate<-lm(z~-1+X,weight=w)
eta<-X%*%betas
# eta<-lmUpdate$fitted
# betas<-lmUpdate$coef
old.betas <- betas
betas <- solve(t(X) %*% diag(w) %*% X, t(X) %*% diag(w) %*% z)

#criterion for convergence 
devianceOld<-deviance
deviance<-2*sum(ylogy(y)-y*eta-(y-exp(eta)))
cat("iteration",iteration,"Deviance Old",devianceOld,"Deviance", deviance,"\n")

diff <- norm(betas - old.betas, type="2")
}