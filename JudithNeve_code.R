# necessary libraries

library(bain)
library(scales)

################### Preliminary analysis ###################

# Model 0: preliminary model to select predictors

red <- read.csv("winequality-red.csv", sep = ";") # read in data

# center all predictors - so we will be able to compare betas
for (i in 1:ncol(red)){
  red[,i] <- (red[,i] - mean(red[,i]))/sd(red[,i])
}

mod.red <- lm(quality ~ .-1, data = red)
summary(mod.red)

# include the parameters with estimates that have absolute value larger than .1

included.pars <- names(coef(mod.red))[abs(coef(mod.red)) > 0.1]

################### Bayesian models ###################

# read in data
white <- read.csv("winequality-white.csv", sep = ";")

# center all predictors - so we will be able to compare betas
for (i in 1:ncol(white)){
  white[,i] <- (white[,i] - mean(white[,i]))/sd(white[,i])
}

# splitting the data into predictors and outcome
X <- as.matrix(white[,included.pars])
Y <- white[,"quality"]


# MCMC sampling function definition

mcmc.sampling <- function(Y, X, n.chains, burnin, iter, # outcome, predictor, number of chains, burnin length, and number of iterations
                          prior.mean.b, prior.sd.b, prior.shape.sigma, prior.scale.sigma, # prior distribution parameters
                          included.pars) { # names of parameters to include
  pred <- length(included.pars)
  prior.var.b <- prior.sd.b^2
  
  mcmc.sampler <- list() # prepare a list that will be filled with a matrix of sampled values for each chain
  
  for (ch in 1:n.chains) {
    # we create a matrix to fill with the sampled values: one row per iteration, one column per predictor
    chain <- matrix(NA, nrow = iter, ncol = pred+2)
    colnames(chain) <- c(included.pars, "Sigma", "Accepted.MH")
    # we define random starting values
    B.start <- runif(n = pred, min = -1, max = 1) # beta are bound between -1 and 1
    sig.start <- runif(n = 1, min = 0, max = 1) # sigma is bound between 0 and 1
    chain[1,] <- c(B.start, sig.start, NA) # the first row is the initial values
    
    # now we sample
    b.current <- B.start # in order to not write a separate step for each beta, we use a temporary vector with the most recent value for each beta within a loop
    for (i in 2:iter) {
      # the conditional posterior of each B
      for (j in 1:(pred-1)) { # all but the last beta are done using gibbs
        posterior.mean.b <- (sum(X[,j]*(Y-apply(b.current[-j]*t(X[,-j]), 2, sum))) / # the sum of the x_ij*(y - yhat), where yhat is estimated using the most recent beta values and without the relevant predictor
                               chain[i-1, pred+1]^2 + # divide the sum by the current estimate of sigma^2
                               prior.mean.b[j]/prior.var.b[j]) / # add the priod mean/prior variance
          (sum(X[,j]^2) / chain[i-1, pred+1]^2 + 1/prior.var.b[j]) # divide the whole by the sum of the squared predictor / sigma^2 + the inverse of the prior variance
        
        posterior.var.b <- 1/((sum(X[,j]^2)/chain[i-1, pred+1]^2 + 1/prior.var.b[j])) # posterior variance
        posterior.sd.b <- sqrt(posterior.var.b) # and we take the sqrt for the rnorm function
        
        chain[i,j] <- rnorm(1, posterior.mean.b, posterior.sd.b) # we draw from the posterior distribution
        b.current[j] <- chain[i,j] # then update the vector that is used within the loop
      }
      
      # the conditional posterior of the last beta - MH
      u <- runif(1, 0, 1) # get u
      
      beta.prop <- rnorm(1, prior.mean.b[pred], prior.sd.b[pred]) # proposal density is the same as the prior: best guess
      
      # calculate the density function - we use the log because the values are very small
      density.data <- 0
      density.data.t <- 0
      for (x in 1:length(Y)) {
        mean.prop <- sum(c(chain[i,1:(pred-1)], beta.prop)*X[x,]) # y.hat
        mean.t <- sum(c(chain[i,1:(pred-1)], b.current[pred])*X[x,])
        
        density.data <- density.data + dnorm(Y[x],
                                          mean = mean.prop,
                                          sd = chain[i-1,pred+1],
                                          log = TRUE)
        density.data.t <- density.data.t + dnorm(Y[x],
                                          mean = mean.t,
                                          sd = chain[i-1,pred+1],
                                          log = TRUE)
      }
      
      # true density: prior times density, or log prior + log density
      p.beta <- dnorm(beta.prop, prior.mean.b[pred], prior.sd.b[pred], log = TRUE) + density.data
      p.beta.t <- dnorm(b.current[pred], prior.mean.b[pred], prior.sd.b[pred], log = TRUE) + density.data.t
      # log pdf of the proposal density
      q.beta <- dnorm(beta.prop, prior.mean.b[pred], prior.sd.b[pred], log = TRUE)
      q.beta.t <- dnorm(b.current[pred], prior.mean.b[pred], prior.sd.b[pred], log = TRUE)
      
      AR <- p.beta - p.beta.t + q.beta.t - q.beta # acceptance ratio
      
      chain[i,pred] <- ifelse(log(u) <= AR, beta.prop, chain[i-1,pred]) # save either sigma or repeat the previous one
      b.current[pred] <- chain[i, pred]
      chain[i,pred+2] <- ifelse(log(u) <= AR, 1, 0) # 1 is accepted 0 is rejected
      
      # the conditional posterior of sigma^2 - back to gibbs
      posterior.shape.sigma <- length(Y)/2 + prior.shape.sigma
      posterior.scale.sigma <- sum((Y - apply(b.current*t(X), 2, sum))^2)/2 + prior.scale.sigma
      
      chain[i,pred+1] <- sqrt(1/rgamma(1, shape = posterior.shape.sigma, rate = posterior.scale.sigma))
    }
    
    mcmc.sampler[[ch]] <- chain[-(1:burnin),] # discard the burnin
  }
  return(mcmc.sampler)
}

# set seed, numbers of iterations, burnin, etc - these will be the same for all models
set.seed(0070661)
n.chains <- 2
burnin <- 500
iter <- 10000

################### Model 1: taste only ###################

included.pars.taste <- included.pars[c(1, 5)]
X.taste <- X[,included.pars.taste]

# set prior parameters - informative prior from the red wine dataset
prior.mean.b.taste <- coef(mod.red)[included.pars.taste] # set the prior mean of betas
prior.sd.b.taste <- 2*summary(mod.red)$coefficients[included.pars.taste,2] # set the prior sd of betas
prior.shape.sigma <- nrow(red)/2
prior.scale.sigma <- var(mod.red$residuals)*prior.shape.sigma

# run the sampler
taste.model <- mcmc.sampling(Y, X.taste, n.chains, burnin, iter,
                             prior.mean.b.taste, prior.sd.b.taste, prior.shape.sigma = prior.shape.sigma, prior.scale.sigma = prior.scale.sigma, included.pars.taste)

################### Model 2: taste excluded ###################

included.pars.taste.exc <- included.pars[-c(1, 5)]
X.taste.exc <- X[,included.pars.taste.exc]

# set prior parameters - informative prior from the red wine dataset
prior.mean.b.taste.exc <- coef(mod.red)[included.pars.taste.exc] # set the prior mean of betas
prior.sd.b.taste.exc <- 2*summary(mod.red)$coefficients[included.pars.taste.exc,2] # set the prior sd of betas
prior.shape.sigma <- nrow(red)/2
prior.scale.sigma <- var(mod.red$residuals)*prior.shape.sigma

# run the sampler
taste.exc.model <- mcmc.sampling(Y, X.taste.exc, n.chains, burnin, iter,
                                 prior.mean.b.taste.exc, prior.sd.b.taste.exc, prior.shape.sigma = prior.shape.sigma, prior.scale.sigma = prior.scale.sigma, included.pars.taste.exc)

################### Model 3: selected predictors ###################

# set prior parameters - informative prior from the red wine dataset
prior.mean.b <- coef(mod.red)[included.pars] # set the prior mean of betas
prior.sd.b <- 2*summary(mod.red)$coefficients[included.pars,2] # set the prior sd of betas
prior.shape.sigma <- nrow(red)/2
prior.scale.sigma <- var(mod.red$residuals)*prior.shape.sigma

# run the sampler
full.model <- mcmc.sampling(Y, X, n.chains, burnin, iter, prior.mean.b, prior.sd.b, prior.shape.sigma = prior.shape.sigma, prior.scale.sigma = prior.scale.sigma, included.pars)

################### Assess convergence of the models ###################

################### Model 1 ###################

### History plots

pred.taste <- ncol(X.taste)

par(mfrow = c(1, 3))
for (i in 1:(pred.taste+1)) {
  # plot one chain
  plot(taste.model[[1]][,i], type = "l", col = alpha("blue", 0.5), ylab = colnames(taste.model[[1]])[i], xlab = "")
  # add lines for the other chain
  lines(1:(iter-burnin), taste.model[[2]][,i], col = alpha("red", 0.5))
}

### Autocorrelations

autocorr.taste <- list() # create an empty list that will be filled with matrices of autocorrelations

for (j in 1:n.chains) {
  autocorr.taste[[j]] <- matrix(1, nrow = 30, ncol = pred.taste+1)
  for (i in 1:(pred.taste+1)) {
    for (k in 1:29) { # calculate the correlation up until lag 29 (30 autocorrelations total)
      autocorr.taste[[j]][k+1, i] <- sum((taste.model[[j]][-((iter-burnin+1-k):(iter-burnin+1)),i] - mean(taste.model[[j]][,i]))*((taste.model[[j]][-(1:k),i]) - mean(taste.model[[j]][,i]))) / sum((taste.model[[j]][,i] - mean(taste.model[[j]][,i]))^2)
    }
  }
}

# plot the autocorrelations
par(mfrow = c(1, 3))
for (i in 1:(pred.taste+1)) {
  barplot(autocorr.taste[[1]][,i], xlab = colnames(taste.model[[1]])[i]) # change 1 to 2 to see for the second chain
}

### Gelman & Rubin statistic

GR.taste <- rep(NA, pred.taste+1) # create empty vector for all parameters of the model
for (i in 1:(pred.taste+1)) {
  N <- iter-burnin+1
  var.chain <- rep(NA, n.chains)
  mean.chain <- rep(NA, n.chains)
  for (m in 1:n.chains) {
    mean.chain[m] <- mean(taste.model[[m]][,i])
    var.chain[m] <- var(taste.model[[m]][,i])
  }
  W <- sum(var.chain)/n.chains
  B <- N * sum((mean.chain - mean(mean.chain))^2)/(n.chains-1)
  V.hat <- (1-1/N)*W + B/N
  GR.taste[i] <- sqrt(V.hat/W)
}

GR.taste

### MC error

MC.taste <- matrix(NA, pred.taste+1, n.chains)
for (j in 1:n.chains) {
  for (i in 1:(pred.taste+1)) {
    MC.taste[i,j] <- sd(taste.model[[j]][,i])/sqrt(iter-burnin)
  }
}

MC.taste

################### Model 2 ###################

### History plots

pred.taste.exc <- ncol(X.taste.exc)

par(mfrow = c(2, 2))
for (i in 1:(pred.taste.exc+1)) {
  plot(taste.exc.model[[1]][,i], type = "l", col = alpha("blue", 0.5), ylab = colnames(taste.exc.model[[1]])[i], xlab = "")
  lines(1:(iter-burnin), taste.exc.model[[2]][,i], col = alpha("red", 0.5))
}

### Autocorrelations

autocorr.taste.exc <- list()

for (j in 1:n.chains) {
  autocorr.taste.exc[[j]] <- matrix(1, nrow = 30, ncol = pred.taste.exc+1)
  for (i in 1:(pred.taste.exc+1)) {
    for (k in 1:29) {
      autocorr.taste.exc[[j]][k+1, i] <- sum((taste.exc.model[[j]][-((iter-burnin+1-k):(iter-burnin+1)),i] - mean(taste.exc.model[[j]][,i]))*((taste.exc.model[[j]][-(1:k),i]) - mean(taste.exc.model[[j]][,i]))) / sum((taste.exc.model[[j]][,i] - mean(taste.exc.model[[j]][,i]))^2)
    }
  }
}

# plot the autocorrelation
par(mfrow = c(2, 2))
for (i in 1:(pred.taste.exc+1)) {
  barplot(autocorr.taste.exc[[1]][,i], xlab = colnames(taste.exc.model[[1]])[i]) # change 1 to 2 to see for the second chain
}

### Gelman & Rubin statistic

GR.taste.exc <- rep(NA, pred.taste.exc+1)
for (i in 1:(pred.taste.exc+1)) {
  N <- iter-burnin+1
  var.chain <- rep(NA, n.chains)
  mean.chain <- rep(NA, n.chains)
  for (m in 1:n.chains) {
    mean.chain[m] <- mean(taste.exc.model[[m]][,i])
    var.chain[m] <- var(taste.exc.model[[m]][,i])
  }
  W <- sum(var.chain)/n.chains
  B <- N * sum((mean.chain - mean(mean.chain))^2)/(n.chains-1)
  V.hat <- (1-1/N)*W + B/N
  GR.taste.exc[i] <- sqrt(V.hat/W)
}

GR.taste.exc

### MC error

MC.taste.exc <- matrix(NA, pred.taste.exc+1, n.chains)
for (j in 1:n.chains) {
  for (i in 1:(pred.taste.exc+1)) {
    MC.taste.exc[i,j] <- sd(taste.exc.model[[j]][,i])/sqrt(iter-burnin)
  }
}

MC.taste.exc

################### Model 3 ###################

### History plots

pred <- ncol(X)

par(mfrow = c(2, 3)) # so we can see all the plots at once
for (i in 1:(pred+1)) {
  # plot one chain first
  plot(full.model[[1]][,i], type = "l", col = alpha("blue", 0.5),
       ylab = colnames(full.model[[1]])[i], xlab = "")
  # then add the other chain
  lines(1:(iter-burnin), full.model[[2]][,i], col = alpha("red", 0.5))
}

### Autocorrelations

autocorr <- list() # make an empty list that will be filled with a matrix of autocorrelations per chain

for (j in 1:n.chains) {
  # create a matrix: each row is the autocorrelation at a given lag, each column is a predictor
  autocorr[[j]] <- matrix(1, nrow = 30, ncol = pred+1)
  for (i in 1:(pred+1)) {
    for (k in 1:29) {
      # autocorrelation formula
      autocorr[[j]][k+1, i] <- sum((full.model[[j]][-((iter-burnin+1-k):(iter-burnin+1)),i] - mean(full.model[[j]][,i]))*((full.model[[j]][-(1:k),i]) - mean(full.model[[j]][,i]))) / sum((full.model[[j]][,i] - mean(full.model[[j]][,i]))^2)
    }
  }
}

par(mfrow = c(2, 3)) # so that we can see all plots at once
for (i in 1:(pred+1)) {
  barplot(autocorr[[1]][,i], xlab = colnames(full.model[[1]])[i]) # change 1 to 2 to see for the second chain
}

### Gelman & Rubin statistic

GR <- rep(NA, pred+1) # empty vector to fill with the values
for (i in 1:(pred+1)) { # for each predictor:
  N <- iter-burnin+1 # kept iterations per sample
  var.chain <- rep(NA, n.chains) # variance of the predictor in each chain
  mean.chain <- rep(NA, n.chains) # mean of the predictor in each chain
  for (m in 1:n.chains) {
    mean.chain[m] <- mean(full.model[[m]][,i])
    var.chain[m] <- var(full.model[[m]][,i])
  }
  W <- sum(var.chain)/n.chains
  B <- N * sum((mean.chain - mean(mean.chain))^2)/(n.chains-1)
  V.hat <- (1-1/N)*W + B/N
  GR[i] <- sqrt(V.hat/W)
}

GR
max(abs(1-GR)) # maxiumum deviation from 1

### MC error

MC <- matrix(NA, pred+1, n.chains)
for (j in 1:n.chains) {
  for (i in 1:(pred+1)) {
    MC[i,j] <- sd(full.model[[j]][,i])/sqrt(iter-burnin)
  }
}

MC

total.sample <- rbind(full.model[[1]], full.model[[2]]) # combine the two chains into a matrix

comp.MC <- apply(total.sample, 2, sd)*0.05 # 5% of the SD

################### Checking there is a linear relationship between each predictor and the outcome ###################

# Step 1 - the null hypothesis: Y_i = sum beta_j*x_j + e_i, e_i ~ N(0, sigma^2)
# Step 2 - sample from the posterior distribution: full.model
# Step 3 - generate a dataset for each sample

### making a function to generate a dataset
generate.dataset <- function(mcmc.sampler, X) { # sampled values from the posterior and predictors
  samples <- rbind(mcmc.sampler[[1]], mcmc.sampler[[2]]) # get all the sampled values into one matrix
  p <- ncol(samples)-2 # number of predictors
  iter <- nrow(samples) # number of samples
  n <- nrow(X) # sample size of data
  
  posterior.predictive.data <- matrix(NA, nrow = iter, ncol = n) # a matrix where each row will be filled with a generated y vector
  
  for (i in 1:iter) {
    for (j in 1:n) {
      y.hat <- sum(X[j,]*samples[i,1:p]) # expected value
      y.sim <- rnorm(1, y.hat, samples[i,p+1]) # draw a value
      posterior.predictive.data[i,j] <- y.sim # save the value
    }
  }
  
  return(posterior.predictive.data)
}

### generate the dataset
full.mod.postpred <- generate.dataset(full.model, X)

# Step 4 - test statistic: the correlation of each predictor with Y (simulated), absolute value, then the number of these > 0.1

### making a function for the test statistic
test.statistic <- function(posterior.predictive, X, Y) {
  n.postpred <- nrow(posterior.predictive) # how many samples do we have?
  n <- nrow(X) # how many observations in each sample?
  
  test.statistics <- rep(NA, n.postpred)
  
  for (i in 1:n.postpred) {
    test.statistics[i] <- sum(abs(cor(X, posterior.predictive[i,])) > 0.2) # if all the predictors are relevant, this is p. then p > in practice, i.e. how often does ideal data fit better than true data
  }
  observed <- sum(abs(cor(X, Y)) > 0.2)
  
  output <- list(observed = observed, test.statistics = test.statistics)
  
  return(output)
}

### getting the test statistic
full.mod.test <- test.statistic(full.mod.postpred, X, Y)

# Step 5 - ppp-value: the proportion of the test statistic that is larger than the observed

### making a function
ppp.value <- function(model.test) {
  sum(model.test$test.statistics > model.test$observed)/length(model.test$test.statistics) 
}

### getting the ppp-value
ppp.full.mod <- ppp.value(full.mod.test)


################### Parameter estimates, credible intervals ###################

par.est <- apply(total.sample, 2, mean) # posterior mean
credible.interval <- apply(total.sample, 2, quantile, probs = c(0.025, 0.975)) # credible interval

################### Model comparison ###################

################### DIC ###################

# make a DIC function
DIC <- function(mcmc.sampler, Y, X){ # the sampled values (as a matrix), outcome and predictors
  p <- ncol(mcmc.sampler)-2
  iter <- nrow(mcmc.sampler)
  beta.bar <- apply(mcmc.sampler, 2, mean)[1:p] # regression coefficient estimates
  sigma.bar <- apply(mcmc.sampler, 2, mean)[p+1] # residual sd estimate
  
  # calculate the density for each observation, using the posterior means and all sampled values
  density.mean <- rep(NA, length(Y)) # using the posterior means
  density.iteration <- matrix(NA, length(Y), iter) # 1 row per observation, 1 column per iteration of the sampler
  for (i in 1:length(Y)) {
    density.mean[i] <- dnorm(Y[i], mean = sum(X[i,]*beta.bar), sd = sigma.bar, log = TRUE)
    for (j in 1:iter) {
      density.iteration[i,j] <- dnorm(Y[i], mean = sum(X[i,]*mcmc.sampler[j,1:p]), sd = mcmc.sampler[j,p+1], log = TRUE)
    }
  }
  logdensity.mean <- sum(density.mean) # get the overall density (using the posterior means)
  logdensity.iteration <- apply(density.iteration, 1, sum) # get the overall density of each run of the sampler
  
  d.bar <- -2*logdensity.mean
  d.theta <- -2*(sum(logdensity.iteration))
  
  pd <- d.theta/iter - d.bar # Spiegelhalter et al., 2002
  
  DIC <- d.bar + 2*pd
  
  return(DIC)
}

# get the DIC of each model
dic.full.mod <- DIC(total.sample, Y, X)
dic.taste.mod <- DIC(rbind(taste.model[[1]], taste.model[[2]]), Y, X.taste)
dic.taste.exc.mod <- DIC(rbind(taste.exc.model[[1]], taste.exc.model[[2]]), Y, X.taste.exc)

################### Bayes factor ###################

BF <- bain(x = abs(par.est[1:pred]), # using the MCMC sampler
           # hypothesis 1; hypothesis 2
           "(volatile.acidity, alcohol) > (chlorides, total.sulfur.dioxide, sulphates); (volatile.acidity, alcohol) < (chlorides, total.sulfur.dioxide, sulphates)",
           n = 4898, Sigma = cov(total.sample)[1:pred, 1:pred])
BF

