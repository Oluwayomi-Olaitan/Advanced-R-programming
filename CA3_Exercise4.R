## Optimising the log-likelihood function for a multivariate linear regression model.

df <- with(mtcars, data.frame(y=mpg, x1=disp, x2=hp, x3=wt))

# Writing a function nll_lm()

nll_lm <- function(data, par) {
  # with 3 predictors, we have 4 beta estimates
  beta <- par[1:(length(par) - 1)]
  sigma2 <- exp(par[length(par)]) #reparameterised the variance
  
  # negative likelihood in terms of the residuals
  
  # constructing a matrix for the predictor variables and a vector for y
  X <- as.matrix(cbind(1, data$x1, data$x2, data$x3))
  y <- data$y
  
  # estimating y_hat
  y_hat <- X %*% beta
  
  # estimating the residuals
  residuals <- y - y_hat
  
  #the negative log-likelihood, given e_i = N(0, sigma2)
  llik <- dnorm(residuals, mean = 0, sd = sqrt(sigma2), log = TRUE)
  return(-sum(llik))
  
  # llik <- (length(data$y) / 2) * log(2 * pi * sigma2) + (1 / (2 * sigma2)) * sum(residuals^2)
  # return(llik)
  
}

# initial guesses
# as instructed the mean(df$y) as the initial guess of the intercept parameter,
#hence the slopes and the initial guess for sigma is 0.
inits <- c(mean(df$y), rep(0,4))

fit <- optim(par = inits, 
             fn = nll_lm, 
             data = df, 
             hessian=TRUE,
             method="L-BFGS-B", 
             lower=c(-Inf, -Inf, -Inf, -Inf, 0), 
             upper=c(Inf, Inf, Inf, Inf, Inf)
)


fit$par


### Comparing the $\hat{\beta}$ 

# matrix multiplication
X <- as.matrix(cbind(1, df$x1, df$x2, df$x3))

y <- df$y

XTX <- solve(crossprod(X))
XTY <- crossprod(X, y)

Beta <- solve(crossprod(X), crossprod(X, y))
Beta

fit$par[1:4]

all.equal(as.vector(Beta), fit$par[1:4])

### Comparing the $\hat{\sigma}$
  
# the residual squared
n <- length(df$y)
p <- length(Beta)
residual_df <- n-p

y_hat <- X %*% Beta
residuals <-  y - y_hat

sigma_squared <- sum(residuals^2) / residual_df
sigma_hat <- sqrt(sigma_squared)
sigma_hat

fit$par[5]


identical(sigma_squared, fit$par[5])

### Standard errors of the regression coefficients
  
sqrt(diag(solve(fit$hessian)))


################################################################################
# Fit the linear model using lm()
lm_fit <- lm(y ~ x1+x2+x3, data = df)
  
# Extracting the estimated coefficients
beta_hat <- coef(lm_fit)
beta_hat

# Extracting the sigma_hat
sigma_hat <- summary(lm_fit)$sigma
sigma_hat

