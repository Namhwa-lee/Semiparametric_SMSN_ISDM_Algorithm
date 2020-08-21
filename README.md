# ISDM_Algorithm_Semiparametric_SMSN
# Intra Simplex Direction Method for calculating NPMLE of Semiparametric Scale Mixture of Skew Normal #

# Install Packages #
library(sn)
library(moments)
library(histogram)

# Make Simulation Data : 2 components Scale Mixture of Skew Normal #
n=1000
set.seed(1)
Z = sample(c(1:2), size=n, prob=c(.6,.4), replace=TRUE)
table(Z)/n

X = matrix(0, nrow=n)
dataset = data.frame(X,Z)

set.seed(1)
for (i in 1:n){
  if (dataset[i,2]==1) {dataset[i,1] = rsn(1, xi=0, omega = 1, alpha = 5)}
  else if (dataset[i,2]==2) {dataset[i,1] = rsn(1, xi=0, omega = 3, alpha = 5)}
}

data <- as.matrix(dataset[,1])

# Plot the data #
hist.data <- histogram(data, type="regular")
hist(data, breaks=hist.data$breaks, ylim=c(0,0.5), freq = FALSE, main="Histogram of Simulation Data", xlab="X")
curve(0.6*dsn(x,xi=0, omega = 1, alpha = 5) + 0.4*dsn(x,xi=0, omega = 3, alpha = 5), col="red", lwd=2, n=1000, add=TRUE)

# Make ISDM algorithm for SNMIX #
snmix.isdm <- function(data, init.lambda=NULL, init.sigma=NULL, init.xi=NULL, grid, constant) {
  
  # Basic Setting #
  
  nr <- nrow(data) ; iter <- 1 # 1st component with initial scale parameter
  
  prob = location = scale = skew <- vector(mode="list", length = 30) # if needed, increase the length of the list
  
  Dp <- 100 ; Dp_value <- c() # Directional Derivative
  
  obs.loglik = function(f) {
    sum(log(rowSums(f)))
  }
  
  # Assign any initial value to the scale parameter # 
  # In this algorithm, by using the moments of skew normal distribution, I solved the equation and assigned #
  
  find_init.lambda <- function(y, lambda) {
    
    gamma <- skewness(y)
    
    numerator <- sqrt(2) * (4-pi) * lambda^3
    denominator <- (pi + (pi-2)*lambda^2)*(3/2)
    
    return(gamma-(numerator/denominator))
    
  }
  
  find_init.sigma <- function(y, sigma, lambda) {
    
    # Be careful : Input sigma not sigma^2 #
    
    var(y) - (1 - (2/pi)*((lambda/sqrt(1+lambda^2))^2)) * sigma^2
    
  }
  
  find_init.xi <- function(y, xi, sigma, lambda) {
    
    mean(y) - xi - sqrt(2/pi) * (lambda/sqrt(1+lambda^2)) * sigma
    
  }
  
  if (is.null(init.lambda) == TRUE)
  { skew[[1]] <- uniroot(find_init.lambda, interval = c(-100000,100000), y=data, tol=1e-10)$root }
  
  if (is.null(init.sigma) == TRUE)
  { scale[[1]] <- uniroot(find_init.sigma, interval = c(0,100000), y=data, lambda=skew[[1]], tol=1e-10)$root }
  
  if (is.null(init.xi) == TRUE)
  { location[[1]] <- uniroot(find_init.xi, interval = c(-100000, 100000), y=data, sigma=scale[[1]], lambda=skew[[1]], tol=1e-10)$root }
  
  # Estimate initial parameters for ISDM algorithm by using EM-Algorithm #
  
  # Initial Setting #
  
  f_init <- matrix(0, nrow=nr)
  
  obs.lik_1st <- c()
  
  # 1st Loop #
  f_init[,1] <- matrix(dsn(x=data, xi=location[[1]], omega=scale[[1]], alpha=skew[[1]]), nrow=nr)
  
  obs.lik_1st <- obs.loglik(f=f_init)
  
  delta <- function(lambda) {lambda/sqrt(1+lambda^2)} # Be careful when making multivariate case : compute matrix form
  
  equation <- function(data, lambda, xi, sigma_2, s1ij, s2ij) {
    y <- data
    delta <- lambda/sqrt(1+lambda^2)
    first <- sigma_2 * delta * (1-delta^2) * nr
    second <- (1+delta^2) * sum((y-xi)*s1ij)
    third <- delta * sum(s2ij)
    fourth <- delta * sum(((y-xi)^2))
    value=first+second-third-fourth
    value
  }
  
  mu_tau <- matrix(0, nrow=nr, ncol=1)
  sigma_tau <- matrix(0, nrow=1)
  numerator = denominator <- matrix(0, nrow=nr, ncol=1)
  s1ij = s2ij <- matrix(0, nrow=nr, ncol=1)
  sigma_2 <- array(0, dim=c(1,1,1))
  
  error <- 100 
  
  # Loop Statement : Find mixing probability by EM-Algorithm #
  while (error > 1e-08 & length(obs.lik_1st) < 10000) {
    
    # 1st Iteration #
    
    ### E-step ###
    # Update mu_tau_j #
    mu_tau[,1] <- delta(lambda=skew[[1]]) * (data - location[[1]])
    
    # Update sigma_tau #
    sigma_tau <- scale[[1]] * sqrt(1-delta(lambda=skew[[1]])^2)
    
    # Update s1ij & s2ij #
    
    numerator[,1] <- dnorm(skew[[1]] * ((data-location[[1]])/scale[[1]] ))
    denominator[,1] <- pnorm(skew[[1]] * ((data-location[[1]])/scale[[1]] ))
    
    s1ij[,1] <- (mu_tau[,1] + (numerator[,1]/denominator[,1]) * sigma_tau)
    s2ij[,1] <- (mu_tau[,1]^2 + sigma_tau^2 + (numerator[,1]/denominator[,1]) * mu_tau[,1] * sigma_tau)
    
    # CM-step1 : Update xi #
    location[[1]] <- (sum(data)-delta(lambda=skew[[1]]) * sum(s1ij[,1]))/nr 
    
    # CM-step2 : Update sigma #
    sigma_2[,,1] <- ((sum(s2ij[,1]) - 2*delta(lambda=skew[[1]])*sum((data-location[[1]])*s1ij[,1]) + sum((data-location[[1]])^2)) / (2*(1-(delta(lambda=skew[[1]]))^2)*nr))
    
    scale[[1]] <- as.vector(sqrt(sigma_2))
    
    # CM-step4 : Update lambda #
    skew[[1]] <- uniroot(equation, interval = c(-100000,100000), data=data, sigma_2=sigma_2[,,1], xi=location[[1]], s1ij=s1ij[,1], s2ij=s2ij[,1], tol=1e-10)$root
    
    # Stack Observed Log-Likelihood
    
    f_init[,1] <- matrix(dsn(data, xi=location[[1]], omega=scale[[1]], alpha=skew[[1]]), nrow=nr)
    
    denom = rowSums(f_init)
    
    zhat = f_init/denom
    
    sum.z <- colSums(zhat)
    
    obs.lik_1st <- c(obs.lik_1st, obs.loglik(f=f_init)) # check that the EM algorithm works well using below footnote code #
    
    # all(diff(obs.lik_1st)>0) #
    
    # Record Error #
    error <- obs.lik_1st[length(obs.lik_1st)] - obs.lik_1st[length(obs.lik_1st)-1]
    
  }
  
  # Print the initial parameters : uni-component skew normal distribution #
  
  print(rbind(location=location[[1]], scale=scale[[1]], skewness=skew[[1]]))
  
  # Initial Distribution : 1st step #
  
  given.dist = NPMLE <- matrix(f_init, nrow=nr) # using skew[[1]] & location[[1]] & scale[[1]]
  
  obs.lik_final <- sum(log(NPMLE)) # observed log-likelihood of initial distribution
  
  # Directional Derivative Function Depends on scale parameter : sigma #
  
  Dp_sigma <- function(xi, sigma, lambda, f_p){
    f_scale <- matrix(dsn(x=data, xi=xi, omega=sigma, alpha=lambda), nrow=nr)
    return((sum(f_scale/f_p) - nr))
  }
  
  # Grid Search & Find Directional Derivative #
  
  sigma.grid <- seq(from=1e-05, to=grid, by=0.01) # if wanted, split the range more densely.
  
  # Function for finding all local maxima #
  
  local.max <- function(x, xi, lambda) {
    rate <- sign(diff(x))
    num <- 0
    for (i in 2:length(rate)) {
      if (rate[i-1]==1 & rate[i-1] != rate[i]) {num <- num+1}
    }
    change <- vector(mode="list", length=num)
    occur <- 0
    for (i in 2:length(rate)) {
      if (rate[i-1]==1 & rate[i-1] != rate[i]) {occur <- occur+1} & {change[[occur]] <- c(i-1, i)}
    }
    loc.max <- matrix(0, ncol=2, nrow=length(change))
    for (k in 1:length(change)) {
      loc.max[k,] <- unlist(optimize(Dp_sigma, lower=sigma.grid[change[[k]][1]], upper=sigma.grid[change[[k]][2]], xi=xi, lambda=lambda, f_p=given.dist, maximum=TRUE))
    }
    colnames(loc.max) <- c("maxima", "dir_deriv")
    loc.max[loc.max[,2]>0,]
  }
  
  # Function for finding location & skewness parameter which make the maximum likelihood #
  
  find.xi_lambda <- function(param, prob, support, data){
    obs.lik_param <- c()
    xi <- param[1] ; lambda <- param[2]
    f_obj <- matrix(0, nrow=nr, ncol=length(support))
    for (i in 1:length(support)){
      f_obj[,i] <- prob[i]*matrix(dsn(x=data, xi=xi, omega=support[i], alpha=lambda))
    }
    obs.lik_param <- obs.loglik(f=f_obj)
    obs.lik_param
  }
  
  # Loop statement : Stopping Rule : Among Local Maxima, Max(Directional Derivative) is less than 1 #
  
  prob[[1]] <- 1
  
  # Estimation of Semi-parametric SMSN by using Directional Derivative #
  
  while (Dp > 1 | iter == nr)  {
    
    iter <- iter + 1 
    
    value <- matrix(0, nrow=length(sigma.grid))
    
    for (i in 1:length(sigma.grid)) {
      value[i,] <- Dp_sigma(sigma=sigma.grid[i], xi=location[[iter-1]], lambda=skew[[iter-1]], f_p=NPMLE)
    }
    dir_deriv <- matrix(cbind(sigma.grid, value), nrow=length(sigma.grid))
    
    plot(dir_deriv, type="l", lwd=2, xlab="scale parameter", ylab='Directional Derivative', main="Directional Derivative w.r.t scale parameter")
    
    maxima <- matrix(local.max(x=value, xi=location[[iter-1]], lambda=skew[[iter-1]]),ncol=2); colnames(maxima) <- c("maxima", "dir_deriv")
    
    for(i in 1:nrow(maxima)){
       points(maxima[i,"maxima"], Dp_sigma(sigma=maxima[i,"maxima"], xi=location[[iter-1]], lambda=skew[[iter-1]], f_p=NPMLE), col="red", cex=1.5, lwd=2)
    }
    
    if (nrow(maxima)>10) {print("Not enough dimension to include all local maximas")} & {break}
    
    scale[[iter]] <- maxima[,1] + constant
    
    support <- unlist(scale)
    
    support <- rbind(support, dir_deriv=0)
    
    for (i in 1:ncol(support)){
      
      support["dir_deriv",i] <- Dp_sigma(sigma=support["support",i], xi=location[[iter-1]], lambda=skew[[iter-1]], f_p=NPMLE)
      
    }
    
    # Find Mixing Probability by Using EM-Algorithm #
    
    # Initial Setting #
    iter_EM <- 1
    
    prob[[iter]] <- rep(1,ncol(support))/ncol(support)
    
    f_stack = f_prob <- matrix(0, nrow=nr, ncol=ncol(support))
    
    obs.lik_EM <- c()
    
    for (j in 1:ncol(support)) {
      f_stack[,j] <- matrix(dsn(x=data, xi=location[[iter-1]], omega=support["support",j], alpha=skew[[iter-1]]), nrow=nr)
      f_prob[,j] <- prob[[iter]][j] * f_stack[,j]
    }
    
    # 1st Loop #
    
    obs.lik_EM <- obs.loglik(f=f_prob)
    
    error <- 100
    
    # Loop Statement : Find mixing probability by EM-Algorithm #
    while (error > 1e-08) {
      
      iter_EM <- iter_EM+1
      
      denom = rowSums(f_prob)
      
      zhat = f_prob/denom
      
      sum.z <- colSums(zhat)
      
      prob[[iter]] <- sum.z/nr
      
      for (i in 1:ncol(support)){
        
        f_prob[,i] <- prob[[iter]][i] * f_stack[,i]
        
      }
      
      obs.lik_EM <- c(obs.lik_EM, obs.loglik(f=f_prob)) # all(diff(obs.lik_EM)>0)
      
      error <- obs.lik_EM[length(obs.lik_EM)] - obs.lik_EM[length(obs.lik_EM)-1]
      
    }
    
    # Estimate location & skewness parameter #
    
    other_param <- optim(par=c(1,1), fn=find.xi_lambda, prob=prob[[iter]], support=support["support",], data=data, control=list(fnscale=-1))$par
    
    location[[iter]] <- other_param[1] ; skew[[iter]] <- other_param[2]
    
    for (j in 1:ncol(support)) {
      f_stack[,j] <- matrix(dsn(x=data, xi=location[[iter]], omega=support["support",j], alpha=skew[[iter]]), nrow=nr)
      f_prob[,j] <- prob[[iter]][j] * f_stack[,j]
    }
    
    NPMLE <- matrix(rowSums(f_prob), nrow=nr)
    
    # Arbitrary bin numbers of histogram just for reference #
    
    hist(data, nclass=50, freq=FALSE, xlim = c(floor(min(data)),ceiling(max(data))),
         ylab='density',main="Scale Mixture of Skew Normal Distribution")
    points(sort(data), NPMLE[order(data)], type="l", col="blue", lwd=2)
    
    # Theoretical stopping rule : directional derivative using max(local_maxima) = 0 #
    # But, due to computational problems, stopping Rule : Among Local Maxima, Max(Directional Derivative) is less than 1 #
    
    if (iter==2) {given.dist <- matrix(dsn(data, xi=location[[2]], omega=scale[[1]], alpha=skew[[2]]), nrow=nr)}
    if (iter==2) {obs.lik_final <- sum(log(given.dist))} 
    
    # observed log-likelihood of initial distribution
    
    Q_sigma <- matrix(dsn(x=data, xi=location[[iter]], omega=support["support",which.max(support["dir_deriv",])], alpha=skew[[iter]]), nrow=nr)
    
    Dp_value <- c(Dp_value, sum(Q_sigma/given.dist)-nr)
    
    Dp <- Dp_value[iter-1]
    
    obs.lik_final <- c(obs.lik_final, obs.loglik(NPMLE))
    
    print(list(Directional_Derivative=Dp, change_likelihood = obs.lik_final[length(obs.lik_final)]-obs.lik_final[length(obs.lik_final)-1]))
    print(rbind(support, mixing_prob=prob[[iter]]))
    print(rbind(location=location[[iter]], skewness=skew[[iter]]))
    
    given.dist <- NPMLE # run after plotting the directional deriviative
    
  }
  
  scale_parameter <- vector(mode="list", length=iter)
  for (i in 1:iter){
    scale_parameter[[i]] <- scale[[i]]
  }
  
  return(list(Iteration = iter, Directional_Derivative=Dp_value, mixing_probability=prob[[iter]], location_parameter=location[[iter]],
              scale_parameter=unlist(scale), skewness_parameter=skew[[iter]], NPMLE=NPMLE, Observed_Loglik=obs.lik_final, support_step=scale_parameter))
  
}

# Check the function #

res <- snmix.isdm(data=data, grid=20, constant=0.1) # how to use 
res$Iteration 
res$Directional_Derivative
res$location_parameter ; res$skewness_parameter
all(diff(res$Observed_Loglik)>0) # Since the observed log_lik should be non-decreasing, output must be "True".

hist.sn <- histogram(data, type="regular") 
hist(data, breaks=hist.sn$breaks, freq = FALSE, ylim=c(0,.5), main="Histogram : Scale Mixture of Skew Normal")
points(sort(data), res$NPMLE[order(data)], col="blue", lwd=2, type="l", lty=2) # NPMLE of the data
curve(0.6*dsn(x, xi=0, omega=1, alpha=5)+0.4*dsn(x, xi=0, omega=3, alpha=5), n=1000, add=TRUE, col='red', lwd=2) # True density
legend("topright", legend = c("True Density", "NPMLE by ISDM"),
       col = c("red", "blue"), lty=c(1,1), lwd=2, cex=1.3)

# You can use this function when estimating the density of the data, especially skewed data. It also covers symmetric data # 
###########################
# Log-Normal Distribution #
###########################
# (1) Log-Normal with meanlog=0, sdlog=1 #
set.seed(100)
data.ln01 <- rlnorm(1000, meanlog = 0, sdlog = 1)
head(data.ln01)
data.ln01 <- as.matrix(data.ln01, nrow=1)
hist.lnorm01 <- histogram(data.ln01, type="regular")

res.ln01 <- snmix.isdm(data=data.ln01, grid=50, constant=0.1)
res.ln01$Iteration
res.ln01$Directional_Derivative
res.ln01$location_parameter ; res.ln01$skewness_parameter
all(diff(res.ln01$Observed_Loglik)>0)

hist(data.ln01, breaks=hist.lnorm01$breaks, freq=FALSE,xlim=c(0,10), ylim=c(0,0.8), xlab="x",main="lognormal distribution with mu=0, sigma=1")
curve(dlnorm(x, 0,1), add=TRUE, col='red', n=1000, lwd=2)
points(sort(data.ln01), res.ln01$NPMLE[order(data.ln01)], col="blue", lwd=2, type="l", lty=2)
legend("topright", legend = c("True Density", "NPMLE by ISDM"),
       col = c("red", "blue"), lty=c(1,2), lwd=2, cex=1.3)

# (2) Log-Normal with meanlog=0, sdlog=0.5 #
data.ln05 <- rlnorm(1000, meanlog = 0, sdlog = 0.5)
data.ln05 <- as.matrix(data.ln05, nrow=1)
hist.lnorm05 <- histogram(data.ln05, type="regular")

res.ln05<- snmix.isdm(data=data.ln05, grid=50, constant=0.1)
res.ln05$Iteration
res.ln05$Directional_Derivative
res.ln05$location_parameter ; res.ln05$skewness_parameter
all(diff(res.ln05$Observed_Loglik)>0)

hist(data.ln05, breaks=hist.lnorm05$breaks, freq=FALSE, xlim=c(0,5), ylim=c(0,1), xlab="x", main="lognormal distribution with mu=0, sigma=0.5")
curve(dlnorm(x, 0,0.5), add=TRUE, col='red', n=1000, lwd=2)
points(sort(data.ln05), res.ln05$NPMLE[order(data.ln05)], col="blue", lwd=2, type="l", lty=2)
legend("topright", legend = c("True Density", "NPMLE by ISDM"),
       col = c("red", "blue"), lty=c(1,2), lwd=2, cex=1.3)

########################
# Weibull Distribution #
########################

# (1) Weibull distribution with shape=1, scale=1 #
set.seed(100)
data.wei1 <- rweibull(1000, shape = 1, scale = 1)

data.wei1 <- as.matrix(data.wei1, nrow=1)

hist.wei1 <- histogram(data.wei1, type="regular")
curve(dweibull(x, 1, 1), add=TRUE, col='red', n=1000, lwd=2)

res.wei1<- snmix.isdm(data=data.wei1, grid=30, constant=0.2)

res.wei1$Iteration
res.wei1$Directional_Derivative
res.wei1$location_parameter ; res.wei1$skewness_parameter
all(diff(res.wei1$Observed_Loglik)>0)
res.wei1$scale_parameter
res.wei1$mixing_probability

hist(data.wei1, breaks=hist.wei1$breaks, ylim=c(0,1.2), freq = FALSE, xlab="x", main="Weibull Distribution with shape=1, scale=1")
curve(dweibull(x, 1, 1), add=TRUE, col='red', n=1000, lwd=2)
points(sort(data.wei1), res.wei1$NPMLE[order(data.wei1)], col="blue", lwd=2, type="l", lty=2)
legend("topright", legend = c("True Density", "NPMLE by ISDM"),
       col = c("red", "blue"), lty=c(1,2), lwd=2, cex=1.3)

# (2) Weibull distribution with shape=1.5, scale=1 #
set.seed(100)
data.wei15 <- rweibull(1000, shape = 1.5, scale = 1)

data.wei15 <- as.matrix(data.wei15, nrow=1)

hist.wei15 <- histogram(data.wei15, type="regular")
curve(dweibull(x, 1.5, 1), add=TRUE, col='red', n=1000, lwd=2)


res.wei15<- snmix.isdm(data=data.wei15, grid=50, constant=0.1)

res.wei15$Iteration
res.wei15$Directional_Derivative
res.wei15$location_parameter ; res.wei15$skewness_parameter
all(diff(res.wei15$Observed_Loglik)>0)
res.wei15$scale_parameter
res.wei15$mixing_probability

hist(data.wei15, breaks=hist.wei15$breaks, ylim=c(0,1),freq = FALSE, xlab="x", main="Weibull Distribution with shape=1.5, scale=1")
curve(dweibull(x, 1.5, 1), add=TRUE, col='red', n=1000, lwd=2)
points(sort(data.wei15), res.wei15$NPMLE[order(data.wei15)], col="blue", lwd=2, type="l", lty=2)
legend("topright", legend = c("True Density", "NPMLE by ISDM"),
       col = c("red", "blue"), lty=c(1,2), lwd=2, cex=1.3)

##################
# t Distribution #
##################

# (1) t-distribution with d.f. 3 #

set.seed(100)
data.t3 <- rt(1000, 3)
data.t3 <- as.matrix(data.t3, nrow=1)
hist.t3 <- histogram(data.t3, type="regular")

res.t3<- snmix.isdm(data=data.t3, grid=50, constant=0.1)
res.t3$Iteration
res.t3$Directional_Derivative
res.t3$location_parameter ; res.t3$skewness_parameter
all(diff(res.t3$Observed_Loglik)>0)
res.t3$scale_parameter
res.t3$mixing_probability
hist(data.t3, breaks=hist.t3$breaks, ylim=c(0,.4),freq = FALSE, xlab="x",main="Student T Distribution with d.f. 3")
curve(dt(x, 3), add=TRUE, col='red', n=1000, lwd=2)
points(sort(data.t3), res.t3$NPMLE[order(data.t3)], col="blue", lwd=2, type="l", lty=2)
legend("topleft", legend = c("True Density", "NPMLE by ISDM"),
       col = c("red", "blue"), lty=c(1,2), lwd=2, cex=1.3)

########################
# Laplace Distribution #
########################

# (1) Laplace distribution with location=0, scale=1 #

library(LaplacesDemon)
set.seed(100)

data.lap <- rlaplace(1000, location=0, scale=1)
data.lap <- as.matrix(data.lap, nrow=1)
hist.lap <- histogram(data.lap, type="regular") 

res.lap1<- snmix.isdm(data=data.lap, grid=50, constant=0.1)
res.lap1$Iteration
res.lap1$Directional_Derivative
res.lap1$location_parameter ; res.lap1$skewness_parameter
all(diff(res.lap1$Observed_Loglik)>0)
res.lap1$scale_parameter
res.lap1$mixing_probability

hist(data.lap, breaks=hist.lap$breaks, ylim=c(0,.5),freq = FALSE, xlab="x", main="Laplace Distribution with location = 0, scale =1")
curve(dlaplace(x, 0, 1), add=TRUE, col='red', n=1000, lwd=2)
points(sort(data.lap), res.lap1$NPMLE[order(data.lap)], col="blue", lwd=2, type="l", lty=2)
legend("topright", legend = c("True Density", "NPMLE by ISDM"),
       col = c("red", "blue"), lty=c(1,2), lwd=2, cex=1.3)
