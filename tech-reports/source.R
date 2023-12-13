#############

# establishes the interpolated empirical CDF
Ftilde <- function(y, t, ystar, component_name){
  y <- sort(y)
  n <- length(y)
  ytilde <- rep(0, n + 1)
  
  if (component_name == 'bWAR' | component_name == 'fWAR' | component_name == 'ERA') {
    ytilde[1] <- y[1] - (y[2] - y[1])
  }
  if (component_name == 'HR'| component_name == 'BB') {
    # since the minimal HR is greater or equal to 0.
    ytilde[1] <- 0
  }
  if (component_name == 'bWAR_p' | component_name == 'fWAR_p') {
    ytilde[1] <- y[1] - (y[2] - y[1])/10
  }
  if (component_name == 'SO' | component_name == 'AVG') {
    ytilde[1] <- ifelse(y[1] - (y[2] - y[1]) < 0, 0, y[1] - (y[2] - y[1]) )
  }
  
  ytilde[n+1] <- y[n] + ystar
  ytilde[2:n] <- unlist(lapply(2:n, function(j){
    (y[j]+y[j-1])/2 
  }))
  
  if (t >= ytilde[n+1]) {
    1 - 0.1^7
  } else if (t <= ytilde[1]) {
    0
  } else {
    j <- length(which(ytilde < t))
    (j - 1) / n + (t - ytilde[j]) / (n*(ytilde[j+1] - ytilde[j]))
  }
  
}

#############

#############

# computes the ystar in the empirical CDF for each season
thresh_fun <- function(component, component_name){
  ## stab means we add a small number to the largest value to avoid numerical stability problems. 
  ## stab depends on the scale of input values. In most cases, the default is 0.01. 
  ## cutoff means we select a certain percentage of systems that include maximal possible components in the tail
  ## rather than include k possible components in the tail based on the adjusted R adjusted square. 
  ## If the distance between the largest value and second largest value in the right tail is relatively large,
  ## adjusted R square may not be a good quantity to represent how well the fit is. 
  ## Then we usually choose the one sixth of all systems that their distances between the largest value and 
  ## second largest value are in the top one sixth. 
  
  if (component_name == 'bWAR') {
    #cutoff <- 1.5e-2 # raw comp
    # cutoff <- 1.17e-2 #BB
    cutoff <- 1.4e-2 # schell
    stab <- 0.01
  }
  if (component_name == 'fWAR') {
    #cutoff <- 1.75e-2
    cutoff <- 1.45e-2 # schell
    stab <- 0.01
  }
  if (component_name == 'HR') {
    cutoff <- 1.70e-2 # schell
    #cutoff <- 2.06e-2
    stab <- 0.01
  }
  if (component_name == 'BB') {
    #cutoff <- 4.05e-2
    cutoff <- 3.30e-2 # schell
    stab <- 0.027
  }
  
  if (component_name == 'AVG') {
    #cutoff <- 2.35e-2
    cutoff <- 2.30e-2 # schell
    stab <- 0.01
  }
  if (component_name == 'ERA') {
    # cutoff <- 0.47 # raw comp
    cutoff <- 0.4 # Schell
    stab <- 0.2
  }
  if (component_name == 'bWAR_p') {
    # cutoff <- 7.55e-3 # raw comp
    # cutoff <- 6e-3 # BB
    cutoff <- 7.09e-3 # Schell
    stab <- 0.01
  }
  if (component_name == 'fWAR_p') {
    #cutoff <- 5.7e-3 # raw comp
    cutoff <- 5.124e-3 # Schell
    stab <- 0.01
  }
  if (component_name == 'SO') {
    # cutoff <- 1.4 # raw comp
    cutoff <- 1.31 # Schell
    stab <- 1
  }
  # obtain initial quantities for linear approximation
  Y <- sort(as.matrix(component))
  n <- length(Y)
  Y[n] <- Y[n] + stab # for stability
  pi <- 1 - (n:1 - 1/3)/(n + 1/3)
  W <- log(pi/(1-pi))
  K1 = max(6, floor(1.3*sqrt(n))); K2 = 2*floor(log10(n)*sqrt(n))
  k <- 6
  
  # use arguments from Scholz section 3 for estimating k
  #
  # this argument is based on model fit and not longest stretch of 
  # contiguous I0
  ind <- NULL
  try({
    k_selector <- do.call(rbind, lapply(6:K2, function(k){
      
      Ytil <- Y - median(Y)
      Ztil <- tail(Ytil, k)
      M1k <- 1/(k-1) * sum( log(Ztil[2:k]/Ztil[1]) )
      M2k <- 1/(k-1) * sum( log(Ztil[2:k]/Ztil[1])^2 )
      ck <- M1k + 1 - 0.5*(1 - M1k^2/M2k)^{-1}
      fck <- ((-n*log(pi))^{-ck} - 1)/ck
      
      Sigma <- matrix(0, k, k)
      for(i in 1:k){
        for(j in 1:i){
          Sigma[i,j] <- i^{-ck-1} * j^{-ck}
        } 
      }
      for(j in 1:k){
        for(i in 1:(j-1)){
          Sigma[i,j] <- j^{-ck-1} * i^{-ck}
        } 
      }
      
      rotate <- function(x) t(apply(x, 2, rev))
      Sigma <- rotate(rotate(Sigma))
      Sigma.inv <-  solve(Sigma)
      eig <- eigen(Sigma.inv)
      C <- eig$vec %*% diag(sqrt(eig$val)) %*% t(eig$vec)
      Zk <- C %*% tail(Y, k)
      Xk <- cbind(1, tail(fck, k))
      Wk <-  C %*% Xk
      # try linear and quadratic model
      m1 <- lm(tail(Y, k) ~ tail(fck, k))
      m2 <- lm(tail(Y, k) ~ tail(fck, k) + I(tail(fck, k)^2))
      m3 <- lm(Zk ~ -1 + Wk)
      delta.sq <- summary(m3)$sigma^2
      Tk <- coef(m3)[2] / summary(m3)$sigma
      
      kappa.sq <- solve(crossprod(Wk))[2,2]
      kappa <- sqrt(kappa.sq)
      I0 <- c(kappa * qt(0.25, df = k - 2, ncp = 1/kappa),
              kappa * qt(0.75, df = k - 2, ncp = 1/kappa))
      I1 <- c(kappa * qt(0.05, df = k - 2, ncp = 1/kappa), 
              kappa * qt(0.95, df = k - 2, ncp = 1/kappa))
      I0int <- ifelse(I0[1] <= Tk && Tk <= I0[2], 1, 0)
      I1int <- ifelse(I1[1] <= Tk && Tk <= I1[2], 1, 0)
      c(k, Tk, I0int, I1int, summary(m1)$adj.r.squared, 
        summary(m2)$adj.r.squared)
      
    }))
    
    #k <- k_selector[max(which(k_selector[, 3] == 1)), 1]
    #k <- k_selector[which.max(k_selector[, 5]), 1]
    k_selector <- as.data.frame(k_selector)
    colnames(k_selector) <- c("k", "Tk", "I0", "I1", "R.sq", "Rquad.sq")
    k_selector_I0 <- k_selector %>% filter(I0 == 1)
    a <- which.max(k_selector_I0$R.sq)
    b <- which.max(k_selector_I0$Rquad.sq)
    ind <- which.max(c(k_selector_I0[a, ]$R.sq, 
                       k_selector_I0[b, ]$Rquad.sq))
    k <- k_selector_I0[c(a,b)[ind] , 1]
    if(diff(Y)[n-1] > cutoff){ 
      k <- max(k_selector_I0$k)
      if(k < 0) k <- K2
    }
    
  }, silent = TRUE)
  
  if(length(k) == 0) k <- round(mean(K1,K2))
  if(is.na(k)) k <- round(mean(K1,K2))
  if(k == 0) k <- round(mean(K1,K2))
  if(k >= n) k <- K2
  
  
  # find probability value using linear tail behavior
  Z <- tail(Y, k)
  m1 <- lm(tail(Y, k) ~ tail(pi, k))
  beta <- m1$coefficients
  ystar <- ub <- 0
  f <- function(x) beta[1] + beta[2] * x - max(Y)
  #delta <- beta[2]
  try({
    foo <- uniroot(f, c(0.0001, 5), tol = 1e-10)
    ub <- foo$root        
  })
  
  # find probability value using logistic tail behavior
  if(ub >= 1){
    m1 <- lm(tail(Y,k) ~ tail(W, k))
    beta <- m1$coefficients
    f <- function(x) beta[1] + beta[2] * log(x/(1-x)) - max(Y)
    try({
      foo <- uniroot(f, c(0.000001, 0.999999), tol = 1e-10)
      ub <- foo$root        
    })  
  }
  
  # if possible, find ystar by tying logistic behavior argument to 
  # our Ftilde function
  if(ub >= Ftilde(y = Y, t = max(Y), ystar = 10, component_name = component_name)){
    try({
      g <- function(ystar) ub - Ftilde(y = Y, t = max(Y), ystar = ystar, component_name = component_name)
      bar <- uniroot(g, c(0, 10), tol = 1e-10)
      ystar <- bar$root
    })
  }
  
  # if the above is not possible, try a similar approach for different 
  # suitable values of k. 
  #
  # The above fails because ub < Ftilde(y = Y, t = max(Y), ystar = 10, component_name = component_name) 
  # suggesting that the largest achiever is performaing much worse than 
  # expected. Thus ystar should be "large". A default large value will 
  # be ystar = 6 (altered to be log(1 + 6) for stability). This will 
  # be used when all else fails.
  flag <- NULL
  if(ub < Ftilde(y = Y, t = max(Y), ystar = 10, component_name = component_name)){
    
    # first try for largest suitable k as dictated by Scholz Section 3
    k <- max(k_selector_I0$k) 
    m1 <- lm(tail(Y,k) ~ tail(W, k) + I(tail(W, k)^2))
    beta <- m1$coefficients
    f <- function(x) beta[1] + beta[2] * log(x/(1-x)) +
      beta[3] * log(x/(1-x))^2 - max(Y)
    flag <- try({
      foo <- uniroot(f, c(0.0001, 0.9999), tol = 1e-10)
      ub <- foo$root
    }, silent = TRUE)
    while(class(flag) == "try-error"){
      k <- k - 1
      # method fails; use ystar = 4
      if(k < 6){
        ystar <- 6
        break
      }
      m1 <- lm(tail(Y,k) ~ tail(W, k) + I(tail(W, k)^2))
      beta <- m1$coefficients
      f <- function(x) beta[1] + beta[2] * log(x/(1-x)) +
        beta[3] * log(x/(1-x))^2 - max(Y)
      flag <- try({
        foo <- uniroot(f, c(0.0001, 0.9999), tol = 1e-10)
        ub <- foo$root
      }, silent = TRUE)
    }
    
    ystar_1 <- NULL
    try({
      g <- function(ystar) ub - Ftilde(y = Y, t = max(Y), ystar = ystar, component_name = component_name)
      bar <- uniroot(g, c(0, 10), tol = 1e-10)
      ystar_1 <- bar$root    	
    }, silent = TRUE)
    if(length(ystar_1) == 0) ystar_1 <- 6
    
    
    # now try for smallest suitable k as dictated by Scholz Section 3
    k <- min(k_selector_I0$k)
    if(length(k) == 0) k <- 6
    m1 <- lm(tail(Y,k) ~ tail(W, k) + I(tail(W, k)^2))
    beta <- m1$coefficients
    f <- function(x) beta[1] + beta[2] * log(x/(1-x)) + 
      beta[3] * log(x/(1-x))^2 - max(Y)
    flag <- try({
      foo <- uniroot(f, c(0.0001, 0.9999), tol = 1e-10)
      ub <- foo$root        
    }, silent = TRUE)  
    while(class(flag) == "try-error"){
      k <- k + 1
      # method fails; use ystar = 4
      if(k > max(k_selector_I0$k)){
        ystar <- 6
        break
      }
      m1 <- lm(tail(Y,k) ~ tail(W, k) + I(tail(W, k)^2))
      beta <- m1$coefficients
      f <- function(x) beta[1] + beta[2] * log(x/(1-x)) + 
        beta[3] * log(x/(1-x))^2 - max(Y)
      flag <- try({
        foo <- uniroot(f, c(0.0001, 0.9999), tol = 1e-10)
        ub <- foo$root        
      }, silent = TRUE)  
    }
    
    ystar_2 <- NULL
    try({
      g <- function(ystar) ub - Ftilde(y = Y, t = max(Y), ystar = ystar, component_name = component_name)
      bar <- uniroot(g, c(0, 10), tol = 1e-10)
      ystar_2 <- bar$root			
    }, silent = TRUE)
    if(length(ystar_2) == 0) ystar_2 <- 6
    
    # take ystar as the average of the lowest working k and 
    # largest working k
    ystar <- mean(c(ystar_1, ystar_2))
    
  }
  
  # if changing k does not work, then try throwing out extreme 
  # observations and computing ystar for the reduced sample (Ytil)
  #
  # then compute ystar = max(Y) - max(Ytil) + ystar*
  #
  # where ystar* is computed with respect to Ytil
  if(ystar == 6){
    k <- k_selector_I0[c(a,b)[ind] , 1]
    if(diff(Y)[n-1] > cutoff){ 
      k <- max(k_selector_I0$k)
      if(k < 0) k <- K2
    }
    if(length(k) == 0) k <- round(mean(K1,K2))
    if(is.na(k)) k <- round(mean(K1,K2))
    if(k == 0) k <- round(mean(K1,K2))
    if(k >= n) k <- K2
    
    m1 <- lm(tail(Y, k) ~ tail(W, k) + I(tail(W, k)^2))
    beta <- m1$coefficients
    f <- function(x) beta[1] + beta[2] * log(x/(1-x)) + 
      beta[3] * log(x/(1-x))^2 - max(Y)
    flag <- flag2 <- try({
      foo <- uniroot(f, c(0.0001, 0.9999), tol = 1e-10)
      ub <- foo$root        
    }, silent = TRUE)
    
    n_lwr <- n 
    Ytil <- Y
    Xtil <- 1 - (n_lwr:1 - 1/3)/(n_lwr + 1/3)
    Wtil <- log(Xtil/(1-Xtil))
    while(class(flag) == "try-error" | class(flag2) == "try-error"){
      Ytil <- Ytil[-n_lwr]
      if(any(tail(Ytil, k) < 0)){
        ystar <- 6
        break
      }
      n_lwr <- n_lwr - 1
      Xtil <- 1 - (n_lwr:1 - 1/3)/(n_lwr + 1/3)
      Wtil <- log(Xtil/(1-Xtil))
      m2 <- lm(tail(Ytil, k) ~ tail(Wtil, k) + I(tail(Wtil, k)^2))
      beta <- m2$coefficients
      f <- function(x) beta[1] + beta[2] * log(x/(1-x)) + 
        beta[3] * log(x/(1-x))^2 - max(Ytil)
      flag <- try({
        foo <- uniroot(f, c(0.0001, 0.9999), tol = 1e-10)
        ub <- foo$root        
      }, silent = TRUE)
      flag2 <- try({
        g <- function(ystar) ub - Ftilde(y = Ytil, t = max(Ytil), ystar = ystar, component_name = component_name)
        bar <- uniroot(g, c(0, 10), tol = 1e-10)
        ystar <- bar$root
      }, silent = TRUE)
    }
    ystar <- max(Y) - Y[n_lwr] + ystar
    
  }
  
  # for stability
  ystar
  
}

#############

#############

# computes their latent talents using non-parametric distribution measuring the components
talent_computing_nonpara <- function(dataset, component_name, year, ystar, alpha){
  
  Ftilde <- function(y, t, ystar, component_name){
    y <- sort(y)
    n <- length(y)
    ytilde <- rep(0, n + 1)
    
    if (component_name == 'bWAR' | component_name == 'fWAR' | component_name == 'ERA') {
      ytilde[1] <- y[1] - (y[2] - y[1])
    }
    if (component_name == 'HR'| component_name == 'BB') {
      # since the minimal HR is greater or equal to 0.
      ytilde[1] <- 0
    }
    if (component_name == 'bWAR_p' | component_name == 'fWAR_p') {
      ytilde[1] <- y[1] - (y[2] - y[1])/10
    }
    if (component_name == 'SO' | component_name == 'AVG') {
      ytilde[1] <- ifelse(y[1] - (y[2] - y[1]) < 0, 0, y[1] - (y[2] - y[1]) )
    }
    
    ytilde[n+1] <- y[n] + ystar
    ytilde[2:n] <- unlist(lapply(2:n, function(j){
      (y[j]+y[j-1])/2
    }))
    
    if (t >= ytilde[n+1]) {
      1 - 0.1^7
    } else if (t <= ytilde[1]) {
      0
    } else {
      j <- length(which(ytilde < t))
      (j - 1) / n + (t - ytilde[j]) / (n*(ytilde[j+1] - ytilde[j]))
    }
    
  }
  
  Aptitude_nonpara <- function(p, alpha, npop){
    
    # converts order stats to their percentiles
    order_pbino <- function(p = 0, k = 1, n = 1e4){
      pbinom(k - 1, prob = p, size = n, lower.tail = FALSE)
    }
    
    # converts a vector of order stats 
    # to their percentiles. This vector should be the entire 
    # sample sorted in increasing order
    p <- sort(p) # just in case
    n <- length(p)
    u = unlist(lapply(1:n, function(j){
      #order_pbino(p[j], k = 251-n+j, n = 251)
      order_pbino(p[j], k = j, n = n)
    }))
    
    # transforms percentiles from order stats (in increasing order)
    # to Pareto values corresponding to the general population 
    # of a greater than or equal to size
    # default alpha is that of the Pareto principle 80-20
    n <- length(u)
    if(length(npop) == 1) npop <- rep(npop, n)
    unlist(lapply(1:n, function(j){
      qPareto(qbeta(u[j], j + npop[j] -n , n + 1 - j), t = 1, alpha = alpha)
      #qPareto(qbeta(u[j], j + npop[j] -n + 251-n , n + 1 - j), t = 1, alpha = alpha)
    }))
  }
  
  foo <- dataset %>% filter(yearID == year) %>% 
    arrange(comp) 
  bar <- foo %>% filter(full_time == 'Y')
  full_comp <- bar$comp
  ## batter WAR talent
  bar <- bar %>% 
    mutate(WAR_talent = 
             Aptitude_nonpara(p = unlist(lapply(comp, function(xx) 
               Ftilde(y = comp, t = xx, ystar = ystar, component_name = component_name))), alpha = alpha, npop = pops))
  
  max_WAR_talent <- max(bar$WAR_talent) - 1
  range <- which(!(foo$playerID %in% bar$playerID))
  
  ## using the distribution from full time players
  bar <- rbind(bar, do.call(rbind, lapply(range, function(j){
    rbind(bar %>% dplyr::select(-WAR_talent), foo[j, ]) %>% arrange(comp) %>%
      mutate(WAR_talent = Aptitude_nonpara(p = unlist(lapply(comp, function(xx) 
        Ftilde(y = full_comp, t = xx, ystar = ystar, component_name = component_name))), alpha = alpha, npop = pops)) %>%
      filter(full_time == 'N') %>% 
      mutate(WAR_talent = ifelse(WAR_talent > max_WAR_talent+1, max_WAR_talent, WAR_talent))
  })))
  bar %>% mutate(ystar = ystar)
}

#############

#############

# maps their latent talents to a common lab environment and get the era-adjusted statistics
# using non-parametric distribution measuring the components
career_talent_nonpara <- function(dataset, component_name, snippet, alpha){
  Rev_Aptitude_nonpara <- function(x, ytilde, alpha, npop){
    
    # transforms ordered Pareto values corresponding to
    # the general population to percentiles from order stats
    # (in increasing order)
    n <- length(x)
    if(length(npop) == 1) npop <- rep(npop, n)
    u = unlist(lapply(1:n, function(j){
      pbeta(pPareto(x[j], t = 1, alpha = alpha), j + npop[j]-n, n + 1 - j)
    }))
    
    
    ## map the quantile to the a predicated sample value
    map_Y <- function(u, ytilde){
      n <- length(ytilde)-1
      seqence <- seq(0, 1, 1/n)
      pos <- findInterval(u, seqence)
      out <- (n*u -pos + 1) * (ytilde[(pos+1)] - ytilde[pos]) + ytilde[pos]
      return(out)
    }
    
    ## map the vector of quantiles to the predicated sample values
    n <- length(u)
    a <- qbeta(u, shape1 = 1:n, shape2 = n:1)
    out <- sapply(1:n, function(x) map_Y(a[x], ytilde = ytilde))
    out
    
  }
  
  talent_all <- sort(c(talent_new, snippet$WAR_talent))
  index_snippet <- max(which(talent_all == snippet$WAR_talent))
  
  snippet %>% 
    mutate(adj_comp = Rev_Aptitude_nonpara(talent_all, 
                                           ytilde = ytilde, alpha = alpha, 
                                           npop = pop_new)[index_snippet])
  
}

#############

#############

# computes their latent talents using parametric distribution measuring the components
talent_computing_para <- function(dataset, year, alpha){
  # dataset should include the component name and corresponding npop 
  Aptitude_para <- function(y, mean = 0, sd = 1, alpha, npop){
    
    # converts normal order stats to their percentiles
    order_pnorm <- function(q = 0, mean = 0, sd = 1, k = 1, n = 1e4){
      p <- pnorm(q, mean = mean, sd = sd)
      pbinom(k - 1, prob = p, size = n, lower.tail = FALSE)
    }
    
    # converts a vector of normal order stats 
    # to their percentiles. This vector should be the entire 
    # sample sorted in increasing order
    q <- sort(y) # just in case
    n <- length(y)
    u = unlist(lapply(1:n, function(j){
      order_pnorm(y[j], k = j, n = n, mean = mean, sd = sd)
    }))
    
    # transforms percentiles from order stats (in increasing order)
    # to Pareto values corresponding to the general population 
    # of a greater than or equal to size
    # default alpha is that of the Pareto principle 80-20
    n <- length(u)
    if(length(npop) == 1) npop <- rep(npop, n)
    unlist(lapply(1:n, function(j){
      qPareto(qbeta(u[j], j + npop[j]-n, n + 1 - j), t = 1, alpha = alpha)
    }))
  }
  
  # consider the full time batters
  foo <- dataset %>% filter(yearID == year) %>% arrange(comp)
  ## average talent
  bar <- foo %>% filter(full_time == 'Y')
  nsys <- nrow(bar)
  mu_comp <- mean(bar$comp)
  sd_comp <- sd(bar$comp)
  # components to build the CDF using parametric distribution
  full_comp <- bar$comp
  bar <- bar %>% mutate(talent = 
                          Aptitude_para(full_comp, mean = mu_comp, sd = sd_comp, alpha = alpha, npop = pops))
  max_talent <- max(bar$talent) - 1
  # consider the non-full time batters
  non_full <- foo %>% filter(full_time == 'N')
  bar <- rbind(bar, do.call(rbind, lapply(1:nrow(non_full), function(j){
    m <- rbind(bar %>% dplyr::select(-talent), non_full[j, ]) %>% arrange(comp)
    m %>% mutate(talent = Aptitude_para(comp, mean = mu_comp, sd = sd_comp, alpha = alpha, npop = pops)) %>%
      filter(full_time == 'N') %>% 
      mutate(talent = ifelse(talent > max_talent+1, max_talent, talent))
  })))
  bar  
}

#############

#############

# maps their latent talents to a common lab environment and get the era-adjusted statistics
# using parametric distribution measuring the components
career_talent_para <- function(dataset, snippet, alpha){
  Rev_Aptitude_para <- function(x, mean = 0, sd = 1, alpha, npop){
    
    # transforms ordered Pareto values corresponding to 
    # the general population to percentiles from order stats 
    # (in increasing order)
    n <- length(x)
    if(length(npop) == 1) npop <- rep(npop, n)  
    u = unlist(lapply(1:n, function(j){
      pbeta(pPareto(x[j], t = 1, alpha = alpha), j + npop[j]-n, n + 1 - j)
    }))
    
    # take the ordered percentiles and convert them to order statistics 
    # from a normal distribution
    n <- length(u)
    qnorm(qbeta(u, shape1 = 1:n, shape2 = n:1), mean = mean, sd = sd)
    
  }
  
  talent_all <- sort(c(talent_new, snippet$talent))
  index_snippet <- max(which(talent_all == snippet$talent))
  
  snippet %>% 
    mutate(adj_comp = Rev_Aptitude_para(talent_all, mean = mean_new, alpha = alpha, sd = sd_new, 
                                        npop = pop_new)[index_snippet])
  
}

#############

#############

# computes their latent talents using non-parametric distribution measuring the components
talent_computing_nonpara_fixed <- function(dataset, component_name, year, ystar, alpha){
  
  Ftilde_fixed <- function(y, t, ystar, component_name){
    y <- sort(y)
    n <- length(y)
    ytilde <- rep(0, n + 1)
    
    if (component_name == 'bWAR' | component_name == 'fWAR' | component_name == 'ERA') {
      ytilde[1] <- y[1] - (y[2] - y[1])
    }
    if (component_name == 'HR'| component_name == 'BB') {
      # since the minimal HR is greater or equal to 0.
      ytilde[1] <- 0
    }
    if (component_name == 'bWAR_p' | component_name == 'fWAR_p') {
      ytilde[1] <- y[1] - (y[2] - y[1])/10
    }
    if (component_name == 'SO' | component_name == 'AVG') {
      ytilde[1] <- ifelse(y[1] - (y[2] - y[1]) < 0, 0, y[1] - (y[2] - y[1]) )
    }
    
    ytilde[n+1] <- y[n] + ystar
    ytilde[2:n] <- unlist(lapply(2:n, function(j){
      (y[j]+y[j-1])/2
    }))
    
    if (t >= ytilde[n+1]) {
      1 - 0.1^7
    } else if (t <= ytilde[1]) {
      1 - n/251
    } else {
      j <- length(which(ytilde < t))
      (250 - n + j) / 251 + (t - ytilde[j]) / (251*(ytilde[j+1] - ytilde[j]))
    }
    
  }
  
  Aptitude_nonpara <- function(p, alpha, npop){
    
    # converts order stats to their percentiles
    order_pbino <- function(p = 0, k = 1, n = 1e4){
      pbinom(k - 1, prob = p, size = n, lower.tail = FALSE)
    }
    
    # converts a vector of order stats 
    # to their percentiles. This vector should be the entire 
    # sample sorted in increasing order
    p <- sort(p) # just in case
    n <- length(p)
    u = unlist(lapply(1:n, function(j){
      order_pbino(p[j], k = 251-n+j, n = 251)
      #order_pbino(p[j], k = j, n = n)
    }))
    
    # transforms percentiles from order stats (in increasing order)
    # to Pareto values corresponding to the general population 
    # of a greater than or equal to size
    # default alpha is that of the Pareto principle 80-20
    n <- length(u)
    if(length(npop) == 1) npop <- rep(npop, n)
    unlist(lapply(1:n, function(j){
      qPareto(qbeta(u[j], j + npop[j] -n , n + 1 - j), t = 1, alpha = alpha)
      #qPareto(qbeta(u[j], j + npop[j] -n + 251-n , n + 1 - j), t = 1, alpha = alpha)
    }))
  }
  
  foo <- dataset %>% filter(yearID == year) %>% 
    arrange(comp) 
  bar <- foo %>% filter(full_time == 'Y')
  full_comp <- bar$comp
  ## batter WAR talent
  bar <- bar %>% 
    mutate(WAR_talent = 
             Aptitude_nonpara(p = unlist(lapply(comp, function(xx) 
               Ftilde_fixed(y = comp, t = xx, ystar = ystar, component_name = component_name))), alpha = alpha, npop = pops))
  
  max_WAR_talent <- max(bar$WAR_talent) - 1
  range <- which(!(foo$playerID %in% bar$playerID))
  
  ## using the distribution from full time players
  bar <- rbind(bar, do.call(rbind, lapply(range, function(j){
    rbind(bar %>% dplyr::select(-WAR_talent), foo[j, ]) %>% arrange(comp) %>%
      mutate(WAR_talent = Aptitude_nonpara(p = unlist(lapply(comp, function(xx) 
        Ftilde_fixed(y = full_comp, t = xx, ystar = ystar, component_name = component_name))), alpha = alpha, npop = pops)) %>%
      filter(full_time == 'N') %>% 
      mutate(WAR_talent = ifelse(WAR_talent > max_WAR_talent+1, max_WAR_talent, WAR_talent))
  })))
  bar %>% mutate(ystar = ystar)
}

#############

#############

# computes their latent talents using parametric distribution measuring the components
talent_computing_para_fixed <- function(dataset, year, alpha){
  # dataset should include the component name and corresponding npop 
  Aptitude_para <- function(y, mean = 0, sd = 1, alpha, npop){
    
    # converts normal order stats to their percentiles
    order_pnorm <- function(q = 0, mean = 0, sd = 1, k = 1, n = 1e4){
      p <- pnorm(q, mean = mean, sd = sd)
      pbinom(k - 1, prob = p, size = n, lower.tail = FALSE)
    }
    
    # converts a vector of normal order stats 
    # to their percentiles. This vector should be the entire 
    # sample sorted in increasing order
    q <- sort(y) # just in case
    n <- length(y)
    u = unlist(lapply(1:n, function(j){
      order_pnorm(y[j], k = 251 - n + j, n = 251, mean = mean, sd = sd)
    }))
    
    # transforms percentiles from order stats (in increasing order)
    # to Pareto values corresponding to the general population 
    # of a greater than or equal to size
    # default alpha is that of the Pareto principle 80-20
    n <- length(u)
    if(length(npop) == 1) npop <- rep(npop, n)
    unlist(lapply(1:n, function(j){
      qPareto(qbeta(u[j], j + npop[j]-n, n + 1 - j), t = 1, alpha = alpha)
    }))
  }
  
  # consider the full time batters
  foo <- dataset %>% filter(yearID == year) %>% arrange(comp)
  ## average talent
  bar <- foo %>% filter(full_time == 'Y')
  nsys <- nrow(bar)
  mu_comp <- mean(bar$comp)
  sd_comp <- sd(bar$comp)
  mu_fixed <- min(bar$comp) + nsys * (mu_comp - min(bar$comp)) / 251
  sd_fixed <- sqrt((sum((bar$comp - mu_fixed)^2) + (nsys * (mu_comp - min(bar$comp)) / 251)^2*(251-nsys))/250)
  # components to build the CDF using parametric distribution
  full_comp <- bar$comp
  bar <- bar %>% mutate(talent = 
                          Aptitude_para(full_comp, mean = mu_fixed, sd = sd_fixed, alpha = alpha, npop = pops))
  max_talent <- max(bar$talent) - 1
  # consider the non-full time batters
  non_full <- foo %>% filter(full_time == 'N')
  bar <- rbind(bar, do.call(rbind, lapply(1:nrow(non_full), function(j){
    m <- rbind(bar %>% dplyr::select(-talent), non_full[j, ]) %>% arrange(comp)
    m %>% mutate(talent = Aptitude_para(comp, mean = mu_fixed, sd = sd_fixed, alpha = alpha, npop = pops)) %>%
      filter(full_time == 'N') %>% 
      mutate(talent = ifelse(talent > max_talent+1, max_talent, talent))
  })))
  bar  
}

#############

#############

# computes their latent talents using non-parametric distribution measuring the components
# the talent follows folded normal distribution
talent_computing_nonpara_folded <- function(dataset, component_name, year, ystar){
  
  Ftilde <- function(y, t, ystar, component_name){
    y <- sort(y)
    n <- length(y)
    ytilde <- rep(0, n + 1)
    
    if (component_name == 'bWAR' | component_name == 'fWAR' | component_name == 'ERA') {
      ytilde[1] <- y[1] - (y[2] - y[1])
    }
    if (component_name == 'HR'| component_name == 'BB') {
      # since the minimal HR is greater or equal to 0.
      ytilde[1] <- 0
    }
    if (component_name == 'bWAR_p' | component_name == 'fWAR_p') {
      ytilde[1] <- y[1] - (y[2] - y[1])/10
    }
    if (component_name == 'SO' | component_name == 'AVG') {
      ytilde[1] <- ifelse(y[1] - (y[2] - y[1]) < 0, 0, y[1] - (y[2] - y[1]) )
    }
    
    ytilde[n+1] <- y[n] + ystar
    ytilde[2:n] <- unlist(lapply(2:n, function(j){
      (y[j]+y[j-1])/2
    }))
    
    if (t >= ytilde[n+1]) {
      1 - 0.1^7
    } else if (t <= ytilde[1]) {
      0
    } else {
      j <- length(which(ytilde < t))
      (j - 1) / n + (t - ytilde[j]) / (n*(ytilde[j+1] - ytilde[j]))
    }
    
  }
  
  Aptitude_nonpara <- function(p, npop){
    
    # converts order stats to their percentiles
    order_pbino <- function(p = 0, k = 1, n = 1e4){
      pbinom(k - 1, prob = p, size = n, lower.tail = FALSE)
    }
    
    # converts a vector of order stats 
    # to their percentiles. This vector should be the entire 
    # sample sorted in increasing order
    p <- sort(p) # just in case
    n <- length(p)
    u = unlist(lapply(1:n, function(j){
      #order_pbino(p[j], k = 251-n+j, n = 251)
      order_pbino(p[j], k = j, n = n)
    }))
    
    # transforms percentiles from order stats (in increasing order)
    # to Pareto values corresponding to the general population 
    # of a greater than or equal to size
    # default alpha is that of the Pareto principle 80-20
    n <- length(u)
    if(length(npop) == 1) npop <- rep(npop, n)
    unlist(lapply(1:n, function(j){
      qfoldnorm(qbeta(u[j], j + npop[j] -n , n + 1 - j),mean = 0, sd = 1)
    }))
  }
  
  foo <- dataset %>% filter(yearID == year) %>% 
    arrange(comp) 
  bar <- foo %>% filter(full_time == 'Y')
  full_comp <- bar$comp
  ## batter WAR talent
  bar <- bar %>% 
    mutate(WAR_talent = 
             Aptitude_nonpara(p = unlist(lapply(comp, function(xx) 
               Ftilde(y = comp, t = xx, ystar = ystar, component_name = component_name))), npop = pops))
  
  max_WAR_talent <- max(bar$WAR_talent) - 1
  range <- which(!(foo$playerID %in% bar$playerID))
  
  ## using the distribution from full time players
  bar <- rbind(bar, do.call(rbind, lapply(range, function(j){
    rbind(bar %>% dplyr::select(-WAR_talent), foo[j, ]) %>% arrange(comp) %>%
      mutate(WAR_talent = Aptitude_nonpara(p = unlist(lapply(comp, function(xx) 
        Ftilde(y = full_comp, t = xx, ystar = ystar, component_name = component_name))), npop = pops)) %>%
      filter(full_time == 'N') %>% 
      mutate(WAR_talent = ifelse(WAR_talent > max_WAR_talent+1, max_WAR_talent, WAR_talent))
  })))
  bar %>% mutate(ystar = ystar)
}

#############

#############

# maps their latent talents to a common lab environment and get the era-adjusted statistics
# using non-parametric distribution measuring the components
# talent follows folded normal distribution
career_talent_nonpara_folded <- function(dataset, component_name, snippet){
  Rev_Aptitude_nonpara <- function(x, ytilde, alpha = 1.16, npop){
    
    # transforms ordered Pareto values corresponding to
    # the general population to percentiles from order stats
    # (in increasing order)
    n <- length(x)
    if(length(npop) == 1) npop <- rep(npop, n)
    u = unlist(lapply(1:n, function(j){
      pbeta(pfoldnorm(x[j], mean = 0, sd = 1), j + npop[j]-n, n + 1 - j)
    }))
    
    
    ## map the quantile to the a predicated sample value
    map_Y <- function(u, ytilde){
      n <- length(ytilde)-1
      seqence <- seq(0, 1, 1/n)
      pos <- findInterval(u, seqence)
      out <- (n*u -pos + 1) * (ytilde[(pos+1)] - ytilde[pos]) + ytilde[pos]
      return(out)
    }
    
    ## map the vector of quantiles to the predicated sample values
    n <- length(u)
    a <- qbeta(u, shape1 = 1:n, shape2 = n:1)
    out <- sapply(1:n, function(x) map_Y(a[x], ytilde = ytilde))
    out
    
  }
  
  talent_all <- sort(c(talent_new, snippet$WAR_talent))
  index_snippet <- max(which(talent_all == snippet$WAR_talent))
  
  snippet %>% 
    mutate(adj_comp = Rev_Aptitude_nonpara(talent_all, 
                                           ytilde = ytilde, 
                                           npop = pop_new)[index_snippet])
  
}

#############

#############

# computes their latent talents using non-parametric distribution measuring the components
# the talent follows normal distribution
talent_computing_nonpara_normal <- function(dataset, component_name, year, ystar){
  
  Ftilde <- function(y, t, ystar, component_name){
    y <- sort(y)
    n <- length(y)
    ytilde <- rep(0, n + 1)
    
    if (component_name == 'bWAR' | component_name == 'fWAR' | component_name == 'ERA') {
      ytilde[1] <- y[1] - (y[2] - y[1])
    }
    if (component_name == 'HR'| component_name == 'BB') {
      # since the minimal HR is greater or equal to 0.
      ytilde[1] <- 0
    }
    if (component_name == 'bWAR_p' | component_name == 'fWAR_p') {
      ytilde[1] <- y[1] - (y[2] - y[1])/10
    }
    if (component_name == 'SO' | component_name == 'AVG') {
      ytilde[1] <- ifelse(y[1] - (y[2] - y[1]) < 0, 0, y[1] - (y[2] - y[1]) )
    }
    
    ytilde[n+1] <- y[n] + ystar
    ytilde[2:n] <- unlist(lapply(2:n, function(j){
      (y[j]+y[j-1])/2
    }))
    
    if (t >= ytilde[n+1]) {
      1 - 0.1^7
    } else if (t <= ytilde[1]) {
      0
    } else {
      j <- length(which(ytilde < t))
      (j - 1) / n + (t - ytilde[j]) / (n*(ytilde[j+1] - ytilde[j]))
    }
    
  }
  
  Aptitude_nonpara <- function(p, npop){
    
    # converts order stats to their percentiles
    order_pbino <- function(p = 0, k = 1, n = 1e4){
      pbinom(k - 1, prob = p, size = n, lower.tail = FALSE)
    }
    
    # converts a vector of order stats 
    # to their percentiles. This vector should be the entire 
    # sample sorted in increasing order
    p <- sort(p) # just in case
    n <- length(p)
    u = unlist(lapply(1:n, function(j){
      #order_pbino(p[j], k = 251-n+j, n = 251)
      order_pbino(p[j], k = j, n = n)
    }))
    
    # transforms percentiles from order stats (in increasing order)
    # to Pareto values corresponding to the general population 
    # of a greater than or equal to size
    # default alpha is that of the Pareto principle 80-20
    n <- length(u)
    if(length(npop) == 1) npop <- rep(npop, n)
    unlist(lapply(1:n, function(j){
      qnorm(qbeta(u[j], j + npop[j] -n , n + 1 - j),mean = 0, sd = 1)
    }))
  }
  
  foo <- dataset %>% filter(yearID == year) %>% 
    arrange(comp) 
  bar <- foo %>% filter(full_time == 'Y')
  full_comp <- bar$comp
  ## batter WAR talent
  bar <- bar %>% 
    mutate(WAR_talent = 
             Aptitude_nonpara(p = unlist(lapply(comp, function(xx) 
               Ftilde(y = comp, t = xx, ystar = ystar, component_name = component_name))), npop = pops))
  
  max_WAR_talent <- max(bar$WAR_talent) - 1
  range <- which(!(foo$playerID %in% bar$playerID))
  
  ## using the distribution from full time players
  bar <- rbind(bar, do.call(rbind, lapply(range, function(j){
    rbind(bar %>% dplyr::select(-WAR_talent), foo[j, ]) %>% arrange(comp) %>%
      mutate(WAR_talent = Aptitude_nonpara(p = unlist(lapply(comp, function(xx) 
        Ftilde(y = full_comp, t = xx, ystar = ystar, component_name = component_name))), npop = pops)) %>%
      filter(full_time == 'N') %>% 
      mutate(WAR_talent = ifelse(WAR_talent > max_WAR_talent+1, max_WAR_talent, WAR_talent))
  })))
  bar %>% mutate(ystar = ystar)
}

#############

#############

# maps their latent talents to a common lab environment and get the era-adjusted statistics
# using non-parametric distribution measuring the components
# talent follows normal distribution
career_talent_nonpara_normal <- function(dataset, component_name, snippet){
  Rev_Aptitude_nonpara <- function(x, ytilde, npop){
    
    # transforms ordered Pareto values corresponding to
    # the general population to percentiles from order stats
    # (in increasing order)
    n <- length(x)
    if(length(npop) == 1) npop <- rep(npop, n)
    u = unlist(lapply(1:n, function(j){
      pbeta(pnorm(x[j], mean = 0, sd = 1), j + npop[j]-n, n + 1 - j)
    }))
    
    
    ## map the quantile to the a predicated sample value
    map_Y <- function(u, ytilde){
      n <- length(ytilde)-1
      seqence <- seq(0, 1, 1/n)
      pos <- findInterval(u, seqence)
      out <- (n*u -pos + 1) * (ytilde[(pos+1)] - ytilde[pos]) + ytilde[pos]
      return(out)
    }
    
    ## map the vector of quantiles to the predicated sample values
    n <- length(u)
    a <- qbeta(u, shape1 = 1:n, shape2 = n:1)
    out <- sapply(1:n, function(x) map_Y(a[x], ytilde = ytilde))
    out
    
  }
  
  talent_all <- sort(c(talent_new, snippet$WAR_talent))
  index_snippet <- max(which(talent_all == snippet$WAR_talent))
  
  snippet %>% 
    mutate(adj_comp = Rev_Aptitude_nonpara(talent_all, 
                                           ytilde = ytilde, 
                                           npop = pop_new)[index_snippet])
  
}

#############