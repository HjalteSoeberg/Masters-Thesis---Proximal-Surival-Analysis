library(MASS)
library(survival)

set.seed(1234)
n <- 100

est <- rep(0,1000) #the PCE
truth <- rep(0,1000) #true p(T>0.5)
est_PDRE <- rep(0,1000) # PDRE
est_IPCW <- rep(0,1000) # FIE
est_DRE <- rep(0,1000) # DRE



for (g in 1:1000){
  # need to specify the increments, which will be based on the number of digits in TT and CC
  
  X <- pmax(rnorm(n, mean = 0.6, sd = 0.45),0)
  U <- pmax(rnorm(n, mean = 0.6, sd = 0.45),0)
  Z <- rnorm(n, mean = 1.4+0.3*X-0.9*U, sd = 0.25)
  W <- rnorm(n, mean = 0.6-0.2*X+0.5*U, sd = 0.25)
  
  
  #setup 1
  TT <- rexp(n, 0.25+0.3*X+0.6*U)
  
  CC <- rexp(n, 0.1+0.25*X+1*U)
  CC <- pmin(CC,3)
  
  #setup 2
  # TT <- rep(NA, n)
  # for (i in 1:n){
  #   u <- runif(1, min = 0, max = 1)
  #   TT[i] <- 1*(-log(u)*exp(-(0.3*X[i]+0.6*U[i])))
  # }
  # 
  # CC <- rexp(n, 0.1+0.25*X+1*U)
  # CC <- pmin(CC,3)
  
  
  
  #setup3
  # TT <- rexp(n, 0.25+0.3*X+0.6*U)
  # 
  # CC <- rep(NA, n)
  # for (i in 1:n){
  #   u <- runif(1, min = 0, max = 1)
  #   CC[i] <- 3*(-log(u)*exp(-(0.25*X[i]+1*U[i]))) #S(t)=1-F(t)
  #   if (CC[i] >=3){
  #     CC[i] <- 3
  #   }
  # }
  
  
  
  T_tilde <- pmin(TT,CC)
  delta <- rep(0, n)
  delta[which(T_tilde == TT)] <- 1
  
  delta <- delta[order(T_tilde)]
  W <- W[order(T_tilde)]
  Z <- Z[order(T_tilde)]
  X <- X[order(T_tilde)]
  U <- U[order(T_tilde)]
  TT <- TT[order(T_tilde)]
  CC <- CC[order(T_tilde)]
  T_tilde <- T_tilde[order(T_tilde)]
  
  
  
  
  weights <- rep(1,n)
  
  T_C <- T_tilde[which(delta == 0)]
  T_C <- T_C[order(T_C)]
  
  
  A_0 <- rep(0, length(T_C)+1) #A_hat only jumps when dN_C=1, i.e. at all T_tilde with delta == 0
  A_Z <- rep(0, length(T_C)+1) #The +1 stems from the fact that the first entry has to be 0 as per the initial condition
  A_X <- rep(0, length(T_C)+1) #We only need to encode the jumps in dN_C as all teh entries in A_hat are the same until we reach a jump
  
  M_A <- function(t){
    t_A <- match(t,T_C)
    WW <- W[which(T_tilde >= t)]
    XX <- X[which(T_tilde >= t)]
    ZZ <- Z[which(T_tilde >= t)]
    weightss <- weights[which(T_tilde >= t)]
    
    k <- list(weightss[1]*exp(A_0[t_A]+A_Z[t_A]*ZZ[1]+A_X[t_A]*XX[1])*matrix(c(1,WW[1],XX[1]), nrow = 3)%*%matrix(c(1,ZZ[1],XX[1]), nrow = 1))
    for (i in 2:length(WW)){
      k <- append(k,list(weightss[i]*exp(A_0[t_A]+A_Z[t_A]*ZZ[i]+A_X[t_A]*XX[i])*matrix(c(1,WW[i],XX[i]), nrow = 3)%*%matrix(c(1,ZZ[i],XX[i]), nrow = 1)))
    }
    K <- Reduce('+', k)
    
    if (length(WW) == 1){
      k <- weightss[1]*exp(A_0[t_A]+A_Z[t_A]*ZZ[1]+A_X[t_A]*XX[1])*matrix(c(1,WW[1],XX[1]), nrow = 3)%*%matrix(c(1,ZZ[1],XX[1]), nrow = 1)
      K <- k
    }
    if (length(WW) < 1){
      K <- matrix(rep(0,9),nrow = 3)
    }
    return(K/n)
  }
  
  
  A_sum <- function(t){
    t_A <- match(t,T_C)
    WW <- W[which(T_tilde == t & delta == 0)]
    XX <- X[which(T_tilde == t & delta == 0)]
    ZZ <- Z[which(T_tilde == t & delta == 0)]
    weightss <- weights[which(T_tilde == t & delta == 0)]
    K <- weightss[1]*exp(A_0[t_A]+A_Z[t_A]*ZZ[1]+A_X[t_A]*XX[1])*matrix(c(1,WW[1],XX[1]), nrow = 3)
    return(K)
  }
  
  A_hat <- function(t){
    tt <- T_C[which(T_C <= t)] #all the T_C's that satisfy dN_C ==1 in the interval [0,t]
    
    k <- list(ginv(M_A(tt[1]))%*%A_sum(tt[1]))
    for (I in 1:length(tt)){
      k <- append(k,list(ginv(M_A(tt[I]))%*%A_sum(tt[I])))
    }
    k <- Reduce('+', k)
    if (length(tt)==1){
      k <- (ginv(M_A(tt[1]))%*%A_sum(tt[1]))
    }
    return(k/n)
  }
  
  
  #estimating all the A_hats
  for (l in 1:length(T_C)){
    a_hats <- A_hat(T_C[l])
    A_0[l+1] <- a_hats[1]
    A_Z[l+1] <- a_hats[2]
    A_X[l+1] <- a_hats[3]
  }
  
  
  #here t is the t form the indicator function 1(T_tilde_i > t)
  PCE <- function(t){
    T_PCE <- T_tilde[which(delta == 1)]
    W_PCE <- W[which(delta == 1)]
    Z_PCE <- Z[which(delta == 1)]
    X_PCE <- X[which(delta == 1)]
    weights_PCE <- weights[which(delta == 1)]
    T_PCE_up <- T_PCE[which(T_PCE > t)]
    W_PCE_up <- W_PCE[which(T_PCE > t)]
    Z_PCE_up <- Z_PCE[which(T_PCE > t)]
    X_PCE_up <- X_PCE[which(T_PCE > t)] 
    weights_PCE_UP <- weights_PCE[which(T_PCE > t)] 
    
    pce_up <- 0
    for (i in 1:length(T_PCE_up)){
      t_A <- T_C[which(T_C-T_PCE_up[i] < 0)]
      t_A <- t_A[which.min(abs(t_A - T_PCE_up[i]))]
      t_A <- which(T_C==t_A) + 1
      if (length(t_A) < 1){
        t_A <- 1
      }
      pce_up <- pce_up + exp(A_0[t_A]+A_Z[t_A]*Z_PCE_up[i]+A_X[t_A]*X_PCE_up[i])
    }
    
    pce_down <- 0
    for (i in 1:length(T_PCE)){
      t_A <- T_C[which(T_C-T_PCE[i] < 0)]
      t_A <- t_A[which.min(abs(t_A - T_PCE[i]))]
      t_A <- which(T_C==t_A) + 1
      if (length(t_A) < 1){
        t_A <- 1
      }
      pce_down <- pce_down + exp(A_0[t_A]+A_Z[t_A]*Z_PCE[i]+A_X[t_A]*X_PCE[i])
    }
    
    
    pce <- pce_up/pce_down
    return(pce)
  }
  
  
  est[g] <- PCE(0.5) #estimator
  
  
  # finding all T_tilde where delta == 0 and T_tilde <= 0.5 as we are interrested  i P(T>0.5)
  T_T <- T_tilde[which(delta == 1 & T_tilde <= 0.5)]
  T_T <- T_T[order(T_T)]
  
  
  B_0 <- rep(0, length(T_T)+1) #B_hat only jumps when dN_T=1, i.e. at all T_tilde with delta == 1
  B_W <- rep(0, length(T_T)+1) #The +1 stems from the fact that the last entry has to be 0 as per the initial condition
  B_X <- rep(0, length(T_T)+1) #We only need to encode the jumps in dN_T as all the entries in B_hat are the same until we reach a jump. Note that the order from A_hat and B_hat are reversed.
  #in A_hat the values from jump 1 to jump 2 was the same as the value in jump 1, but in B_hat the value would be the same as jump 2, i.e. if there are jumps at t=1 and t=2 then b_hat(1.5)=b_hat(2)
  
  
  
  M_B <- function(t){
    t_B <- match(t,T_T) + 1
    if (length(t_B) < 1){
      t_B <- length(T_T)+1
    }
    WW <- W[which(T_tilde >= t)]
    XX <- X[which(T_tilde >= t)]
    ZZ <- Z[which(T_tilde >= t)]
    weightss <- weights[which(T_tilde >= t)]
    
    k <- list( weightss[1]*exp(B_0[t_B]+B_W[t_B]*WW[1]+B_X[t_B]*XX[1])*(matrix(c(1,ZZ[1],XX[1]), nrow = 3)%*%matrix(c(1,WW[1],XX[1]), nrow = 1)))
    for (i in 2:length(WW)){
      k <- append(k,list( weightss[i]*exp(B_0[t_B]+B_W[t_B]*WW[i]+B_X[t_B]*XX[i])*(matrix(c(1,ZZ[i],XX[i]), nrow = 3)%*%matrix(c(1,WW[i],XX[i]), nrow = 1))))
    }
    K <- Reduce('+', k)
    
    if (length(WW) == 1){
      k <-  weightss[1]*exp(B_0[t_B]+B_W[t_B]*WW[1]+B_X[t_B]*XX[1])*(matrix(c(1,ZZ[1],XX[1]), nrow = 3)%*%matrix(c(1,WW[1],XX[1]), nrow = 1))
      K <- k
    }
    if (length(WW) < 1){
      K <- matrix(rep(0,9),nrow = 3)
    }
    return(K/n)
  }
  
  
  B_sum <- function(t){
    t_B <- match(t,T_T) + 1
    if (length(t_B) < 1){
      t_B <- length(T_T)+1
    }
    WW <- W[which(T_tilde == t & delta == 1)]
    XX <- X[which(T_tilde == t & delta == 1)]
    ZZ <- Z[which(T_tilde == t & delta == 1)]
    weightss <- weights[which(T_tilde == t & delta == 1)]
    K <- weightss[1]*exp(B_0[t_B]+B_W[t_B]*WW[1]+B_X[t_B]*XX[1])*matrix(c(1,ZZ[1],XX[1]), nrow = 3)
    if (length(WW) == 0){
      K <- matrix(c(0,0,0),nrow = 3)
    }
    return(K)
  }
  
  
  B_hat <- function(t){
    tt <- T_T[which(T_T <= t)] #all the T_T's that satisfy dN_T ==1 in the interval [0,t]
    k <- list(ginv(M_B(tt[1]))%*%B_sum(tt[1]))
    if (length(tt) > 1){
      for (I in 2:length(tt)){
        k <- append(k,list(ginv(M_B(tt[I]))%*%B_sum(tt[I])))
      }}
    k <- Reduce('+', k)
    k <- k
    return(k/n)
  }
  
  
  #estimating all the B_hats
  
  for (j in 1:6){
    for (l in 1:length(T_T)){
      L <- length(T_T)-l+1
      
      B_hats <- B_hat(T_T[L])
      B_0[L] <- B_hats[1]
      B_W[L] <- B_hats[2]
      B_X[L] <- B_hats[3]
    }
  }
  
  
  
  
  PDRE <- function(t){
    T_PCE <- T_tilde[which(delta == 1)]
    W_PCE <- W[which(delta == 1)]
    Z_PCE <- Z[which(delta == 1)]
    X_PCE <- X[which(delta == 1)]
    weights_PCE <- weights[which(delta == 1)]
    
    T_PCE_up <- T_PCE[which(T_PCE > t)]
    W_PCE_up <- W_PCE[which(T_PCE > t)]
    Z_PCE_up <- Z_PCE[which(T_PCE > t)]
    X_PCE_up <- X_PCE[which(T_PCE > t)]
    weights_PCE_up <- weights_PCE[which(T_PCE > t)]
    
    W_C <- W[which(delta == 0)]
    Z_C <- Z[which(delta == 0)]
    X_C <- X[which(delta == 0)]
    weights_C <- weights[which(delta == 0)]
    
    # weights_PCE <- rep(1,length(weights_PCE))
    # weights_PCE_up <- rep(1,length(weights_PCE_up))
    # weights_C <- rep(1,length(weights_C))
    
    pce_up <- 0
    for (i in 1:length(T_PCE_up)){
      t_A <- T_C[which(T_C-T_PCE_up[i] < 0)]
      t_A <- t_A[which.min(abs(t_A - T_PCE_up[i]))]
      t_A <- which(T_C==t_A) + 1
      if (length(t_A) < 1){
        t_A <- 1
      }
      pce_up <- pce_up + weights_PCE_up[i]*exp(A_0[t_A]+A_Z[t_A]*Z_PCE_up[i]+A_X[t_A]*X_PCE_up[i])
    }
    
    pce_down <- 0
    for (i in 1:length(T_PCE)){
      t_A <- T_C[which(T_C-T_PCE[i] < 0)]
      t_A <- t_A[which.min(abs(t_A - T_PCE[i]))]
      t_A <- which(T_C==t_A) + 1
      if (length(t_A) < 1){
        t_A <- 1
      }
      pce_down <- pce_down + weights_PCE[i]*exp(A_0[t_A]+A_Z[t_A]*Z_PCE[i]+A_X[t_A]*X_PCE[i])
    }
    
    int_sum_up <- 0
    for (i in 1:length(T_tilde)) {
      int_sum <- 0
      T_CC <- T_C[which(T_C <= T_tilde[i])]
      
      if (length(T_CC) > 0){
        for (j in 1:length(T_CC)) {
          t_A <- T_C[which(T_C-T_CC[j] < 0)]
          t_A <- t_A[which.min(abs(t_A - T_CC[j]))]
          t_A <- which(T_C==t_A) + 1
          if (length(t_A) < 1){
            t_A <- 1
          }
          t_B <- T_T[which(T_T-T_CC[j] > 0)]
          t_B <- t_B[which.min(abs(t_B - T_CC[j]))]
          t_B <- which(T_T==t_B)
          if (length(t_B) < 1){
            t_B <- length(T_T)+1
          }
          ifelse(t_A > 1,
                 int_sum <- int_sum + weights[i]*(exp(A_0[t_A]+A_Z[t_A]*Z[i]+A_X[t_A]*X[i])-exp(A_0[t_A-1]+A_Z[t_A-1]*Z[i]+A_X[t_A-1]*X[i]))*exp((B_0[t_B]+B_W[t_B]*W[i]+B_X[t_B]*X[i])),
                 int_sum <- int_sum + 0
          )
        }
      }
      
      int_sum_up <- int_sum_up + int_sum
    }
    
    
    
    
    
    int_sum_down <- 0
    for (i in 1:length(T_tilde)) {
      t_A <- T_C[which(T_C-T_tilde[i] < 0)]
      t_A <- t_A[which.min(abs(t_A - T_tilde[i]))]
      t_A <- which(T_C==t_A) + 1
      if (length(t_A) < 1){
        t_A <- 1
      }
      int_sum_down <- int_sum_down + weights[i]*exp(A_0[t_A]+A_Z[t_A]*Z[i]+A_X[t_A]*X[i])-1
    }
    
    
    
    
    QdN_C_sum <- 0
    HQdN_C_sum <- 0
    for (i in 1:length(T_C)){
      t_B <- T_T[which(T_T-T_C[i] > 0)]
      t_B <- t_B[which.min(abs(t_B - T_C[i]))]
      t_B <- which(T_T==t_B)
      if (length(t_B) < 1){
        t_B <- length(T_T)+1
      }
      t_A <- T_C[which(T_C-T_C[i] < 0)]
      t_A <- t_A[which.min(abs(t_A - T_C[i]))]
      t_A <- which(T_C==t_A) + 1
      if (length(t_A) < 1){
        t_A <- 1
      }
      QdN_C_sum <- QdN_C_sum + exp(A_0[t_A]+A_Z[t_A]*Z_C[i]+A_X[t_A]*X_C[i])
      ifelse(T_C[i] <= t,
             HQdN_C_sum <- HQdN_C_sum + weights_C[i]*exp((B_0[t_B]+B_W[t_B]*W_C[i]+B_X[t_B]*X_C[i]))*exp(A_0[t_A]+A_Z[t_A]*Z_C[i]+A_X[t_A]*X_C[i]),
             HQdN_C_sum <- HQdN_C_sum + weights_C[i]*exp(A_0[t_A]+A_Z[t_A]*Z_C[i]+A_X[t_A]*X_C[i])
      )
    }
    
    
    pce <- (pce_up-int_sum_up+HQdN_C_sum)/(pce_down-int_sum_down+QdN_C_sum)
    #pce1 <- pce_up/pce_down
    return(pce)
    
    
    
    
    
    
    
    
  }
  
  
  
  est_PDRE[g] <-  PDRE(0.5)
  
  
  ##DRE
  
  {
    T_C <- T_tilde[which(delta == 0)]
    T_C <- T_C[order(T_C)]
    
    
    A_0 <- rep(0, length(T_C)+1) #A_hat only jumps when dN_C=1, i.e. at all T_tilde with delta == 0
    A_Z <- rep(0, length(T_C)+1) #The +1 stems from the fact that the first entry has to be 0 as per the initial condition
    A_X <- rep(0, length(T_C)+1) #We only need to encode the jumps in dN_C as all teh entries in A_hat are the same until we reach a jump
    A_W <- rep(0, length(T_C)+1)
    
    M_A <- function(t){
      t_A <- T_C[which(T_C-t < 0)]
      t_A <- t_A[which.min(abs(t_A - t))]
      t_A <- which(T_C==t_A)
      ifelse(length(t_A) < 1,
             t_A <- 1,
             t_A <- t_A+1)
      WW <- W[which(T_tilde >= t)]
      XX <- X[which(T_tilde >= t)]
      ZZ <- Z[which(T_tilde >= t)]
      weightss <- weights[which(T_tilde >= t)]
      
      k <- list(weightss[1]*exp(A_0[t_A]+A_Z[t_A]*ZZ[1]+A_X[t_A]*XX[1]+A_W[t_A]*WW[1])*matrix(c(1,WW[1],ZZ[1],XX[1]), nrow = 4)%*%matrix(c(1,WW[1],ZZ[1],XX[1]), nrow = 1))
      for (i in 2:length(WW)){
        k <- append(k,list(weightss[i]*exp(A_0[t_A]+A_Z[t_A]*ZZ[i]+A_X[t_A]*XX[i]+A_W[t_A]*WW[i])*matrix(c(1,WW[i],ZZ[i],XX[i]), nrow = 4)%*%matrix(c(1,WW[i],ZZ[i],XX[i]), nrow = 1)))
      }
      K <- Reduce('+', k)
      
      if (length(WW) == 1){
        k <- exp(weightss[1]*A_0[t_A]+A_Z[t_A]*ZZ[1]+A_X[t_A]*XX[1]+A_W[t_A]*WW[1])*matrix(c(1,WW[1],ZZ[1],XX[1]), nrow = 4)%*%matrix(c(1,WW[1],ZZ[1],XX[1]), nrow = 1)
        K <- k
      }
      if (length(WW) < 1){
        K <- matrix(rep(0,16),nrow = 4)
      }
      return(K/n)
    }
    
    
    A_sum <- function(t){
      t_A <- T_C[which(T_C-t < 0)]
      t_A <- t_A[which.min(abs(t_A - t))]
      t_A <- which(T_C==t_A)
      ifelse(length(t_A) < 1,
             t_A <- 1,
             t_A <- t_A+1)
      WW <- W[which(T_tilde == t & delta == 0)]
      XX <- X[which(T_tilde == t & delta == 0)]
      ZZ <- Z[which(T_tilde == t & delta == 0)]
      weightss <- weights[which(T_tilde == t & delta == 0)]
      K <- weightss[1]*exp(A_0[t_A]+A_Z[t_A]*ZZ[1]+A_X[t_A]*XX[1]+A_W[t_A]*WW[1])*matrix(c(1,WW[1],ZZ[1],XX[1]), nrow = 4)
      return(K)
    }
    
    A_hat <- function(t){
      tt <- T_C[which(T_C <= t)] #all the T_C's that satisfy dN_C ==1 in the interval [0,t]
      
      k <- list(ginv(M_A(tt[1]))%*%A_sum(tt[1]))
      for (i in 1:length(tt)){
        k <- append(k,list(ginv(M_A(tt[i]))%*%A_sum(tt[i])))
      }
      k <- Reduce('+', k)
      if (length(tt)==1){
        k <- (ginv(M_A(tt[1]))%*%A_sum(tt[1]))
      }
      return(k/n)
    }
    
    
    #estimating all the A_hats
    for (i in 1:length(T_C)){
      a_hats <- A_hat(T_C[i])
      A_0[i+1] <- a_hats[1]
      A_W[i+1] <- a_hats[2]
      A_Z[i+1] <- a_hats[3]
      A_X[i+1] <- a_hats[4]
    }
    
    
    
    # finding all T_tilde where delta == 0 and T_tilde <= 0.5 as we are interrested  i P(T>0.5)
    T_T <- T_tilde[which(delta == 1 & T_tilde <= 0.5)]
    T_T <- T_T[order(T_T)]
    
    
    B_0 <- rep(0, length(T_T)+1) #B_hat only jumps when dN_T=1, i.e. at all T_tilde with delta == 1
    B_W <- rep(0, length(T_T)+1) #The +1 stems from the fact that the last entry has to be 0 as per the initial condition
    B_X <- rep(0, length(T_T)+1) #We only need to encode the jumps in dN_T as all the entries in B_hat are the same until we reach a jump. Note that the order from A_hat and B_hat are reversed.
    B_Z <- rep(0, length(T_T)+1)                           #in A_hat the values from jump 1 to jump 2 was the same as the value in jump 1, but in B_hat the value would be the same as jump 2, i.e. if there are jumps at t=1 and t=2 then b_hat(1.5)=b_hat(2)
    
    
    
    M_B <- function(t){
      t_B <- T_T[which(T_T-t > 0)]
      t_B <- t_B[which.min(abs(t_B - t))]
      t_B <- which(T_T==t_B)
      if (length(t_B) < 1){
        t_B <- length(T_T)+1
      }
      WW <- W[which(T_tilde >= t)]
      XX <- X[which(T_tilde >= t)]
      ZZ <- Z[which(T_tilde >= t)]
      weightss <- weights[which(T_tilde >= t)]
      
      k <- list(weightss[1]*exp(B_0[t_B]+B_W[t_B]*WW[1]+B_X[t_B]*XX[1]+B_Z[t_B]*ZZ[1])*(matrix(c(1,WW[1],ZZ[1],XX[1]), nrow = 4)%*%matrix(c(1,WW[1],ZZ[1],XX[1]), nrow = 1)))
      for (i in 2:length(WW)){
        k <- append(k,list(weightss[i]*exp(B_0[t_B]+B_W[t_B]*WW[i]+B_X[t_B]*XX[i]+B_Z[t_B]*ZZ[i])*(matrix(c(1,WW[i],ZZ[i],XX[i]), nrow = 4)%*%matrix(c(1,WW[i],ZZ[i],XX[i]), nrow = 1))))
      }
      K <- Reduce('+', k)
      
      if (length(WW) == 1){
        k <- weightss[1]*exp(B_0[t_B]+B_W[t_B]*WW[1]+B_X[t_B]*XX[1]+B_Z[t_B]*ZZ[1])*(matrix(c(1,WW[1],ZZ[1],XX[1]), nrow = 4)%*%matrix(c(1,WW[1],ZZ[1],XX[1]), nrow = 1))
        K <- k
      }
      if (length(WW) < 1){
        K <- matrix(rep(0,16),nrow = 4)
      }
      return(K/n)
    }
    
    
    B_sum <- function(t){
      t_B <- T_T[which(T_T-t > 0)]
      t_B <- t_B[which.min(abs(t_B - t))]
      t_B <- which(T_T==t_B)
      if (length(t_B) < 1){
        t_B <- length(T_T)+1
      }
      WW <- W[which(T_tilde == t & delta == 1)]
      XX <- X[which(T_tilde == t & delta == 1)]
      ZZ <- Z[which(T_tilde == t & delta == 1)]
      weightss <- weights[which(T_tilde == t & delta == 1)]
      K <- weightss[1]*exp(B_0[t_B]+B_W[t_B]*WW[1]+B_X[t_B]*XX[1]+B_Z[t_B]*ZZ[1])*matrix(c(1,WW[1],ZZ[1],XX[1]), nrow = 4)
      if (length(WW) == 0){
        K <- matrix(c(0,0,0,0),nrow = 4)
      }
      return(K)
    }
    
    
    B_hat <- function(t){
      tt <- T_T[which(T_T <= t)] #all the T_T's that satisfy dN_T ==1 in the interval [0,t]
      k <- list(ginv(M_B(tt[1]))%*%B_sum(tt[1]))
      if (length(tt) > 1){
        for (i in 2:length(tt)){
          k <- append(k,list(ginv(M_B(tt[i]))%*%B_sum(tt[i])))
        }}
      k <- Reduce('+', k)
      k <- k
      return(k/n)
    }
    
    
    #estimating all the B_hats
    
    for (j in 1:6){
      for (l in 1:length(T_T)){
        L <- length(T_T)-l+1
        
        B_hats <- B_hat(T_T[L])
        B_0[L] <- B_hats[1]
        B_W[L] <- B_hats[2]
        B_Z[L] <- B_hats[3]
        B_X[L] <- B_hats[4]
      }
    }
    
    DRE <- function(t){
      T_PCE <- T_tilde[which(delta == 1)]
      W_PCE <- W[which(delta == 1)]
      Z_PCE <- Z[which(delta == 1)]
      X_PCE <- X[which(delta == 1)]
      T_PCE_up <- T_PCE[which(T_PCE > t)]
      W_PCE_up <- W_PCE[which(T_PCE > t)]
      Z_PCE_up <- Z_PCE[which(T_PCE > t)]
      X_PCE_up <- X_PCE[which(T_PCE > t)]
      W_C <- W[which(delta == 0)]
      Z_C <- Z[which(delta == 0)]
      X_C <- X[which(delta == 0)]
      
      pce_up <- 0
      for (i in 1:length(T_PCE_up)){
        t_A <- T_C[which(T_C-T_PCE_up[i] < 0)]
        t_A <- t_A[which.min(abs(t_A - T_PCE_up[i]))]
        t_A <- which(T_C==t_A)
        if (length(t_A) < 1){
          t_A <- 1
        }
        pce_up <- pce_up + exp(A_0[t_A]+A_Z[t_A]*Z_PCE_up[i]+A_X[t_A]*X_PCE_up[i]+A_W[t_A]*W_PCE_up[i])
      }
      
      pce_down <- 0
      for (i in 1:length(T_PCE)){
        t_A <- T_C[which(T_C-T_PCE[i] < 0)]
        t_A <- t_A[which.min(abs(t_A - T_PCE[i]))]
        t_A <- which(T_C==t_A)
        if (length(t_A) < 1){
          t_A <- 1
        }
        pce_down <- pce_down + exp(A_0[t_A]+A_Z[t_A]*Z_PCE[i]+A_X[t_A]*X_PCE[i]+A_W[t_A]*W_PCE[i])
      }
      
      int_sum_up <- 0
      for (i in 1:length(T_tilde)) {
        int_sum <- 0
        T_CC <- T_C[which(T_C <= T_tilde[i])]
        
        if (length(T_CC) > 0){
          for (j in 1:length(T_CC)) {
            t_A <- T_C[which(T_C-T_CC[j] < 0)]
            t_A <- t_A[which.min(abs(t_A - T_CC[j]))]
            t_A <- which(T_C==t_A)
            if (length(t_A) < 1){
              t_A <- 1
            }
            t_B <- T_T[which(T_T-T_CC[j] > 0)]
            t_B <- t_B[which.min(abs(t_B - T_CC[j]))]
            t_B <- which(T_T==t_B)
            if (length(t_B) < 1){
              t_B <- length(T_T)+1
            }
            ifelse(t_A > 1,
                   int_sum <- int_sum + (exp(A_0[t_A]+A_Z[t_A]*Z[i]+A_X[t_A]*X[i]+A_W[t_A]*W[i])-exp(A_0[t_A-1]+A_Z[t_A-1]*Z[i]+A_X[t_A-1]*X[i]+A_W[t_A-1]*W[i]))*exp((B_0[t_B]+B_W[t_B]*W[i]+B_X[t_B]*X[i]+B_Z[t_B]*Z[i])),
                   int_sum <- int_sum + 0
            )
          }
        }
        
        int_sum_up <- int_sum_up + int_sum
      }
      
      
      
      
      
      int_sum_down <- 0
      for (i in 1:length(T_tilde)) {
        t_A <- T_C[which(T_C-T_tilde[i] < 0)]
        t_A <- t_A[which.min(abs(t_A - T_tilde[i]))]
        t_A <- which(T_C==t_A)
        if (length(t_A) < 1){
          t_A <- 1
        }
        int_sum_down <- int_sum_down + exp(A_0[t_A]+A_Z[t_A]*Z[i]+A_X[t_A]*X[i]+A_W[t_A]*W[i])-1
      }
      
      
      
      
      QdN_C_sum <- 0
      HQdN_C_sum <- 0
      for (i in 1:length(T_C)){
        t_B <- T_T[which(T_T-T_C[i] > 0)]
        t_B <- t_B[which.min(abs(t_B - T_C[i]))]
        t_B <- which(T_T==t_B)
        if (length(t_B) < 1){
          t_B <- length(T_T)+1
        }
        t_A <- T_C[which(T_C-T_C[i] < 0)]
        t_A <- t_A[which.min(abs(t_A - T_C[i]))]
        t_A <- which(T_C==t_A)
        if (length(t_A) < 1){
          t_A <- 1
        }
        QdN_C_sum <- QdN_C_sum + exp(A_0[t_A]+A_Z[t_A]*Z_C[i]+A_X[t_A]*X_C[i])
        ifelse(T_C[i] <= t,
               HQdN_C_sum <- HQdN_C_sum + exp((B_0[t_B]+B_W[t_B]*W_C[i]+B_X[t_B]*X_C[i]+B_Z[t_B]*Z_C[i]))*exp(A_0[t_A]+A_Z[t_A]*Z_C[i]+A_X[t_A]*X_C[i]+A_W[t_A]*W_C[i]),
               HQdN_C_sum <- HQdN_C_sum + exp(A_0[t_A]+A_Z[t_A]*Z_C[i]+A_X[t_A]*X_C[i]+A_W[t_A]*W_C[i])
        )
      }
      
      
      pce <- (pce_up-int_sum_up+HQdN_C_sum)/(pce_down-int_sum_down+QdN_C_sum)
      #pce1 <- pce_up/pce_down
      return(pce)
    }
    
    
    
    est_DRE[g] <- DRE(0.5)
    
    
    #IPCW
    PCE <- function(t){
      T_PCE <- T_tilde[which(delta == 1)]
      W_PCE <- W[which(delta == 1)]
      Z_PCE <- Z[which(delta == 1)]
      X_PCE <- X[which(delta == 1)]
      weights_PCE <- weights[which(delta == 1)]
      T_PCE_up <- T_PCE[which(T_PCE > t)]
      W_PCE_up <- W_PCE[which(T_PCE > t)]
      Z_PCE_up <- Z_PCE[which(T_PCE > t)]
      X_PCE_up <- X_PCE[which(T_PCE > t)] 
      weights_PCE_UP <- weights_PCE[which(T_PCE > t)] 
      
      pce_up <- 0
      for (i in 1:length(T_PCE_up)){
        t_A <- T_C[which(T_C-T_PCE_up[i] < 0)]
        t_A <- t_A[which.min(abs(t_A - T_PCE_up[i]))]
        t_A <- which(T_C==t_A) + 1
        if (length(t_A) < 1){
          t_A <- 1
        }
        pce_up <- pce_up + exp(A_0[t_A]+A_Z[t_A]*Z_PCE_up[i]+A_X[t_A]*X_PCE_up[i]+A_W[t_A]*W_PCE_up[i])
      }
      
      pce_down <- 0
      for (i in 1:length(T_PCE)){
        t_A <- T_C[which(T_C-T_PCE[i] < 0)]
        t_A <- t_A[which.min(abs(t_A - T_PCE[i]))]
        t_A <- which(T_C==t_A) + 1
        if (length(t_A) < 1){
          t_A <- 1
        }
        pce_down <- pce_down + exp(A_0[t_A]+A_Z[t_A]*Z_PCE[i]+A_X[t_A]*X_PCE[i]+A_W[t_A]*W_PCE[i])
      }
      
      
      pce <- pce_up/pce_down
      return(pce)
    }
    
  }
  
  est_IPCW[g] <- PCE(0.5)
  
  
  
  
  
  truth[g] <- mean(TT>0.5)
  
  print(g)
  
  
}


##add observations to dataframe
#df <- data.frame()


set.seed(1234)
n <- 100

truth <- rep(0,1000) #true p(T>0.5
est_IPCW_truth <- rep(0,1000) # FIE-truth
est_DRE_truth <- rep(0,1000) # DRE-truth
est_KM <- rep(0,1000) # KM


for (g in 1:1000){
  # need to specify the increments, which will be based on the number of digits in TT and CC
  
  X <- pmax(rnorm(n, mean = 0.6, sd = 0.45),0)
  U <- pmax(rnorm(n, mean = 0.6, sd = 0.45),0)
  Z <- rnorm(n, mean = 1.4+0.3*X-0.9*U, sd = 0.25)
  W <- rnorm(n, mean = 0.6-0.2*X+0.5*U, sd = 0.25)
  
  
  
  #setup 1
  TT <- rexp(n, 0.25+0.3*X+0.6*U)
  
  CC <- rexp(n, 0.1+0.25*X+1*U)
  CC <- pmin(CC,3)
  
  #setup 2
  # TT <- rep(NA, n)
  # for (i in 1:n){
  #   u <- runif(1, min = 0, max = 1)
  #   TT[i] <- 1*(-log(u)*exp(-(0.3*X[i]+0.6*U[i])))
  # }
  # 
  # CC <- rexp(n, 0.1+0.25*X+1*U)
  # CC <- pmin(CC,3)
  
  
  
  #setup3
  # TT <- rexp(n, 0.25+0.3*X+0.6*U)
  # 
  # CC <- rep(NA, n)
  # for (i in 1:n){
  #   u <- runif(1, min = 0, max = 1)
  #   CC[i] <- 3*(-log(u)*exp(-(0.25*X[i]+1*U[i]))) #S(t)=1-F(t)
  #   if (CC[i] >=3){
  #     CC[i] <- 3
  #   }
  # }
  
  
  
  T_tilde <- pmin(TT,CC)
  delta <- rep(0, n)
  delta[which(T_tilde == TT)] <- 1
  
  delta <- delta[order(T_tilde)]
  W <- W[order(T_tilde)]
  Z <- Z[order(T_tilde)]
  X <- X[order(T_tilde)]
  U <- U[order(T_tilde)]
  TT <- TT[order(T_tilde)]
  CC <- CC[order(T_tilde)]
  T_tilde <- T_tilde[order(T_tilde)]
  
  
  
  weights <- rep(1,n)
  
  df <- data.frame(T_tilde,delta,X,Z,W)
  km <- with(df, Surv(T_tilde, delta))
  km_fit <- survfit(Surv(T_tilde, delta) ~ 1, data=df)
  k <- summary(km_fit, times = c(0.5))
  est_KM[g] <- k$surv
  
  ## DRE-truth
  
  # finding all T_tilde where delta == 0
  T_C <- T_tilde[which(delta == 0)]
  T_C <- T_C[order(T_C)]
  
  
  A_0 <- rep(0, length(T_C)+1) #A_hat only jumps Uhen dN_C=1, i.e. at all T_tilde Uith delta == 0
  A_U <- rep(0, length(T_C)+1) #The +1 stems from the fact that the first entry has to be 0 as per the initial condition
  A_X <- rep(0, length(T_C)+1) #Ue only need to encode the jumps in dN_C as all teh entries in A_hat are the same until we reach a jump
  
  M_A <- function(t){
    t_A <- T_C[which(T_C-t < 0)]
    t_A <- t_A[which.min(abs(t_A - t))]
    t_A <- which(T_C==t_A)
    ifelse(length(t_A) < 1,
           t_A <- 1,
           t_A <- t_A+1)
    UU <- U[which(T_tilde >= t)]
    XX <- X[which(T_tilde >= t)]
    UU <- U[which(T_tilde >= t)]
    weightss <- weights[which(T_tilde >= t)]
    
    k <- list(weightss[1]*exp(A_0[t_A]+A_U[t_A]*UU[1]+A_X[t_A]*XX[1])*matrix(c(1,UU[1],XX[1]), nrow = 3)%*%matrix(c(1,UU[1],XX[1]), nrow = 1))
    for (i in 2:length(UU)){
      k <- append(k,list(weightss[i]*exp(A_0[t_A]+A_U[t_A]*UU[i]+A_X[t_A]*XX[i])*matrix(c(1,UU[i],XX[i]), nrow = 3)%*%matrix(c(1,UU[i],XX[i]), nrow = 1)))
    }
    K <- Reduce('+', k)
    
    if (length(UU) == 1){
      k <- weightss[1]*exp(A_0[t_A]+A_U[t_A]*UU[1]+A_X[t_A]*XX[1])*matrix(c(1,UU[1],XX[1]), nrow = 3)%*%matrix(c(1,UU[1],XX[1]), nrow = 1)
      K <- k
    }
    if (length(UU) < 1){
      K <- matrix(rep(0,9),nrow = 3)
    }
    return(K/n)
  }
  
  
  A_sum <- function(t){
    t_A <- T_C[which(T_C-t < 0)]
    t_A <- t_A[which.min(abs(t_A - t))]
    t_A <- which(T_C==t_A)
    ifelse(length(t_A) < 1,
           t_A <- 1,
           t_A <- t_A+1)
    UU <- U[which(T_tilde == t & delta == 0)]
    XX <- X[which(T_tilde == t & delta == 0)]
    UU <- U[which(T_tilde == t & delta == 0)]
    weightss <- weights[which(T_tilde == t & delta == 0)]
    K <- weightss[1]*exp(A_0[t_A]+A_U[t_A]*UU[1]+A_X[t_A]*XX[1])*matrix(c(1,UU[1],XX[1]), nrow = 3)
    return(K)
  }
  
  A_hat <- function(t){
    tt <- T_C[which(T_C <= t)] #all the T_C's that satisfy dN_C ==1 in the interval [0,t]
    
    k <- list(ginv(M_A(tt[1]))%*%A_sum(tt[1]))
    for (i in 1:length(tt)){
      k <- append(k,list(ginv(M_A(tt[i]))%*%A_sum(tt[i])))
    }
    k <- Reduce('+', k)
    if (length(tt)==1){
      k <- (ginv(M_A(tt[1]))%*%A_sum(tt[1]))
    }
    return(k/n)
  }
  
  
  #estimating all the A_hats
  for (i in 1:length(T_C)){
    a_hats <- A_hat(T_C[i])
    A_0[i+1] <- a_hats[1]
    A_U[i+1] <- a_hats[2]
    A_X[i+1] <- a_hats[3]
  }
  
  # finding all T_tilde where delta == 0 and T_tilde <= 0.5 as we are interrested  i P(T>0.5)
  T_T <- T_tilde[which(delta == 1 & T_tilde <= 0.5)]
  T_T <- T_T[order(T_T)]
  
  
  B_0 <- rep(0, length(T_T)+1) #B_hat only jumps when dN_T=1, i.e. at all T_tilde with delta == 1
  B_U <- rep(0, length(T_T)+1) #The +1 stems from the fact that the last entry has to be 0 as per the initial condition
  B_X <- rep(0, length(T_T)+1) #Ue only need to encode the jumps in dN_T as all the entries in B_hat are the same until we reach a jump. Note that the order from A_hat and B_hat are reversed.
  #in A_hat the values from jump 1 to jump 2 was the same as the value in jump 1, but in B_hat the value would be the same as jump 2, i.e. if there are jumps at t=1 and t=2 then b_hat(1.5)=b_hat(2)
  
  
  
  M_B <- function(t){
    t_B <- T_T[which(T_T-t > 0)]
    t_B <- t_B[which.min(abs(t_B - t))]
    t_B <- which(T_T==t_B)
    if (length(t_B) < 1){
      t_B <- length(T_T)+1
    }
    UU <- U[which(T_tilde >= t)]
    XX <- X[which(T_tilde >= t)]
    UU <- U[which(T_tilde >= t)]
    weightss <- weights[which(T_tilde >= t)]
    
    k <- list(weightss[1]*exp(B_0[t_B]+B_U[t_B]*UU[1]+B_X[t_B]*XX[1])*(matrix(c(1,UU[1],XX[1]), nrow = 3)%*%matrix(c(1,UU[1],XX[1]), nrow = 1)))
    for (i in 2:length(UU)){
      k <- append(k,list(weightss[i]*exp(B_0[t_B]+B_U[t_B]*UU[i]+B_X[t_B]*XX[i])*(matrix(c(1,UU[i],XX[i]), nrow = 3)%*%matrix(c(1,UU[i],XX[i]), nrow = 1))))
    }
    K <- Reduce('+', k)
    
    if (length(UU) == 1){
      k <- weightss[1]*exp(B_0[t_B]+B_U[t_B]*UU[1]+B_X[t_B]*XX[1])*(matrix(c(1,UU[1],XX[1]), nrow = 3)%*%matrix(c(1,UU[1],XX[1]), nrow = 1))
      K <- k
    }
    if (length(UU) < 1){
      K <- matrix(rep(0,9),nrow = 3)
    }
    return(K/n)
  }
  
  
  B_sum <- function(t){
    t_B <- T_T[which(T_T-t > 0)]
    t_B <- t_B[which.min(abs(t_B - t))]
    t_B <- which(T_T==t_B)
    if (length(t_B) < 1){
      t_B <- length(T_T)+1
    }
    UU <- U[which(T_tilde == t & delta == 1)]
    XX <- X[which(T_tilde == t & delta == 1)]
    UU <- U[which(T_tilde == t & delta == 1)]
    weightss <- weights[which(T_tilde == t & delta == 1)]
    K <- weightss[1]*exp(B_0[t_B]+B_U[t_B]*UU[1]+B_X[t_B]*XX[1])*matrix(c(1,UU[1],XX[1]), nrow = 3)
    if (length(UU) == 0){
      K <- matrix(c(0,0,0),nrow = 3)
    }
    return(K)
  }
  
  
  B_hat <- function(t){
    tt <- T_T[which(T_T <= t)] #all the T_T's that satisfy dN_T ==1 in the interval [0,t]
    k <- list(ginv(M_B(tt[1]))%*%B_sum(tt[1]))
    if (length(tt) > 1){
      for (i in 2:length(tt)){
        k <- append(k,list(ginv(M_B(tt[i]))%*%B_sum(tt[i])))
      }}
    k <- Reduce('+', k)
    k <- k
    return(k/n)
  }
  
  
  #estimating all the B_hats
  
  for (j in 1:6){
    for (l in 1:length(T_T)){
      L <- length(T_T)-l+1
      
      B_hats <- B_hat(T_T[L])
      B_0[L] <- B_hats[1]
      B_U[L] <- B_hats[2]
      B_X[L] <- B_hats[3]
    }
  }
  
  
  PDRE <- function(t){
    T_PCE <- T_tilde[which(delta == 1)]
    U_PCE <- U[which(delta == 1)]
    U_PCE <- U[which(delta == 1)]
    X_PCE <- X[which(delta == 1)]
    T_PCE_up <- T_PCE[which(T_PCE > t)]
    U_PCE_up <- U_PCE[which(T_PCE > t)]
    U_PCE_up <- U_PCE[which(T_PCE > t)]
    X_PCE_up <- X_PCE[which(T_PCE > t)]
    U_C <- U[which(delta == 0)]
    U_C <- U[which(delta == 0)]
    X_C <- X[which(delta == 0)]
    
    pce_up <- 0
    for (i in 1:length(T_PCE_up)){
      t_A <- T_C[which(T_C-T_PCE_up[i] < 0)]
      t_A <- t_A[which.min(abs(t_A - T_PCE_up[i]))]
      t_A <- which(T_C==t_A)
      if (length(t_A) < 1){
        t_A <- 1
      }
      pce_up <- pce_up + exp(A_0[t_A]+A_U[t_A]*U_PCE_up[i]+A_X[t_A]*X_PCE_up[i])
    }
    
    pce_down <- 0
    for (i in 1:length(T_PCE)){
      t_A <- T_C[which(T_C-T_PCE[i] < 0)]
      t_A <- t_A[which.min(abs(t_A - T_PCE[i]))]
      t_A <- which(T_C==t_A)
      if (length(t_A) < 1){
        t_A <- 1
      }
      pce_down <- pce_down + exp(A_0[t_A]+A_U[t_A]*U_PCE[i]+A_X[t_A]*X_PCE[i])
    }
    
    int_sum_up <- 0
    for (i in 1:length(T_tilde)) {
      int_sum <- 0
      T_CC <- T_C[which(T_C <= T_tilde[i])]
      
      if (length(T_CC) > 0){
        for (j in 1:length(T_CC)) {
          t_A <- T_C[which(T_C-T_CC[j] < 0)]
          t_A <- t_A[which.min(abs(t_A - T_CC[j]))]
          t_A <- which(T_C==t_A)
          if (length(t_A) < 1){
            t_A <- 1
          }
          t_B <- T_T[which(T_T-T_CC[j] > 0)]
          t_B <- t_B[which.min(abs(t_B - T_CC[j]))]
          t_B <- which(T_T==t_B)
          if (length(t_B) < 1){
            t_B <- length(T_T)+1
          }
          ifelse(t_A > 1,
                 int_sum <- int_sum + (exp(A_0[t_A]+A_U[t_A]*U[i]+A_X[t_A]*X[i])-exp(A_0[t_A-1]+A_U[t_A-1]*U[i]+A_X[t_A-1]*X[i]))*exp((B_0[t_B]+B_U[t_B]*U[i]+B_X[t_B]*X[i])),
                 int_sum <- int_sum + 0
          )
        }
      }
      
      int_sum_up <- int_sum_up + int_sum
    }
    
    
    
    
    
    int_sum_down <- 0
    for (i in 1:length(T_tilde)) {
      t_A <- T_C[which(T_C-T_tilde[i] < 0)]
      t_A <- t_A[which.min(abs(t_A - T_tilde[i]))]
      t_A <- which(T_C==t_A)
      if (length(t_A) < 1){
        t_A <- 1
      }
      int_sum_down <- int_sum_down + exp(A_0[t_A]+A_U[t_A]*U[i]+A_X[t_A]*X[i])-1
    }
    
    
    
    
    QdN_C_sum <- 0
    HQdN_C_sum <- 0
    for (i in 1:length(T_C)){
      t_B <- T_T[which(T_T-T_C[i] > 0)]
      t_B <- t_B[which.min(abs(t_B - T_C[i]))]
      t_B <- which(T_T==t_B)
      if (length(t_B) < 1){
        t_B <- length(T_T)+1
      }
      t_A <- T_C[which(T_C-T_C[i] < 0)]
      t_A <- t_A[which.min(abs(t_A - T_C[i]))]
      t_A <- which(T_C==t_A)
      if (length(t_A) < 1){
        t_A <- 1
      }
      QdN_C_sum <- QdN_C_sum + exp(A_0[t_A]+A_U[t_A]*U_C[i]+A_X[t_A]*X_C[i])
      ifelse(T_C[i] <= t,
             HQdN_C_sum <- HQdN_C_sum + exp((B_0[t_B]+B_U[t_B]*U_C[i]+B_X[t_B]*X_C[i]))*exp(A_0[t_A]+A_U[t_A]*U_C[i]+A_X[t_A]*X_C[i]),
             HQdN_C_sum <- HQdN_C_sum + exp(A_0[t_A]+A_U[t_A]*U_C[i]+A_X[t_A]*X_C[i])
      )
    }
    
    
    pce <- (pce_up-int_sum_up+HQdN_C_sum)/(pce_down-int_sum_down+QdN_C_sum)
    #pce1 <- pce_up/pce_down
    return(pce)
  }
  
  
  
  est_DRE_truth[g] <- PDRE(0.5)
  
  
  PCE <- function(t){
    T_PCE <- T_tilde[which(delta == 1)]
    W_PCE <- W[which(delta == 1)]
    Z_PCE <- Z[which(delta == 1)]
    X_PCE <- X[which(delta == 1)]
    U_PCE <- U[which(delta == 1)]
    T_PCE_up <- T_PCE[which(T_PCE > t)]
    W_PCE_up <- W_PCE[which(T_PCE > t)]
    Z_PCE_up <- Z_PCE[which(T_PCE > t)]
    X_PCE_up <- X_PCE[which(T_PCE > t)] 
    U_PCE_up <- U_PCE[which(T_PCE > t)] 
    pce_up <- 0
    
    for (i in 1:length(T_PCE_up)){
      t_A <- T_C[which(T_C-T_PCE_up[i] < 0)]
      t_A <- t_A[which.min(abs(t_A - T_PCE_up[i]))]
      t_A <- which(T_C==t_A) + 1
      if (length(t_A) < 1){
        t_A <- 1
      }
      pce_up <- pce_up + exp(A_0[t_A]+A_U[t_A]*U_PCE_up[i]+A_X[t_A]*X_PCE_up[i])
    }
    
    pce_down <- 0
    for (i in 1:length(T_PCE)){
      t_A <- T_C[which(T_C-T_PCE[i] < 0)]
      t_A <- t_A[which.min(abs(t_A - T_PCE[i]))]
      t_A <- which(T_C==t_A) + 1
      if (length(t_A) < 1){
        t_A <- 1
      }
      pce_down <- pce_down + exp(A_0[t_A]+A_U[t_A]*U_PCE[i]+A_X[t_A]*X_PCE[i])
    }
    pce <- pce_up/pce_down
    return(pce)
  }
  
  est_IPCW_truth[g] <- PCE(0.5)
  
  
  
  
  
  
  truth[g] <- mean(TT>0.5)
  
  
  print(g)
  
  
}

##add observations to dataframe
#df <- data.frame()


