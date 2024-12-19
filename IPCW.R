library(timereg)

est <- rep(0,1000)
tet <- rep(0,1000)
aipcw <- rep(0,1000)

set.seed(12345)
n <- 100

# need to specify the increments, which will be based on the number of digits in TT and CC
for (i in 1:1000){
X <- pmax(rnorm(n, mean = 0.6, sd = 0.45),0)
U <- pmax(rnorm(n, mean = 0.6, sd = 0.45),0)
Z <- rnorm(n, mean = 1.4+0.3*X-0.9*U, sd = 0.25)
W <- rnorm(n, mean = 0.6-0.2*X+0.5*U, sd = 0.25)



TT <- rexp(n, 0.25+0.3*X+0.6*U)

CC <- rexp(n, 0.1+0.25*X+1*U)
CC <- pmin(CC,3)



T_tilde <- pmin(TT,CC)
delta <- rep(0, n)
delta_c <- rep(0,n)
delta[which(T_tilde == TT)] <- 1
delta_c[which(T_tilde == CC)] <- 1

delta <- delta[order(T_tilde)]
W <- W[order(T_tilde)]
Z <- Z[order(T_tilde)]
X <- X[order(T_tilde)]
U <- U[order(T_tilde)]
TT <- TT[order(T_tilde)]
CC <- CC[order(T_tilde)]
T_tilde <- T_tilde[order(T_tilde)]


### fitting model

m1 <- aalen(Surv(T_tilde,delta_c) ~ U+X)
m2 <- aalen(Surv(T_tilde,delta) ~ U+X)

A <- function(t){
  
  x <- m1$cum[,1]
  xx <- x[which( x-t  < 0 )]
  return(match(x[which(   abs(xx - t) == min(abs(xx - t)))],x))
  
}



T_T <- T_tilde[which(delta==1)]
T_C <- T_tilde[which(delta==0)]
WW <- W[which(delta==1)]
XX <- X[which(delta==1)]
ZZ <- Z[which(delta==1)]
UU <- U[which(delta==1)]




tet[i] <- mean(TT>0.5)


up <- 0
down <- 0
for (j in 1:length(T_T)) {
  up <-  up+ ((1*(T_T[j]>0.5))/(exp(m1$cum[A(T_T[j]),2]+m1$cum[A(T_T[j]),3]*UU[j]+m1$cum[A(T_T[j]),4]*XX[j])))
  down <- down + (1/(exp(m1$cum[A(T_T[j]),2]+m1$cum[A(T_T[j]),3]*UU[j]+m1$cum[A(T_T[j]),4]*XX[j])))
}

est[i] <- 1 - up/down




print(i)



}


