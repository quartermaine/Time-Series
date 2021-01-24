## b) -------------------------------------------------
library(astsa)
## script given for the lab ---------------------------
# generate data
set.seed(12345); num=50
w=rnorm(num + 1, 0, 1);
v=rnorm( num, 0, 1)
mu=cumsum(w) # state : mu [ 0 ] , mu [ 1 ] ,... , mu [ 5 0 ]
y = mu[-1] + v # obs : y [ 1 ] ,... , y [ 5 0 ]
Time = 1:num

# filter and smooth ( Ksmooth 0 does both )
ks=Ksmooth0(num, y, A = 1, mu0 = 0, Sigma0 = 1, Phi = 1, cQ = 1, cR = 1)

# moving average smoother for comparison
mysmoother <- function(x, n){filter(as.vector(x),rep(1 / n, n), sides = 2)} # x=vector n=order

# start figure
png(filename="images/plot01.png", width = 1000, height = 300)
par(mfrow = c( 1, 1)); Time = 1:num
# plot(Time, mu[-1], main = 'Predict', ylim = c(-5, 10))
# lines(Time ,y , col = 'green' )
# lines(ks$xp) # one-step-ahead prediction of the state
# lines(ks$xp + 2 * sqrt (ks$Pp), lty = 2, col = 4)
# lines(ks$xp - 2 * sqrt (ks$Pp), lty = 2, col = 4)
plot(Time, mu[-1], main = 'Filter and Smoother', ylim = c(-5, 13))
lines(Time, y, col='green')
lines(ks$xf)
lines(ks$xf + 2 * sqrt (ks$Pf) ,lty = 2, col = 4)
lines(ks $ xf - 2 * sqrt (ks$Pf) , lty = 2 , col = 4 )
lines(Time, mysmoother(y,5), col='red')
legend(1,13, legend=c('Hidden States', 'Prediction', 'Confidence', 'Smoother^5'), col=c('black', 'green', 'blue', 'red'), lty=c(1, 1, 2, 1), cex=0.8)
# plot(Time, mu[-1] , main = 'Smooth', ylim = c (-5 ,10))
# lines(Time, y, col = 'green')
# lines(ks$xs) # state smoothers
# lines(ks$xs + 2 * sqrt(ks$Ps), lty = 2, col = 4)
# lines(ks$xs - 2 * sqrt(ks$Ps), lty = 2, col = 4)
# dev.off()
# mu[1]; ks$x0n; sqrt(ks$P0n) # initial value info



## c) -------------------------------------------------
# start figure
png(filename="images/plot02.png", width = 1000, height = 600)
par(mfrow = c( 2, 1)); Time = 1:num
# plot(Time, mu[-1], main = 'Predict', ylim = c(-5, 10))
# lines(Time ,y , col = 'green' )
# lines(ks$xp)
# lines(ks$xp + 2 * sqrt (ks$Pp), lty = 2, col = 4)
# lines(ks$xp - 2 * sqrt (ks$Pp), lty = 2, col = 4)
ks=Ksmooth0(num, y, A = 1, mu0 = 0, Sigma0 = 1, Phi = 1, cQ = 1, cR = 1)
plot(Time, mu[-1], main = 'Filter Q=1 R=1', ylim = c(-5, 13))
lines(Time ,y , col = 'green')
lines(ks$xf)
lines(ks$xf + 2 * sqrt (ks$Pf) ,lty = 2, col = 4)
lines(ks $ xf - 2 * sqrt (ks$Pf) , lty = 2 , col = 4 )
legend(1,13, legend=c('Hidden States', 'Prediction', 'Confidence'), col=c('black', 'green', 'blue'), lty=c(1, 1, 2), cex=0.8)
ks=Ksmooth0(num, y, A = 1, mu0 = 0, Sigma0 = 1, Phi = 1, cQ = 10, cR = 1)
plot(Time, mu[-1], main = 'Filter Q=10 R= 1', ylim = c(-5, 13))
lines(Time ,y , col = 'green')
lines(ks$xf)
lines(ks$xf + 2 * sqrt (ks$Pf) ,lty = 2, col = 4)
lines(ks $ xf - 2 * sqrt ( ks$Pf ) , lty = 2 , col = 4 )
# plot(Time, mu[-1] , main = 'Smooth', ylim = c (-5 ,10))
# lines(Time, y, col = 'green')
# lines(ks$xs)
# lines(ks$xs + 2 * sqrt(ks$Ps), lty = 2, col = 4)
# lines(ks$xs - 2 * sqrt(ks$Ps), lty = 2, col = 4)
dev.off()
mu[1]; ks$x0n; sqrt(ks$P0n) # initial value info



## d) -------------------------------------------------
# start figure
png(filename="images/plot03.png", width = 1000, height = 600)
par(mfrow = c( 2, 1))
# plot(Time, mu[-1], main = 'Predict', ylim = c(-5, 10))
# lines(Time ,y , col = 'green' )
# lines(ks$xp)
# lines(ks$xp + 2 * sqrt (ks$Pp), lty = 2, col = 4)
# lines(ks$xp - 2 * sqrt (ks$Pp), lty = 2, col = 4)
ks=Ksmooth0(num, y, A = 1, mu0 = 0, Sigma0 = 1, Phi = 1, cQ = 1, cR = 1)
plot(Time, mu[-1], main = 'Filter Q=1 R=1', ylim = c(-5, 13))
lines(Time ,y , col = 'green')
lines(ks$xf)
lines(ks$xf + 2 * sqrt (ks$Pf) ,lty = 2, col = 4)
lines(ks $ xf - 2 * sqrt (ks$Pf) , lty = 2 , col = 4 )
legend(1,13, legend=c('Hidden States', 'Prediction', 'Confidence'), col=c('black', 'green', 'blue'), lty=c(1, 1, 2), cex=0.8)
ks=Ksmooth0(num, y, A = 1, mu0 = 0, Sigma0 = 1, Phi = 1, cQ = 1, cR = 10)
plot(Time, mu[-1], main = 'Filter Q=1 R=10', ylim = c(-5, 13))
lines(Time ,y , col = 'green')
lines(ks$xf)
lines(ks$xf + 2 * sqrt (ks$Pf) ,lty = 2, col = 4)
lines(ks $ xf - 2 * sqrt ( ks$Pf ) , lty = 2 , col = 4 )
# plot(Time, mu[-1] , main = 'Smooth', ylim = c (-5 ,10))
# lines(Time, y, col = 'green')
# lines(ks$xs)
# lines(ks$xs + 2 * sqrt(ks$Ps), lty = 2, col = 4)
# lines(ks$xs - 2 * sqrt(ks$Ps), lty = 2, col = 4)
dev.off()
mu[1]; ks$x0n; sqrt(ks$P0n) # initial value info



## e) -------------------------------------------------
my_kalman_filter <- function(A_t, C_t, Q_t, R_t, m_0, P_0, x_T ){
  # initializaction
  bigT <- length(x_T) + 1
  m_t <-rep(0, bigT)
  m_t[1] <- m_0
  P_t <- -rep(0, bigT)
  P_t[1] <- P_0
  for(t in 2:(bigT - 1)){
    K_t <- P_t[(t-1)] * C_t * (1/((C_t * P_t[(t-1)]) * C_t + R_t))
    m_t[t] <- m_t[(t - 1)] + K_t * (x_T[(t - 1)] - C_t * m_t[(t - 1)])
    P_t[t] <- (1 - K_t * C_t) + P_t[(t - 1)]
    
    m_t[(t + 1)] <- A_t * m_t[t]
    P_t[(t + 1)] <- A_t * P_t[t] * A_t + Q_t 
  }
  return(list(predictions = m_t, P_t = P_t))
}

png(filename="images/plot04.png", width = 1000, height = 900)
par(mfrow = c( 3, 1))
# Comparison Q and R one to one ratio
ourkf <- my_kalman_filter(A_t = 1, C_t = 1, Q_t = 1, R_t = 1, m_0=0, P_0=1, x_T = y)
ks=Ksmooth0(num, y, A = 1, mu0 = 0, Sigma0 = 1, Phi = 1, cQ = 1, cR = 1)
plot(Time, mu[-1], main = 'Filter Q=1 R=1', ylim = c(-5, 13))
lines(Time ,y , col = 'green')
lines(ks$xf)
lines(ks$xf + 2 * sqrt (ks$Pf) ,lty = 2, col = 4)
lines(ks $ xf - 2 * sqrt (ks$Pf) , lty = 2 , col = 4 )
lines(ourkf$predictions[-1], col="red")
legend(1,13, legend=c('Hidden States', 'Prediction', 'Confidence', 'myprediction'), col=c('black', 'green', 'blue', 'red'), lty=c(1, 1, 2), cex=0.8)
# comparison with grater R ratio
ourkf <- my_kalman_filter(A_t = 1, C_t = 1, Q_t = 1, R_t = 10, m_0=0, P_0=1, x_T = y)
ks=Ksmooth0(num, y, A = 1, mu0 = 0, Sigma0 = 1, Phi = 1, cQ = 1, cR = 10)
plot(Time, mu[-1], main = 'Filter Q=1 R=10', ylim = c(-5, 13))
lines(Time ,y , col = 'green')
lines(ks$xf)
lines(ks$xf + 2 * sqrt (ks$Pf) ,lty = 2, col = 4)
lines(ks $ xf - 2 * sqrt (ks$Pf) , lty = 2 , col = 4 )
lines(ourkf$predictions[-1], col="red")
legend(1,13, legend=c('Hidden States', 'Prediction', 'Confidence', 'myprediction'), col=c('black', 'green', 'blue', 'red'), lty=c(1, 1, 2), cex=0.8)
# comparison greater Q ratio
ourkf <- my_kalman_filter(A_t = 1, C_t = 1, Q_t = 10, R_t = 1, m_0=0, P_0=1, x_T = y)
ks=Ksmooth0(num, y, A = 1, mu0 = 0, Sigma0 = 1, Phi = 1, cQ = 10, cR = 1)
plot(Time, mu[-1], main = 'Filter Q=10 R=1', ylim = c(-5, 13))
lines(Time ,y , col = 'green')
lines(ks$xf)
lines(ks$xf + 2 * sqrt (ks$Pf) ,lty = 2, col = 4)
lines(ks $ xf - 2 * sqrt (ks$Pf) , lty = 2 , col = 4 )
lines(ourkf$predictions[-1], col="red")
legend(1,13, legend=c('Hidden States', 'Prediction', 'Confidence', 'myprediction'), col=c('black', 'green', 'blue', 'red'), lty=c(1, 1, 2), cex=0.8)
par(mfrow=c(1,1))
dev.off()
