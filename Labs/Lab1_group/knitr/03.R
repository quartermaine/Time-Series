col3<-sample(pal,1) ; col4<-sample(pal,1)

plot.ts(gas,col=col3,lwd=2,main="Time series of oil and gas",
        panel.first=grid(25,25),ylim=c(0,350))
lines(oil,col=col4,lwd=2,
      panel.first=grid(25,25))
legend(2000,350, legend=c("gas", "oil"),
       col=c(col3, col4), lty=1, cex=0.8,
       title="TS types", text.font=4, bg='white',
       lwd=2)

logOil<-log(oil) ; logGas<-log(gas)

plot.ts(logGas,col=col3,lwd=2,main="Time series of log(oil) and log(gas)",
        panel.first=grid(25,25),ylim=c(0,6))
lines(logOil,col=col4,lwd=2,
      panel.first=grid(25,25))
legend(2000,2, legend=c("log(gas)", "log(oil)"),
       col=c(col3, col4), lty=1, cex=0.8,
       title="TS types", text.font=4, bg='white',
       lwd=2)

x_t<-diff(log(oil)) ; y_t<-diff(log(gas))
col5<-sample(pal,1) ; col6<-"#A6CEE3"

par(mfrow=c(1,2))
plot.ts(x_t, main="First difference Oil",
        col=col5,panel.first=grid(25,25),lwd=2)
plot.ts(y_t,col=col6,main="First difference Gas",
        panel.first=grid(25,25),lwd=2)

par(mfrow=c(1,2))
acf(x_t,main="ACF Oil",col=sample(pal,1),panel.first=grid(25,25))
acf(y_t,main="ACF Gas",col=sample(pal,1),panel.first=grid(25,25))

foo = expression("y"[t]*"~x"[t],"y"[t]*"~x"[t-1],
                 "y"[t]*"~x"[t-2],"y"[t]*"~x"[t-3])

par(mfrow=c(2,2))
for(i in 1:4){
  plot(y_t,lag(x_t,i),main=bquote("Scatterplot of"~.(foo[[i]])),
       col=sample(pal,1),panel.first=grid(25,25),pch=18)
  lines(lowess(y_t, lag(x_t, i)),lwd=2)
}

dummy=ifelse(x_t>0,0,1)
D<-ts.intersect(y_t,dL1=lag(x_t),x_t,dL2=lag(x_t,2),dframe=T)
summary(fit<-lm(y_t~dL1+x_t+dL2,data=D,na.action=NULL))

p_values<-summary(fit)$coefficients[,4]
p_values<0.05

plot(fit)

par(mfrow=c(1,2))
plot.ts(fit$residuals,main="Plot of the residuals",
        col=sample(pal,1),panel.first=grid(25,25),ylab="Residuals",lwd=2)
acf(fit$residuals,main="ACF of REsiduals",col=sample(pal,1),
    panel.first=grid(25,25))
