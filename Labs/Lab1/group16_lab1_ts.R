
# Assignment 1 ---------------------------------------------------------------------

# a)

library(ggplot2)
library(ggfortify)
library(forecast)

set.seed(12345)
require(RColorBrewer)
#pal<-brewer.pal(9,"Paired")
pal<-colors()

xt=arima.sim(list(order=c(2,0,0),ar=c(0,-0.8)),n=100,start.innov=c(0,0),n.start=2)


t<-1:100
cs=cos( (2*pi*t)/5 )

# apply the filter 

f_coefs=rep(0.2,5)
xt_filtered=filter(xt,filter=f_coefs,sides=1)
cs_filtered=filter(cs,filter = f_coefs,sides=1)

col1<-sample(pal,1)
col2<-sample(pal,1)

par(mfrow=c(2,1))
plot.ts(xt,col=col1,lwd=2,main=expression("x"[t]*"=-0.8x"[t-2]*"+w"[t]),
        panel.first=grid(25,25))
lines(xt_filtered)
plot.ts(cs,col=col2,lwd=2,main=expression("x"[t]*"=cos(2"*pi*"t/5)"),
        panel.first=grid(25,25))
lines(cs_filtered)

#autoplot(xt,col=sample(pal,1))+autolayer(xt_filtered,col=samle(pal,1))

# b)

check_roots<-function(z){
  return(sqrt(Re(z)^2+Im(z)^2 ))
}

inv_caus_func<-function(AR_operator,MA_operator){
  n1<-length(AR_operator)-1
  n2<-length(MA_operator)-1
  
  res<-polyroot(AR_operator)
  res2<-polyroot(MA_operator)
  casuality<-sapply(res,function(y){check_roots(y)})
  #print(casuality>1)
  invert<-sapply(res2,function(y){check_roots(y)})
  #print(invert>1)
  #print(sum(casuality>1))
  #print(sum(invert>1))
  if( (sum(casuality>1)==n1) & (sum(invert>1)==n2) ){
    print("The series is invertible and casual")
  }else if(sum(casuality>1)==n1 ){
    print("The series is casal only!")
  }else if( sum(invert>1)==n2){
    print("The series is invertible only !")
  }else{
    print("Not casual or invertible ")
  }
}


ar_operator= c(1,-4,2,0,0,1) ; ma_operator= c(1,0,3,0,1,0,-1)

inv_caus_func(ar_operator,ma_operator)

# c)
set.seed(54321)

sim.ts<-arima.sim(model = list(order = c(1, 0, 2), ar=-3/4 ,ma=c(0,-1/9) ), n = 100)


plot.ts(sim.ts,col=sample(pal,1),panel.first=grid(25,25),lwd=2)
acf(sim.ts,main="ACF Plot",panel.first=grid(25,25),col=sample(pal,1),
    lwd=2)
acf(ARMAacf(ar=-3/4 ,ma=c(0,-1/9),lag=20),main="Theoretical ACF",
    col=sample(pal,1),lwd=2)



# Asiignment 2 ---------------------------------------------------------------------

# a)

library(astsa)

rhine<-read.csv2("Rhine.csv",sep=";")

rhine_ts<-ts(rhine[,4],start=c(rhine[1,1],1),end=c(2002,12),frequency=12)


plot.ts(rhine_ts,col=sample(pal,1),main="Rhine River TS",
        lwd=2,panel.first=grid(25,25))

lag1.plot(rhine_ts,12)
lag.plot(rhine_ts,12)


# b)
fit<-lm(rhine_ts~time(rhine_ts),na.action = NULL)

summary(fit)


col3<-sample(pal,1) ; col4<-sample(pal,1)

plot.ts(rhine_ts,col=col3,
        main="TS with fitted LR",panel.first=grid(25,25),lwd=2)
abline(fit,col=col4,lwd=2)
legend(130,7, legend=c("TS Rhine", "LR"),
       col=c(col3, col4), lty=1, cex=0.8,
       text.font=4, bg='white',lwd=2)


par(mfrow=c(1,2))
plot(residuals(fit),main="Detrended Rhine\n with LR",col=sample(pal,1),
     panel.first=grid(15,25),lwd=2)
acf(residuals(fit), 28,main="ACF with LR",panel.first=grid(25,25))

# c)

kernelSmooth<-ksmooth(rhine$Time, rhine_ts, "normal", bandwidth=4)

col5<-sample(pal,1)
plot.ts(rhine_ts,col=col5,
        main="TS with Kernel Smooth",panel.first=grid(25,25),lwd=2)
lines(kernelSmooth,col="red",lwd=2,lty=2)
legend(130,7, legend=c("TS Rhine", "KS"),
       col=c(col5, "red"), lty=1, cex=0.8,
       text.font=4, bg='white',lwd=2)

res<-(rhine_ts-kernelSmooth$y)

par(mfrow=c(1,2))
plot(res,main="Detrended Rhine\n with KS",col=sample(pal,1),
     panel.first=grid(15,25),lwd=2)
acf(res,main="ACF with KS",panel.first=grid(25,25))


# d)

new_rhine<-cbind(rhine,month_enc=c("January","February","March","April","May","June","July","August",
                           "September","October","November","December"))
new_rhine$month_enc<-as.factor(new_rhine$month_enc)
str(new_rhine)

seasonal.means<-lm(TotN_conc~month_enc+Time,data=new_rhine)

par(mfrow=c(1,2))
plot(resid(seasonal.means),main="Detrended Rhine\n with Seasonal Mean",col=sample(pal,1),
     panel.first=grid(15,25),lwd=2,type = "l")
acf(seasonal.means$residuals,main="ACF with Seasonal Means",panel.first=grid(25,25))

extractAIC(seasonal.means)

# e)

temp <- step(seasonal.means,direction = "both",trace = 0)
temp$anova$AIC
summary(temp)
extractAIC(temp)

library(MASS)
temp1<-stepAIC(seasonal.means,direction = "both")

#temp1$effects


#library(leaps)
#leaps<-regsubsets(TotN_conc~month_enc+Time,data=new_rhine,nbest=10)
# view results
#summary(leaps)
# plot a table of models showing variables in each model.
# models are ordered by the selection statistic.
#plot(leaps,scale="r2")
# plot statistic by subset size


# Assignment 3 ------------------------------------------------------------

# a)

library(astsa)

col3<-sample(pal,1) ; col4<-sample(pal,1)

plot.ts(gas,col=col3,lwd=2,main="Time series of oil and gas",
        panel.first=grid(25,25),ylim=c(0,350))
lines(oil,col=col4,lwd=2,
      panel.first=grid(25,25))
legend(2000,350, legend=c("gas", "oil"),
       col=c(col3, col4), lty=1, cex=0.8,
       title="TS types", text.font=4, bg='white',
       lwd=2)


# b)

logOil<-log(oil) ; logGas<-log(gas)

plot.ts(logGas,col=col3,lwd=2,main="Time series of log(oil) and log(gas)",
        panel.first=grid(25,25),ylim=c(0,6))
lines(logOil,col=col4,lwd=2,
      panel.first=grid(25,25))
legend(2000,2, legend=c("log(gas)", "log(oil)"),
       col=c(col3, col4), lty=1, cex=0.8,
       title="TS types", text.font=4, bg='white',
       lwd=2)

# c)

x_t<-diff(log(oil)) ; y_t<-diff(log(gas))
col5<-sample(pal,1) ; col6<-"#A6CEE3"

par(mfrow=c(1,2))
plot.ts(x_t, main="First difference Oil",
        col=col5,panel.first=grid(25,25),lwd=2)
plot.ts(y_t,col=col6,main="First difference Gas",
        panel.first=grid(25,25),lwd=2)

par(mfrow=c(1,2))
acf(x_t,main="ACF Oil",col=sample(pal,1),panel.first=grid(25,25),
    lwd=2)
acf(y_t,main="ACF Gas",col=sample(pal,1),panel.first=grid(25,25),
    lwd=2)

# d)
foo = expression("y"[t]*"~x"[t],"y"[t]*"~x"[t-1],
                 "y"[t]*"~x"[t-2],"y"[t]*"~x"[t-3])

par(mfrow=c(2,2))
for(i in 1:4){
  plot(lag(x_t,i),y_t,main=bquote("Scatterplot of"~.(foo[[i]])),
       col=sample(pal,1),panel.first=grid(25,25),pch=18)
  lines(ksmooth(lag(x_t,i),y_t,bandwidth=0.04,kernel="normal"),lwd=2,
        col="red",lty=2)
  #lines(lowess(y_t,lag(x_t,i),f=0.2),lwd=2)
  }

# e)

dummy=ifelse(x_t>0,0,1)
D<-ts.intersect(y_t,dL1=lag(x_t),x_t,dL2=lag(x_t,2),dframe=T)
summary(fit<-lm(y_t~dL1+x_t+dL2,data=D,na.action=NULL))

p_values<-summary(fit)$coefficients[,4]
p_values<0.05

plot(fit)

par(mfrow=c(1,2)) #,bg="lightcyan1")
plot.ts(fit$residuals,main="Plot of the residuals",
        col=sample(pal,1),panel.first=grid(25,25),ylab="Residuals",lwd=2)
acf(fit$residuals,main="ACF of REsiduals",col=sample(pal,1),
    panel.first=grid(25,25))



