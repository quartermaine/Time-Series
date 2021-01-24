
# Note  -------------------------------------------------------------------

# The code contains parts that I didn't use in the lab report just having for 
# I don't really know why...




# Assignment 1 ------------------------------------------------------------


#  a) ------------------------------------


set.seed(12345)

AR3.sim=arima.sim(list(ar=c(0.8,-0.2,0.1)), n=1000)

data1=ts.intersect(x=AR3.sim, x1=lag(AR3.sim,-1), x2=lag(AR3.sim,-2), x3=lag(AR3.sim,-3), dframe = T)

res1=lm(x~x1+x2,data=data1)

res2=lm(x3~x1+x2,data=data1)

resids1=residuals(res1)

resids2=residuals(res2)

estimatedCorr=cor(cbind(resids1,resids2))

theoreticalPacf=pacf(AR3.sim) 

sum_tab=cbind(estimatedCorr[1,2],theoreticalPacf[3][[1]]) 
colnames(sum_tab)=c("Estimated Corr","PACF")
sum_tab # they are close 

#  b) ------------------------------------

AR2.sim=arima.sim(list(ar=c(0.8,0.1)), n=100)

# Methods
AR_yw=ar(AR2.sim,order.max = 2,aic=F) # Yule-Walker
AR_ols=ar(AR2.sim,method="ols",order.max=2,aic=F) # OLS
AR_mle=ar(AR2.sim,method="mle",order.max=2,aic=F) # ML

d=rbind(AR_yw$ar,AR_ols$ar,AR_mle$ar,c(0.8,0.1)) ;rownames(d)=c("YW","OLS","MLE","TRUE")
colnames(d)=c("phi1","phi2") ; d
# check iF phi2 in CI for ML estimate
0.1%in%(c(AR_mle$ar[2]-sqrt(AR_mle$asy.var.coef[2,2])*1.96,AR_mle$ar[2]+sqrt(AR_mle$asy.var.coef[2,2])*1.96))

PACF = ARMAacf(ar=c(0.8,0.1), pacf=TRUE) ; PACF # theoretical PACF

# some plots of the series 
# plot.ts(AR2.sim,col="blue",lwd=2)
# lines(arima.sim(list(ar=AR_yw$ar),n=100),col="red",lwd=2)
# lines(arima.sim(list(ar=AR_mle$ar),n=100),col="green",lwd=2)
# lines(arima.sim(list(ar=AR_ols$ar),n=100),col="black",lwd=2)

#  c) ------------------------------------
library(astsa)

seasonal.sim=arima.sim(list(order=c(0,0,13),
                            ma=c(0.3,rep(0,10),0.6,(0.6*0.3))), n=200)
par(mfrow=c(1,2))
acf(seasonal.sim,main="Sample ACF",panel.first=grid(25,25))
pacf(seasonal.sim,main="Sample PACF",panel.first=grid(25,25))

ACF = ARMAacf(ma=c(0.3,rep(0,10),0.6,(0.6*0.3)))
PACF = ARMAacf(ma=c(0.3,rep(0,10),0.6,(0.6*0.3)), pacf=TRUE)

par(mfrow=c(1,2))
plot(ACF, type="h", xlab="Lag", ylim=c(-.4,.8),
     main="Theoretical ACF",panel.first=grid(25,25)) ;abline(h=0)
plot(PACF, type="h", xlab="Lag", ylim=c(-.4,.8),
     main="Theoretical PACF",panel.first=grid(25,25)) ; abline(h=0)

#  d) ------------------------------------
seasonal.sim1=arima.sim(list(order=c(0,0,13),
                             ma=c(0.3,rep(0,10),0.6,(0.6*0.3))), n=200)

# forecast using the sarima.for
fore.sar=sarima.for(seasonal.sim1,n.ahead=30,0,0,1,0,0,1,12)
# forecast using the 
fore = predict(arima(seasonal.sim1,order=c(0,0,1),
                     seasonal = list(order=c(0,0,1),period=12)), n.ahead=30)
cols=sample(colors(),2)
ts.plot(seasonal.sim1, fore$pred, col=cols,
        ylab="ARMA(0,0,1)x(0,0,1)"[12],lwd=3)
U = fore$pred+fore$se; L = fore$pred-fore$se
xx = c(time(U), rev(time(U))); yy = c(L, rev(U))
polygon(xx, yy, border = 8, col = gray(.6, alpha = .2))
points(fore$pred, pch=20, col="red")
legend("topleft",legend=c("Sim TS","Predicted 30"),
       col=c(cols[1],cols[2]),lty=1,lwd=2)

# using kernlab now

kernel_fit=kernlab::gausspr(x=c(1:200),seasonal.sim1)

kernels_preds=kernlab::predict(kernel_fit,c(200:230))

cols2=sample(colors(),2)
ts.plot(seasonal.sim1, ts(kernels_preds,start = 200,end=230), col=cols2,
        gpars=list(main="ARMA(0,0,1)x(0,0,1)"[12]),lwd=3)
#U4 = kernels_preds+fore$se; L4= kernels_preds-fore$se
#xx4 = c(time(U4), rev(time(U4))); yy4 = c(L4, rev(U4))
#polygon(xx4, yy4, border = 8, col = gray(.6, alpha = .2))
points(ts(kernels_preds,start=200,end=230), pch=20, col="red",cex=0.5)


#  e) ------------------------------------

arma.sim<-arima.sim(list(order=c(1,0,1),ar=0.7,ma=0.5),n=50)

sarima.for(arma.sim[1:40],n.ahead=10,1,0,1,0,0,0,0,no.constant = T)
# predictions with the sarima.for
fore1=predict(arima(arma.sim[1:40],order=c(1,0,1),include.mean = F),n.ahead = 10)

col1=sample(colors(),1) ; col2=sample(colors(),1) ; col3=sample(colors(),1)

ts.plot(as.ts(arma.sim[1:41]),fore1$pred,col=c(col1,col2),
        lwd=2,main="ARMA(1,1) with predictions")
lines(ts(arma.sim[41:50],start=41,end=50),col=col3,lwd=2,type="o")
U1 = fore1$pred+fore1$se; L1 = fore1$pred-fore1$se
xx1 = c(time(U1), rev(time(U1))); yy1 = c(L1, rev(U1))
polygon(xx1, yy1, border = 8, col = gray(.6, alpha = .2))
points(fore1$pred, pch=20, col=2,cex=0.5)
legend("bottomleft",legend=c("First 40 values","Predicted 10","True 10"),
       col=c(col1,col2,col3),lty=1,lwd=2)

mat<-cbind(as.vector(U1),as.vector(L1),as.vector(arma.sim[41:50]))
mat<-as.data.frame(mat) # ; mat


library(dplyr)

mat=mat%>%mutate(res=ifelse(( (mat$V3>mat$V1) | (mat$V3<mat$V2)),1,0 )) ; mat
sum(mat$res)

# Assignment 2 ------------------------------------------------------------


library(astsa)
require(RColorBrewer)


# data chicken -------------------

data(chicken)

diff.xt<-diff(chicken)

par(mfrow=c(2,2),bg = 'whitesmoke')
acf(chicken,60,col=sample(colors(),1),
    main=expression("ACF x"[t]),lwd=2)
pacf(chicken,60,col=sample(colors(),1),
     main=expression("PACFx"[t]),lwd=2)
acf(diff.xt,60,col=sample(colors(),1),
    main=expression("ACF diff x"[t]),lwd=2)
pacf(diff.xt,60,col=sample(colors(),1),
     main=expression("PACF diff x"[t]),lwd=2)
mtext(c("ACF~PACF Plots chicken", "ACF~PACF PLots for diff(chicken)"), side = 3,
      line = c(-2,-18), outer = TRUE)
# data so2 ----------------------

data(so2)

diff.xt<-diff(so2)

par(mfrow=c(2,2),bg = 'whitesmoke')
acf(chicken,40,col=sample(colors(),1),
    main=expression("ACF x"[t]))
pacf(chicken,40,col=sample(colors(),1),
     main=expression("PACFx"[t]))
acf(diff.xt,40,col=sample(colors(),1),
    main=expression("ACF diff x"[t]))
pacf(diff.xt,40,col=sample(colors(),1),
     main=expression("PACF diff x"[t]))
mtext(c("ACF~PACF Plots so2", "ACF~PACF PLots for diff(so2)"), side = 3,
      line = c(-2,-18), outer = TRUE)

# data EQcount ------------------

data("EQcount")


diff.xt<-diff(EQcount)

par(mfrow=c(2,2),bg = 'whitesmoke')
acf(chicken,40,col=sample(colors(),1),
    main=expression("ACF x"[t]))
pacf(chicken,40,col=sample(colors(),1),
     main=expression("PACFx"[t]))
acf(diff.xt,40,col=sample(colors(),1),
    main=expression("ACF diff x"[t]))
pacf(diff.xt,40,col=sample(colors(),1),
     main=expression("PACF diff x"[t]))
mtext(c("ACF~PACF Plots EQcount", "ACF~PACF PLots for diff(EQcount)"), side = 3,
      line = c(-2,-18), outer = TRUE)

# HCT ---------------------------

data(HCT)

diff.xt<-diff(HCT)

par(mfrow=c(2,2),bg = 'whitesmoke')
acf(chicken,40,col=sample(colors(),1),
    main=expression("ACF x"[t]),lwd=2)
pacf(chicken,40,col=sample(colors(),1),
     main=expression("PACFx"[t]),lwd=2)
acf(diff.xt,40,col=sample(colors(),1),
    main=expression("ACF diff x"[t]),lwd=2)
pacf(diff.xt,40,col=sample(colors(),1),
     main=expression("PACF diff x"[t]),lwd=2)
mtext(c("ACF~PACF Plots EQcount", "ACF~PACF PLots for diff(EQcount)"), side = 3,
      line = c(-2,-35), outer = TRUE)

# Assignment 3 ------------------------------------------------------------

#  a) ------------------------------------

library(ggplot2)

data(oil)
forecast::auto.arima(oil) # to automatically find the best model 

autoplot(oil)

oil%>%forecast::ggtsdisplay(main="TS-ACF-PACF for oil TS",
                              col=sample(colors(),1),theme=theme_gray())
oil%>%diff()%>%forecast::ggtsdisplay(main="TS-ACF-PACF for diff(oil) TS",
                                     col=sample(colors(),1),theme=theme_gray())
diff(log(oil))%>%forecast::ggtsdisplay(main="TS-ACF-PACF for log(diff(oil)) TS",
                                     col=sample(colors(),1),theme=theme_gray())

plot.ts(oil,col=sample(colors(),1),lwd=3,
        panel.first=grid(25,25),main="PLot of oil TS")
plot.ts(diff(oil),col=sample(colors(),1),lwd=3,
        panel.first=grid(25,25),main="PLot of diff(oil) TS")
plot.ts(diff(log(oil)),col=sample(colors(),1),lwd=3,
        panel.first=grid(25,25),main="PLot of diff(log(oil)) TS")


doil=diff(log(oil))
tseries::adf.test(doil)
plots=acf2(doil)

TSA::eacf(doil) # 2 choices AR(0,3) ARMA(1,1)

# start with ARMA(0,3)
fit=sarima(doil,0,0,3)
fit
#tsdiag(fit$fit)
resids1=residuals(fit$fit)

resids1%>%forecast::ggtsdisplay(main="TS-ACF-PACF for residuals",
                            col=sample(colors(),1),theme=theme_gray())

forecast::gghistogram(resids1,add.normal = T,add.kde = T) # plot of residuals with gghistogram 
hist(resids1,30,col=sample(colors(),1),border="red",
     panel.first=grid(25,25),freq=F,main="Plot of the residuals")
lines(density(resids1),lwd=2)

TSA::runs(resids1)

# then for ARMA(1,1)

fit1=sarima(doil,1,0,1)
fit1
resids2=residuals(fit1$fit)

resids2%>%forecast::ggtsdisplay(main="TS-ACF-PACF for residuals",
                            col=sample(colors(),1),theme=theme_gray())

hist(resids2,30,col=sample(colors(),1),border="red",
     panel.first=grid(25,25),freq=F,main="Plot of the residuals")
lines(density(resids2),lwd=2)

TSA::runs(resids2)

AIC_mat=cbind(c(fit$AIC,fit1$AIC),c(fit$BIC,fit1$BIC)) # ; AIC_mat
colnames(AIC_mat)=c("AIC","BIC") ; rownames(AIC_mat)=c("ARMA(0,3)","ARMA(1,1)")
AIC_mat

# we choose the ARMA(1,1) and make predictions 

sarima.for(oil,n.ahead=20,1,0,1,0,0,0,0)

fore2=predict(arima(oil,order=c(1,0,1),include.mean = F),n.ahead = 20)

cols4=sample(colors(),2)

ts.plot(oil,fore2$pred,col=cols4,
        lwd=3,main="ARMA(1,1) with predictions")
U2 = fore2$pred+fore2$se; L2 = fore2$pred-fore2$se
xx2 = c(time(U2), rev(time(U2))); yy2 = c(L2, rev(U2))
polygon(xx2, yy2, border = 8, col = gray(.6, alpha = .2))
points(fore2$pred, pch=19, col=2,cex=0.3)
legend("topleft",legend=c("Oil TS","Predicted Oil 20 ahead"),
       col=cols1,lty=1,lwd=2)


#  b) ------------------------------------


data(unemp)

forecast::auto.arima(unemp)

plot(decompose(unemp),col=sample(colors(),1))

# we take first diffrence 

plot(decompose(diff(unemp))) # the series looks more stationary

dunemp=diff(unemp)
plot(decompose(dunemp))
acf2(dunemp)

dx=diff(dunemp,12)
plot(decompose(dx))
acf2(dx)

tseries::adf.test(dx)
plots=acf2(dx) # arma(1,1)

TSA::eacf(dx) # ARMA(1,1) ARMA(0,3) with seasonal =12

# start with ARMA(2,2)  
fit3=sarima(dx,2,0,2)
fit3
#tsdiag(fit$fit)
resids3=residuals(fit3$fit)

resids3%>%forecast::ggtsdisplay(main="TS-ACF-PACF for residuals",
                                col=sample(colors(),1),theme=theme_gray())

forecast::gghistogram(resids3,add.normal = T,add.kde = T) # plot of residuals with gghistogram 
hist(resids3,30,col=sample(colors(),1),border="red",
     panel.first=grid(25,25),freq=F,main="Plot of the residuals")
lines(density(resids3),lwd=2)

TSA::runs(resids3)

# next with ARMA(0,3)12

fit4=sarima(dx,0,0,3)
fit4
#tsdiag(fit$fit)
resids4=residuals(fit4$fit)

resids4%>%forecast::ggtsdisplay(main="TS-ACF-PACF for residuals",
                                col=sample(colors(),1),theme=theme_gray())

forecast::gghistogram(resids4,add.normal = T,add.kde = T) # plot of residuals with gghistogram 

hist(resids4,30,col=sample(colors(),1),border="red",
     panel.first=grid(25,25),freq=F,main="Plot of the residuals")
lines(density(resids4),lwd=2)
rug(resids4,col="red")

TSA::runs(resids4)

AIC_mat1=cbind(c(fit3$AIC,fit4$AIC),c(fit3$BIC,fit4$BIC)) # ; AIC_mat
colnames(AIC_mat1)=c("AIC","BIC") ; rownames(AIC_mat1)=c("ARMA(1,1)","ARMA(0,3)")
AIC_mat1

# we choose the SARMA(0,1,0)x(2,0,2,)12 and make predictions 

sarima.for(unemp,n.ahead=20,1,1,1,2,1,2,12)

sarima.for(unemp,n.ahead = 20,1,1,1,0,1,3,12)

# fore4=predict(arima(unemp,order=c(0,0,0),seasonal=list(order=c(2,0,12)),include.mean = F),n.ahead = 20)

cols5=sample(colors(),2)

ts.plot(unemp,fore4$pred,col=cols5,
        lwd=3,main="ARMA(1,1) with predictions")
U2 = fore4$pred+fore4$se; L2 = fore4$pred-fore4$se
xx2 = c(time(U2), rev(time(U2))); yy2 = c(L2, rev(U2))
polygon(xx2, yy2, border = 8, col = gray(.6, alpha = .2))
points(fore4$pred, pch=19, col=2,cex=0.3)
legend("topleft",legend=c("Oil TS","Predicted Oil 20 ahead"),
       col=cols1,lty=1,lwd=2)




