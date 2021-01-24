## ----setup, include=FALSE------------------------------------------------
# knitr::opts_chunk$set(echo = F, message = F,warning = F)


# Libraries ---------------------------------------------------------------
library(astsa)
library(ggplot2)
library(knitr)
# -------------------------------------------------------------------------

# Assignment 1 ------------------------------------------------------------

#  a) ------------------------------------
set.seed(12345)
# sim AR(3) n=1000
AR3.sim=arima.sim(list(ar=c(0.8,-0.2,0.1)), n=1000)
# bind ts for up to 3 lags
data1=ts.intersect(x=AR3.sim, x1=lag(AR3.sim,-1), x2=lag(AR3.sim,-2), x3=lag(AR3.sim,-3), dframe = T)
# linear regressions
res1=lm(x~x1+x2,data=data1)
res2=lm(x3~x1+x2,data=data1)
# calculate residuals
resids1=residuals(res1)
resids2=residuals(res2)
# calculate correlation for the residuals
estimatedCorr=cor(cbind(resids1,resids2))
# theoretical pacf
theoreticalPacf=ARMAacf(ar=c(0.8,-0.2,0.1),pacf=T)
simulatedPacf=pacf(AR3.sim)
# make a table with the corr and pacf
sum_tab=cbind(estimatedCorr[1,2],simulatedPacf[3][[1]],theoreticalPacf[3]) 
colnames(sum_tab)=c("Est.Corr","Sim.PACF","Theo.PACF")



## ----------------------------
kable(sum_tab)
# ----------------------------


## ---------------------------
#  b) ------------------------------------
# simulate AR(2) n=100
AR2.sim=arima.sim(list(ar=c(0.8,0.1)), n=100)
# perform each method 
AR_yw=ar(AR2.sim,order.max = 2,aic=F) # Yule-Walker
AR_ols=ar(AR2.sim,method="ols",order.max=2,aic=F) # OLS
AR_mle=ar(AR2.sim,method="mle",order.max=2,aic=F) # ML
# make a table with the coefficients
d=rbind(AR_yw$ar,AR_ols$ar,AR_mle$ar,c(0.8,0.1)) ;rownames(d)=c("YW","OLS","MLE","TRUE")
colnames(d)=c("phi1","phi2") 
kable(d) 



## --------------------------------
# check iF phi2 in CI for ML estimate
# compute the CI  
lower_CI=AR_mle$ar[2]-sqrt(AR_mle$asy.var.coef[2,2])*1.96 # lower limit
upper_CI=AR_mle$ar[2]+sqrt(AR_mle$asy.var.coef[2,2])*1.96 # upper limit
# check if the phi2=0.1 in the CI 
cat("The esitmated interval is :",c(lower_CI,upper_CI),"\n")
cat("The value of phi_2 is within the estimated interval?\n")
0.1%in%round(seq(as.numeric(lower_CI[1]),as.numeric(upper_CI[1]),length.out = 1000),2) # returns T,F
# 
# PACF = ARMAacf(ar=c(0.8,0.1), pacf=TRUE) ; PACF # theoretical PACF
# -------------------------------------------------------------------------

#  c) ------------------------------------
# simulate ARIMA(0,0,1)x(0,0,1)12 n=200
seasonal.sim=arima.sim(list(order=c(0,0,13),
                            ma=c(0.3,rep(0,10),0.6,(0.6*0.3))), n=200)
# plot of the simulated PACF and ACF for simulation
par(mfrow=c(2,2),bg="whitesmoke")
acf(seasonal.sim,main="Sample ACF",panel.first=grid(25,25),lag.max=40)
pacf(seasonal.sim,main="Sample PACF",panel.first=grid(25,25),lag.max=40)
# compute theoretical ACF and PACF
ACF = ARMAacf(ma=c(0.3,rep(0,10),0.6,(0.6*0.3)),lag.max = 40)
PACF = ARMAacf(ma=c(0.3,rep(0,10),0.6,(0.6*0.3)), pacf=TRUE,lag.max = 40)
# plot of the theoretical ACF and PACF  
plot(ACF, type="h", xlab="Lag", ylim=c(-.4,.8),
     main="Theoretical ACF",panel.first=grid(25,25)) ;abline(h=0)
plot(PACF, type="h", xlab="Lag", ylim=c(-.4,.8),
     main="Theoretical PACF",panel.first=grid(25,25)) ; abline(h=0)

# -------------------------------------------------------------------------

#  d) ------------------------------------

# using simulation ----------------

# simulate ARIMA(O,0,1)x(0,0,1)12
seasonal.sim1=arima.sim(list(order=c(0,0,13),
                             ma=c(0.3,rep(0,10),0.6,(0.6*0.3))), n=200)
# forecast using the sarima.for and auto plot 
# fore.sar=sarima.for(seasonal.sim1,n.ahead=30,p=0,d=0,q=1,P=0,D=0,Q=1,S=12)

# the same results but using predict
# forecast using the predict and plot of the forecasts
fore = predict(arima(seasonal.sim1,order=c(0,0,1),
                     seasonal = list(order=c(0,0,1),period=12)), n.ahead=30)
# plot of ts and predictions and band
cols=c("mediumpurple3" ,"darkslateblue")
par(mfrow=c(1,1),oma = c( 0, 0, 3, 0 ))
ts.plot(seasonal.sim1, fore$pred, col=cols,lwd=3)
U = fore$pred+fore$se; L = fore$pred-fore$se
xx = c(time(U), rev(time(U))); yy = c(L, rev(U))
polygon(xx, yy, border = 8, col = gray(.6, alpha = .2))
points(fore$pred, pch=20, col="red")
legend("topright",legend=c("Sim TS","Predicted 30"),
         col=c(cols[1],cols[2]),lty=1,lwd=2)
title(expression("ARMA(0,0,1)x(0,0,1)"[12]))

# using the kernelab --------------
# kenel fit not import the kernelab because of the predict 
kernel_fit=kernlab::gausspr(x=c(1:200),seasonal.sim1)
# make predictions
kernels_preds=kernlab::predict(kernel_fit,c(200:230))
# plot the ts and predictions
cols1=c("blue2","black")
par(mfrow=c(1,1),oma=c(0,0,3,0))
ts.plot(seasonal.sim1, ts(kernels_preds,start = 200,end=230), col=cols1,lwd=3)
points(ts(kernels_preds,start=200,end=230), pch=20, col="azure3",cex=0.2)
legend("bottomright",legend=c("Sim TS","Predicted 30"),
       col=c(cols1[1],cols[2]),lty=1,lwd=2)
title(expression("ARMA(0,0,1)x(0,0,1)"[12]*" kernelab"))

# -------------------------------------------------------------------------

#  e) ------------------------------------
# simulate ARMA(1,1) n=50
arma.sim<-arima.sim(list(order=c(1,0,1),ar=0.7,ma=0.5),n=50)
# make predictions with the sarima.for
# sarima.for(arma.sim[1:40],n.ahead=10,1,0,1,0,0,0,0,no.constant = T)

# the same results but using predict
# forecast using the predict and plot of the forecasts
fore1=predict(arima(arma.sim[1:40],order=c(1,0,1),include.mean = F),n.ahead = 10)
col1="forestgreen" ; col2="magenta1" ; col3="orange1"
ts.plot(as.ts(arma.sim[1:41]),fore1$pred,col=c(col1,col2),
        lwd=2,main="ARMA(1,1) with predictions")
lines(ts(arma.sim[41:50],start=41,end=50),col=col3,lwd=2,type="o")
U1 = fore1$pred+fore1$se; L1 = fore1$pred-fore1$se
xx1 = c(time(U1), rev(time(U1))); yy1 = c(L1, rev(U1))
polygon(xx1, yy1, border = 8, col = gray(.6, alpha = .2))
points(fore1$pred, pch=20, col="black",cex=0.5)
legend("topleft",legend=c("First 40 values","Predicted 10","True 10"),
       col=c(col1,col2,col3),lty=1,lwd=2)


## -------------------
# Note
require(dplyr) # we load the library here beacaus it masks the lag and had problems in Ass1.a
# count how many points are not in the band
mat<-cbind(as.vector(U1),as.vector(L1),as.vector(arma.sim[41:50]))
mat<-as.data.frame(mat) # ; mat
mat=mat%>%mutate(res=ifelse(( (mat$V3>mat$V1) | (mat$V3<mat$V2)),1,0 )) 
cat("The number of the values outside the prediction band is:",sum(mat$res))

# -------------------------------------------------------------------------

# Assignment 2 ------------------------------------------------------------
# data chicken -------------------
data(chicken)
# girst diffrence
diff.xt<-diff(chicken)
my_colors=c("seagreen4","red4")
# ACF and PACF plots
par(mfrow=c(2,2),bg = 'whitesmoke',oma=c(0,0,3,0))
acf(chicken,60,col=my_colors[1],
    main=expression("ACF x"[t]),lwd=2)
pacf(chicken,60,col=my_colors[1],
     main=expression("PACFx"[t]),lwd=2)
acf(diff.xt,60,col=my_colors[2],
    main=expression("ACF diff x"[t]),lwd=2)
pacf(diff.xt,60,col=my_colors[2],
     main=expression("PACF diff x"[t]),lwd=2)
title("\nACF and PACF Plots for chicken \nwith first diffrences",outer = TRUE )
#mtext(c("ACF~PACF Plots chicken", "ACF~PACF PLots for diff(chicken)"), side = 3,
#     line = c(-2,-18), outer = TRUE)

# -------------------------------

# data so2 ----------------------
data(so2)
# first diffrence
diff.xt<-diff(so2)
# ACF and PACF plots
par(mfrow=c(2,2),bg = 'whitesmoke',oma=c(0,0,3,0))
acf(chicken,40,col=my_colors[1],
    main=expression("ACF x"[t]),lwd=2)
pacf(chicken,40,col=my_colors[1],
     main=expression("PACFx"[t]),lwd=2)
acf(diff.xt,40,col=my_colors[2],
    main=expression("ACF diff x"[t]),lwd=2)
pacf(diff.xt,40,col=my_colors[2],
     main=expression("PACF diff x"[t]),lwd=2)
title("\nACF and PACF Plots for so2 \nwith first diffrence",outer = TRUE )
#mtext(c("ACF~PACF Plots so2", "ACF~PACF PLots for diff(so2)"), side = 3,
#      line=c(0,-19),outer = T)

# -------------------------------

# data EQcount ------------------
data("EQcount")
# first diffrence
diff.xt<-diff(EQcount)
# plots 
par(mfrow=c(2,2),bg = 'whitesmoke',oma=c(0,0,3,0))
acf(chicken,40,col=my_colors[1],
    main=expression("ACF x"[t]),lwd=2)
pacf(chicken,40,col=my_colors[1],
     main=expression("PACFx"[t]),lwd=2)
acf(diff.xt,40,col=my_colors[2],
    main=expression("ACF diff x"[t]),lwd=2)
pacf(diff.xt,40,col=my_colors[2],
     main=expression("PACF diff x"[t]),lwd=2)
title("\nACF and PACF Plots for EQcount \nwith first diffrence",outer = TRUE )
#mtext(c("ACF~PACF Plots EQcount", "ACF~PACF PLots for diff(EQcount)"), side = 3,
#      line = c(-2,-18), outer = TRUE)

# ---------------------------------

# HCT -----------------------------
data(HCT)
# first diffrence
diff.xt<-diff(HCT)
# plots
par(mfrow=c(2,2),bg = 'whitesmoke',oma=c(0,0,3,0))
acf(chicken,40,col=my_colors[1],
    main=expression("ACF x"[t]),lwd=2)
pacf(chicken,40,col=my_colors[1],
     main=expression("PACFx"[t]),lwd=2)
acf(diff.xt,40,col=my_colors[2],
    main=expression("ACF diff x"[t]),lwd=2)
pacf(diff.xt,40,col=my_colors[2],
     main=expression("PACF diff x"[t]),lwd=2)
title("\nACF and PACF Plots for HTC \nwith first diffrence",outer = TRUE )
#mtext(c("ACF~PACF Plots EQcount", "ACF~PACF PLots for diff(EQcount)"), side = 3,
#      line = c(-2,-35), outer = TRUE)

# -----------------------------------

# Assignment 3 ------------------------------------------------------------

#  a) ------------------------------------
data(oil)
# plots 
par(mfrow=c(2,2),oma=c(0,0,3,0),bg="whitesmoke")
plot.ts(oil,col="darkblue",lwd=2,
        panel.first=grid(25,25),main="PLot of oil TS")
plot.ts(diff(oil),col="darkred",lwd=2,
        panel.first=grid(25,25),main="PLot of diff(oil) TS")
plot.ts(diff(log(oil)),col="purple3",lwd=2,
        panel.first=grid(25,25),main="PLot of diff(log(oil)) TS")
title("TS Plots",outer=T)



## ---------------------------------------------
par(mfrow=c(2,2),bg="whitesmoke",oma=c(0,0,3,0))
acf(oil);acf(diff(oil)) ; acf(diff(log(oil)))
title("ACF Plots",outer = T)


## ---------------------------------------------
par(mfrow=c(2,2),bg="whitesmoke",oma=c(0,0,3,0))
pacf(oil);pacf(diff(oil));pacf(diff(log(oil)))
title("PACF Plots",outer = T)

## ---------------------------------------------
# we work with diff(log(oil)) for the analysis
doil=diff(log(oil))
# test the p-value 
tseries::adf.test(doil)

plots=acf2(doil)

## ------------------------------------------------------------------------
# eacf test 
TSA::eacf(doil) # 2 choices AR(0,3) ARMA(1,1)

# Start with ARMA(0,3) ---------

# fit the model
arma.fit=sarima(doil,p=0,d=0,q=3,details=T)
# tsdiag(arma.fit$fit) diagnostic plots 
resids1=residuals(arma.fit$fit) # calculate residuals


# plot the residuals withe the forecast package 
#resids1%>%forecast::ggtsdisplay(main="TS-ACF-PACF for residuals",
#                                col=sample(colors(),1),theme=theme_gray())

#forecast::gghistogram(resids1,add.normal = T,add.kde = T) # plot of residuals with gghistogram 
# plot the histtogram of the residuals with basic
hist(resids1,30,col=sample(colors(),1),border="red",
     panel.first=grid(25,25),freq=F,main="Plot of the residuals for ARMA(0,3)")
lines(density(resids1),lwd=2)
rug(resids1,col="red4")



## --------------------------------------------
# Test the independence of a sequence of random variables
kable(unlist(TSA::runs(resids1)))

# Then with ARMA(1,1) ---------
# fit model 
arma.fit1=sarima(doil,1,0,1)
resids2=residuals(arma.fit1$fit)

# plot of residuals with the forecast package
#resids2%>%forecast::ggtsdisplay(main="TS-ACF-PACF for residuals",
#                                col=sample(colors(),1),theme=theme_gray())
# plot of the histogram with the basic
hist(resids2,30,col=sample(colors(),1),border="red",
     panel.first=grid(25,25),freq=F,main="Plot of the residuals for ARMA(1,1)")
lines(density(resids2),lwd=2)
rug(resids2,col="red4")



## ----------------------------------
kable(unlist(TSA::runs(resids2)))


## -------------------------------------
# ckeck the AIC and BIC for each model
AIC_BIC_mat=cbind(c(arma.fit$AIC,arma.fit1$AIC),c(arma.fit$BIC,arma.fit1$BIC)) # ; AIC_mat
colnames(AIC_BIC_mat)=c("AIC","BIC") ; rownames(AIC_BIC_mat)=c("ARMA(0,3)","ARMA(1,1)")
kable(AIC_BIC_mat)

# we choose the ARMA(1,1) and make predictions 
par(mfrow=c(1,1),bg="whitesmoke")
S1<-sarima.for(oil,n.ahead=20,p=1,d=1,q=1,P=0,D=0,Q=0,S=0)

# the same results but using predict 
#fore2=predict(arima(oil,order=c(1,0,1),include.mean = F),n.ahead = 20)

#cols4=sample(colors(),2)

#ts.plot(oil,fore2$pred,col=cols4,
#        lwd=3,main="ARMA(1,1) with predictions")
#U2 = fore2$pred+fore2$se; L2 = fore2$pred-fore2$se
#xx2 = c(time(U2), rev(time(U2))); yy2 = c(L2, rev(U2))
#polygon(xx2, yy2, border = 8, col = gray(.6, alpha = .2))
#points(fore2$pred, pch=19, col=2,cex=0.3)
#legend("topleft",legend=c("Oil TS","Predicted Oil 20 ahead"),
#       col=cols1,lty=1,lwd=2)


#  b) ------------------------------------

## -------------------------------------
data(unemp)
# first difference
dunemp=diff(unemp)
# remove seasonal by 12 difference
dx=diff(dunemp,12)

# plots 
par(mfrow=c(2,2),oma=c(0,0,3,0),bg="whitesmoke")
plot.ts(unemp,col="darkblue",lwd=2,
        panel.first=grid(25,25),main="PLot of unemp TS")
plot.ts(diff(oil),col="darkred",lwd=2,
        panel.first=grid(25,25),main="PLot of diff(unemp) TS")
plot.ts(dx,col="purple3",lwd=2,
        panel.first=grid(25,25),main="PLot of diff(diff(unemp)),12) TS")
title("TS Plots ",outer=T)


## ------------------------------------------------------------------------
par(mfrow=c(1,1),bg="whitesmoke")
plot(decompose(dx),col="darkblue")



## ------------------------------------------------------------------------
par(mfrow=c(2,2),bg="whitesmoke",oma=c(0,0,3,0))
acf(unemp);acf(dunemp) ; acf(dx)
title("ACF Plots",outer = T)



## ------------------------------------------------------------------------
par(mfrow=c(2,2),bg="whitesmoke",oma=c(0,0,3,0))
pacf(unemp);pacf(dunemp);pacf(dx)
title("PACF Plots",outer = T)


## ---------------------------------------
# we work with dx for the analysis
# test the p-value 
tseries::adf.test(dx)



## ---------------------------------------
# eacf test 
TSA::eacf(dx) # 2 choices AR(4,0) ARMA(2,2)

# Start with ARMA(4,0) ---------

# fit the model
arma.fit3=sarima(dx,p=4,d=0,q=0,details=T)
# tsdiag(arma.fit$fit) diagnostic plots 
resids3=residuals(arma.fit3$fit) # calculate residuals

# plot the residuals withe the forecast package 
#resids1%>%forecast::ggtsdisplay(main="TS-ACF-PACF for residuals",
#                                col=sample(colors(),1),theme=theme_gray())

#forecast::gghistogram(resids1,add.normal = T,add.kde = T) # plot of residuals with gghistogram 
# plot the histtogram of the residuals with basic
hist(resids3,30,col=sample(colors(),1),border="red",
     panel.first=grid(25,25),freq=F,main="Plot of the residuals for ARMA(4,0)")
lines(density(resids3),lwd=2)
rug(resids3,col="red4")

## ---------------------------------------------
kable(unlist(TSA::runs(resids3)))

# Start with ARMA(0,3) ---------

# fit the model
arma.fit4=sarima(dx,p=2,d=0,q=2,details=T)
# tsdiag(arma.fit$fit) diagnostic plots 
resids4=residuals(arma.fit4$fit) # calculate residuals

# plot the residuals withe the forecast package 
#resids1%>%forecast::ggtsdisplay(main="TS-ACF-PACF for residuals",
#                                col=sample(colors(),1),theme=theme_gray())

#forecast::gghistogram(resids1,add.normal = T,add.kde = T) # plot of residuals with gghistogram 
# plot the histtogram of the residuals with basic
hist(resids4,30,col=sample(colors(),1),border="red",
     panel.first=grid(25,25),freq=F,main="Plot of the residuals for ARMA(0,3)")
lines(density(resids4),lwd=2)
rug(resids4,col="red4")



## -----------------------------------
kable(unlist(TSA::runs(resids4)))


## -------------------------------------
# ckeck the AIC and BIC for each model
AIC_BIC_mat1=cbind(c(arma.fit3$AIC,arma.fit4$AIC),c(arma.fit3$BIC,arma.fit4$BIC)) # ; AIC_mat
colnames(AIC_BIC_mat1)=c("AIC","BIC") ; rownames(AIC_BIC_mat1)=c("ARMA(4,0)","ARMA(2,2)")
kable(AIC_BIC_mat1)

# we choose the ARMA(1,1) and make predictions 
par(mfrow=c(1,1),bg="whitesmoke")
S2<-sarima.for(unemp,n.ahead=20,p=1,d=1,q=1,P=4,D=1,Q=0,S=12)

# the same results but using predict 
#fore2=predict(arima(oil,order=c(1,0,1),include.mean = F),n.ahead = 20)

#cols4=sample(colors(),2)

#ts.plot(oil,fore2$pred,col=cols4,
#        lwd=3,main="ARMA(1,1) with predictions")
#U2 = fore2$pred+fore2$se; L2 = fore2$pred-fore2$se
#xx2 = c(time(U2), rev(time(U2))); yy2 = c(L2, rev(U2))
#polygon(xx2, yy2, border = 8, col = gray(.6, alpha = .2))
#points(fore2$pred, pch=19, col=2,cex=0.3)
#legend("topleft",legend=c("Oil TS","Predicted Oil 20 ahead"),
#       col=cols1,lty=1,lwd=2)
    
## ----ref.label=knitr::all_labels(), echo = T, eval = F-------------------


