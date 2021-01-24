library(astsa)
rhine<-read.csv2("data/Rhine.csv",sep=";")
rhine_ts<-ts(rhine[,4],start=c(rhine[1,1],1),end=c(2002,12),frequency=12)





fit<-lm(rhine_ts~time(rhine_ts),na.action = NULL)

# summary(fit)

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

kernelSmooth<-ksmooth(time(rhine_ts), rhine_ts, "normal", bandwidth=0.2)

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
acf(res,12,main="ACF with KS",panel.first=grid(25,25))



new_rhine<-cbind(rhine,month_enc=c("January","February","March","April","May","June","July","August",
                                   "September","October","November","December"))
new_rhine$month_enc<-as.factor(new_rhine$month_enc)
str(new_rhine)

seasonal.means<-lm(TotN_conc~month_enc+Time,data=new_rhine)

# plot(seasonal.means$)
# 
par(mfrow=c(1,2))
plot(seasonal.means$residuals,main="Detrended Rhine\n with Seasonal Mean",col=sample(pal,1),
     panel.first=grid(15,25),lwd=2,type = "l")
acf(seasonal.means$residuals,12,main="ACF with Seasonal Means",panel.first=grid(25,25))

temp <- step(seasonal.means,direction = "both",trace = 0, steps = 1)
temp$anova$AIC
summary(temp)

library(MASS)
temp1<-stepAIC(seasonal.means,direction = "both")
temp1$effects
