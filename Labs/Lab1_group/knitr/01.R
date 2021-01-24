library(ggplot2)
library(ggfortify)
# autoplot(USAccDeaths)
library(forecast)

set.seed(54321)

require(RColorBrewer)
pal<-brewer.pal(9,"Paired")

xt=arima.sim(list(order=c(2,0,0),ar=c(0,-0.8)),n=100,start.innov=c(0,0),n.start=2)

t<-1:100
cs=cos( (2*pi*t)/5 )

# apply the filter 

f_coefs=rep(0.2,5)
xt_filtered=filter(xt,filter=f_coefs,sides=1)
cs_filtered=filter(cs,filter = f_coefs,sides=1)

col1<-sample(pal,1)
col2<-sample(pal,1)




check_causality<-function(z){
  return(sqrt(Re(z)^2+Im(z)^2 ))
}

inv_caus_func<-function(AR_operator,MA_operator){
  n1<-length(AR_operator)-1
  n2<-length(MA_operator)-1
  
  res<-polyroot(AR_operator)
  res2<-polyroot(MA_operator)
  casuality<-sapply(res,function(y){check_causality(y)})
  #print(casuality>1)
  invert<-sapply(res2,function(y){check_causality(y)})
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

