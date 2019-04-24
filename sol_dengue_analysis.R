sol_dengue_data_ts<-ts(sol_dengue_data,start=c(2010,1),
                       frequency = 12)
plot(sol_dengue_data_ts)

#Exponential Smoothing
library(smooth)

dengue.auto <- es(sol_dengue_data_ts,model="ZZZ",h=12,holdout=FALSE,intervals="parametric",level=0.95,cfType="MSE")
dengue.ses  <- es(sol_dengue_data_ts,model="ANN",h=12,holdout=FALSE,intervals="parametric",level=0.95,cfType="MSE")
dengue.des  <- es(sol_dengue_data_ts,model="AAdN",h=12,holdout=FALSE,intervals="parametric",level=0.95,cfType="MSE")
dengue.hwn  <- es(sol_dengue_data_ts,model="AAN",h=12,holdout=FALSE,intervals="parametric",level=0.95,cfType="MSE")
dengue.hwa  <- es(sol_dengue_data_ts,model="AAA",h=12,holdout=FALSE,intervals="parametric",level=0.95,cfType="MSE")
dengue.hwm  <- es(sol_dengue_data_ts,model="AAM",h=12,holdout=FALSE,intervals="parametric",level=0.95,cfType="MSE")


summary(dengue.auto)

plot(forecast(dengue.hwa,h=68))

#ARIMA

library(forecast)

decomp_dengue<-decompose(sol_dengue_data_ts,"additive")
plot(decomp_dengue)

dengue.auto.arima<-auto.arima(sol_dengue_data_ts,max.order = 15)

dengue.auto.pred<-forecast(dengue.auto.arima, h = 48, level=c(97.5))

plot(dengue.auto.pred)

#Spectral Analysis

sol_dengue_data<-as.numeric(sol_dengue_data_ts)

P = abs(fft(sol_dengue_data_ts)); Fr = 0:39
plot(Fr, P, type="o", xlab="frequency", ylab="periodogram")
f<-abs(fft(sol_dengue_data_ts))
a<-which(f>100)
b<-(a-1)*2

t<-1:40
x=y=matrix(nr=length(t),nc=10)
x[,1]<-cos(2*pi*t/12)
y[,1]<-sin(2*pi*t/12)
x[,2]<-cos(4*pi*t/12)
y[,2]<-sin(4*pi*t/12)
x[,3]<-cos(6*pi*t/12)
y[,3]<-sin(6*pi*t/12)
x[,4]<-cos(8*pi*t/12)
y[,4]<-sin(8*pi*t/12)
x[,5]<-cos(10*pi*t/12)
y[,5]<-sin(10*pi*t/12)

dengue.spec1<-lm(sol_dengue_data~x[,1]+y[,1]+x[,2]+y[,2]
                 +x[,3]+y[,3]+x[,4]+y[,4]+x[,5]+y[,5])


fv<-dengue.spec1$fitted.values
plot(sol_dengue_data,type="l")
lines(1:40,fv,col="red")

summary(dengue.spec1)
