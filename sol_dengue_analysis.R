sol_dengue_data_ts<-ts(sol_dengue_data$X1,start=c(2010,1),
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

#Spectral Analysis and FSHA

sol_dengue_data<-as.numeric(sol_dengue_data_ts)

P = abs(fft(sol_dengue_data_ts)); Fr = 0:39
plot(Fr, P, type="o", xlab="frequency", ylab="periodogram")
f<-abs(fft(sol_dengue_data_ts))
a<-which(f>100)
b<-(a-1)*2

t<-1:40
x=y=matrix(nr=length(t),nc=11)
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
x[,6]<-cos(pi*t)

dengue.fsha.1<-lm(sol_dengue_data~x[,1]+y[,1]+x[,2]+y[,2]
                 +x[,3]+y[,3]+x[,4]+y[,4]+x[,5]+y[,5]+x[,6])
#x[,6]
dengue.fsha.2<-lm(sol_dengue_data~x[,1]+y[,1]+x[,2]+y[,2]
                  +x[,3]+y[,3]+x[,4]+y[,4]+x[,5]+y[,5])
#x[,4]
dengue.fsha.3<-lm(sol_dengue_data~x[,1]+y[,1]+x[,2]+y[,2]
                  +x[,3]+y[,3]+y[,4]+x[,5]+y[,5])
#y[,5]
dengue.fsha.4<-lm(sol_dengue_data~x[,1]+y[,1]+x[,2]+y[,2]
                  +x[,3]+y[,3]+y[,4]+x[,5])
#x[,3]
dengue.fsha.5<-lm(sol_dengue_data~x[,1]+y[,1]+x[,2]+y[,2]
                  +y[,3]+y[,4]+x[,5])
#y[,2]
dengue.fsha.6<-lm(sol_dengue_data~x[,1]+y[,1]+x[,2]+y[,3]
                  +y[,4]+x[,5])
#y[,3]
dengue.fsha.7<-lm(sol_dengue_data~x[,1]+y[,1]+x[,2]+y[,4]
                  +x[,5])
#x[,5]
dengue.fsha.8<-lm(sol_dengue_data~x[,1]+y[,1]+x[,2]+y[,4])
#y[,4]
dengue.fsha.9<-lm(sol_dengue_data~x[,1]+y[,1]+x[,2])
#x[,2]
dengue.fsha.10<-lm(sol_dengue_data~x[,1]+y[,1])
#x[,1]
dengue.fsha.11<-lm(sol_dengue_data~y[,1])

#plot of AIC
AIC<-AIC(dengue.fsha.1,dengue.fsha.2,dengue.fsha.3,dengue.fsha.4,
    dengue.fsha.5,dengue.fsha.6,dengue.fsha.7,dengue.fsha.8,
    dengue.fsha.9,dengue.fsha.10,dengue.fsha.11)$AIC
plot(Models,AIC)

#plot of model and observed values
fv<-dengue.fsha.7$fitted.values
plot(sol_dengue_data,type="l")
lines(1:40,fv,col="red")

summary(dengue.fsha.1)

# plot of complete data

soldeng_complete_data_ts<-ts(soldeng_complete_data$X1,
                             frequency = 12, start = c(2010,1))

plot(soldeng_complete_data_ts,ylab="No of Cases")

#plot of RSE
RSE<-c(summary(dengue.fsha.1)$sigma, summary(dengue.fsha.2)$sigma,
     summary(dengue.fsha.3)$sigma, summary(dengue.fsha.4)$sigma,
     summary(dengue.fsha.5)$sigma, summary(dengue.fsha.6)$sigma,
     summary(dengue.fsha.7)$sigma, summary(dengue.fsha.8)$sigma,
     summary(dengue.fsha.9)$sigma, summary(dengue.fsha.10)$sigma,
     summary(dengue.fsha.11)$sigma)
Models<-c(1:11)

plot(Models,RSE)

#rsq

rsq<-c(summary(dengue.fsha.1)$r.squared, summary(dengue.fsha.2)$r.squared,
       summary(dengue.fsha.3)$r.squared, summary(dengue.fsha.4)$r.squared,
       summary(dengue.fsha.5)$r.squared, summary(dengue.fsha.6)$r.squared,
       summary(dengue.fsha.7)$r.squared, summary(dengue.fsha.8)$r.squared,
       summary(dengue.fsha.9)$r.squared, summary(dengue.fsha.10)$r.squared,
       summary(dengue.fsha.11)$r.squared)

# adj.rsq
adj.rsq<-c(summary(dengue.fsha.1)$adj.r.squared, summary(dengue.fsha.2)$adj.r.squared,
           summary(dengue.fsha.3)$adj.r.squared, summary(dengue.fsha.4)$adj.r.squared,
           summary(dengue.fsha.5)$adj.r.squared, summary(dengue.fsha.6)$adj.r.squared,
           summary(dengue.fsha.7)$adj.r.squared, summary(dengue.fsha.8)$adj.r.squared,
           summary(dengue.fsha.9)$adj.r.squared, summary(dengue.fsha.10)$adj.r.squared,
           summary(dengue.fsha.11)$adj.r.squared)

plot(Models,rsq,ylab="Values",type = "n",ylim = c(0.15,0.41))
lines(Models,rsq,type = "p")
lines(Models,adj.rsq,type="p",col="red")
legend(x="bottomright",legend=c("R Squared","Adjusted R Squared"),col=c("black","red"),
       pch=c(21,21),bty="o",pt.cex=1.2,cex=.8)

#### FSHA on complete data

# imported data to ts data
soldeng_complete_data_ts<-ts(soldeng_complete_data$X1,
                             frequency = 12, start = c(2010,1))

# plot of complete data
plot(soldeng_complete_data_ts,ylab="No of Cases")

soldeng_complete_data<-as.numeric(soldeng_complete_data_ts)

t<-1:108
x=y=matrix(nr=length(t),nc=11)
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
x[,6]<-cos(pi*t)

dengue.fsha.1<-lm(soldeng_complete_data~x[,1]+y[,1]+x[,2]+y[,2]
                  +x[,3]+y[,3]+x[,4]+y[,4]+x[,5]+y[,5]+x[,6])
#x[,3]
dengue.fsha.2<-lm(soldeng_complete_data~x[,1]+y[,1]+x[,2]+y[,2]
                  +y[,3]+x[,4]+y[,4]+x[,5]+y[,5]+x[,6])
#x[,4]
dengue.fsha.3<-lm(soldeng_complete_data~x[,1]+y[,1]+x[,2]+y[,2]
                  +y[,3]+y[,4]+x[,5]+y[,5]+x[,6])
#x[,6]
dengue.fsha.4<-lm(soldeng_complete_data~x[,1]+y[,1]+x[,2]+y[,2]
                  +y[,3]+y[,4]+x[,5]+y[,5])
#y[,2]
dengue.fsha.5<-lm(soldeng_complete_data~x[,1]+y[,1]+x[,2]+y[,3]
                  +y[,4]+x[,5]+y[,5])
#y[,4]
dengue.fsha.6<-lm(soldeng_complete_data~x[,1]+y[,1]+x[,2]+y[,3]
                  +x[,5]+y[,5])
#y[,3]
dengue.fsha.7<-lm(soldeng_complete_data~x[,1]+y[,1]+x[,2]+x[,5]
                  +y[,5])
#x[,5]
dengue.fsha.8<-lm(soldeng_complete_data~x[,1]+y[,1]+x[,2]+y[,5])

#y[,5]
dengue.fsha.9<-lm(soldeng_complete_data~x[,1]+y[,1]+x[,2])

#x[,2]
dengue.fsha.10<-lm(soldeng_complete_data~x[,1]+y[,1])

#x[,1]
dengue.fsha.11<-lm(soldeng_complete_data~y[,1])

summary(dengue.fsha.11)



#plot of AIC
Models<-c(1:11)

AIC<-AIC(dengue.fsha.1,dengue.fsha.2,dengue.fsha.3,dengue.fsha.4,
         dengue.fsha.5,dengue.fsha.6,dengue.fsha.7,dengue.fsha.8,
         dengue.fsha.9,dengue.fsha.10,dengue.fsha.11)$AIC
plot(Models,AIC)





#plot of RSE
RSE<-c(summary(dengue.fsha.1)$sigma, summary(dengue.fsha.2)$sigma,
       summary(dengue.fsha.3)$sigma, summary(dengue.fsha.4)$sigma,
       summary(dengue.fsha.5)$sigma, summary(dengue.fsha.6)$sigma,
       summary(dengue.fsha.7)$sigma, summary(dengue.fsha.8)$sigma,
       summary(dengue.fsha.9)$sigma, summary(dengue.fsha.10)$sigma,
       summary(dengue.fsha.11)$sigma)

plot(Models,RSE)

#rsq

rsq<-c(summary(dengue.fsha.1)$r.squared, summary(dengue.fsha.2)$r.squared,
       summary(dengue.fsha.3)$r.squared, summary(dengue.fsha.4)$r.squared,
       summary(dengue.fsha.5)$r.squared, summary(dengue.fsha.6)$r.squared,
       summary(dengue.fsha.7)$r.squared, summary(dengue.fsha.8)$r.squared,
       summary(dengue.fsha.9)$r.squared, summary(dengue.fsha.10)$r.squared,
       summary(dengue.fsha.11)$r.squared)

# adj.rsq
adj.rsq<-c(summary(dengue.fsha.1)$adj.r.squared, summary(dengue.fsha.2)$adj.r.squared,
           summary(dengue.fsha.3)$adj.r.squared, summary(dengue.fsha.4)$adj.r.squared,
           summary(dengue.fsha.5)$adj.r.squared, summary(dengue.fsha.6)$adj.r.squared,
           summary(dengue.fsha.7)$adj.r.squared, summary(dengue.fsha.8)$adj.r.squared,
           summary(dengue.fsha.9)$adj.r.squared, summary(dengue.fsha.10)$adj.r.squared,
           summary(dengue.fsha.11)$adj.r.squared)

plot(Models,rsq,ylab="Values",type = "n",ylim = c(0.03,0.14))
lines(Models,rsq,type = "p")
lines(Models,adj.rsq,type="p",col="red")
legend(x="bottomright",legend=c("R Squared","Adjusted R Squared"),col=c("black","red"),
       pch=c(21,21),bty="n",pt.cex=1.2,cex=.8)

#plot of model and observed values
fv<-dengue.fsha.11$fitted.values
fv_ts<-ts(fv,start = c(2010,1),frequency = 12)
ts.plot(soldeng_complete_data_ts,fv_ts,lty=c(1,3),ylab = "Number of Cases")
legend(x="topright",legend=c("Observed","Model"),lty = c(1,3),pt.cex=2,cex=.8,bty = "n")
