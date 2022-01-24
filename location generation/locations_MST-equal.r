# This is for generating the MST-equal locations in our simulation.

#The codes is run in R 3.6.0, windows system.

library(R.matlab)
set.seed(20)
x1=vector(length=0)
y1=vector(length=0)
x2=vector(length=0)
y2=vector(length=0)
x3=vector(length=0)
y3=vector(length=0)
x4=vector(length=0)
y4=vector(length=0)
dis1=vector(length=0)
dis_up2=vector(length=0)
dis_down2=vector(length=0)
dis_up3=vector(length=0)
dis_down3=vector(length=0)
dis4=vector(length=0)
tol=0.02
t1=1
t2=1
t3=1
t4=1
tt=0
while(tt<=1000)
{
x=runif(1,0,1)
y=runif(1,0,1)
if(t1<101&y>(0.5+x))
{
dis_temp=abs(0.5+x-y)/sqrt(2)
if(dis_temp<tol)
{next}
x1[t1]=x
y1[t1]=y
dis1[t1]=abs(0.5+x-y)/sqrt(2)
t1=t1+1
}
if(t2<401&y>x&y<(0.5+x))
{
dis_temp=min((abs(0.5+x-y)/sqrt(2)),((abs(x-y)/sqrt(2))))
if(dis_temp<tol)
{next}
x2[t2]=x
y2[t2]=y
dis_up2[t2]=(abs(0.5+x-y)/sqrt(2))
dis_down2[t2]=(abs(x-y)/sqrt(2))
t2=t2+1
}
if(t3<401&y<x&y>(x-0.5))
{
dis_temp=min((abs(x-y-0.5)/sqrt(2)),((abs(x-y)/sqrt(2))))
if(dis_temp<tol)
{next}
x3[t3]=x
y3[t3]=y
dis_up3[t3]=(abs(x-y)/sqrt(2))
dis_down3[t3]=(abs(x-y-0.5)/sqrt(2))
t3=t3+1
}
if(t4<101&y<(x-0.5))
{
dis_temp=abs(x-y-0.5)/sqrt(2)
if(dis_temp<tol)
{next}
x4[t4]=x
y4[t4]=y
dis4[t4]=abs(x-y-0.5)/sqrt(2)
t4=t4+1
}
tt=t1+t2+t3+t4-3
}

plot(c(x1,x2,x3,x4),c(y1,y2,y3,y4),type='p',pch=20)
x4[which.min(x4)]=0.55
x=c(x1,x2,x3,x4)
y=c(y1,y2,y3,y4)

beta_value=matrix(0,nrow=length(x),ncol=3)
beta_value[1:length(x1),]=1
beta_value[(length(x1)+1):(length(x1)+length(x2)),]=-1
beta_value[(length(x1)+length(x2)+1):(length(x1)+length(x2)+length(x3)),]=0.5
beta_value[(length(x1)+length(x2)+length(x3)+1):(length(x1)+length(x2)+length(x3)+length(x4)),]=-0.5

#specify a path for saving results. 
#We have already provided "svc_pattern.mat" in the file "location_generation", 
#following code is not necessary.
writeMat('D:/svc_pattern.mat', lon=cbind(x),lat=cbind(y),beta=beta_value)
