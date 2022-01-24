#This is for empirical analysis.

#SCC* method.

#The codes is run in R 3.6.0, windows system.

rm(list=ls(all=TRUE))
library(RTriangle)
library(R.matlab)
library(caret)
library(fields)
library(rgl)
library(plot3D)
library(Matrix)
library(geoR)
library(gstat)
library(Rlinsolve)

setwd("D:/Program/empirical program")
data=readRDS("data_for_sever_south.rds")

setwd("D:/Program/R functions")
source("Empirical_ADMM_algorithm.r")
source("group_SVC.r")
source("randindex.r")
source("BIC_ADMM_empirical.r")
source("BIC_ADMM_SCC.r")
source("SCC_ADMM_algorithm.r")

H=data$H
latn=data$latn
depthn=data$depthn
z=data$z
x1=data$x1
sample_location=cbind(latn,depthn)
x=latn
y=depthn


n=length(z)
npL=length(z)*2
L=1
nL=length(z)
np=length(z)*2
ne=(nrow(H)-1)

I_B=sparseMatrix(i=1,j=1,x=0,dims=c(n,nL))
S_B=sparseMatrix(i=1,j=1,x=0,dims=c(n,nL))
for(i in 1:n)
{
S_B[i,(L*(i-1)+1):(L*i)]=x1[i]
I_B[i,(L*(i-1)+1):(L*i)]=1
print(i)
}
B=cbind(S_B,I_B)
BB=t(B)%*%B

E=0
EE=sparseMatrix(i=1,j=1,x=0,dims=c(nL,1))
index=matrix(0,nrow=(nrow(H)-1),ncol=2)

Del=sparseMatrix(i=1,j=1,x=0,dims=c((ne*L),(n*L)))

for(i in 1:ne)
{
index[i,]=which(abs(H[i,])==1)

Del[((i-1)*L+1):(i*L),((index[i,1]-1)*L+1):(index[i,1]*L)]=1
Del[((i-1)*L+1):(i*L),((index[i,2]-1)*L+1):(index[i,2]*L)]=-1

B1=sparseMatrix(i=1,j=1,x=0,dims=c(L,nL))
B2=sparseMatrix(i=1,j=1,x=0,dims=c(L,nL))

B1[,((index[i,1]-1)*L+1):(index[i,1]*L)]=1
B2[,((index[i,2]-1)*L+1):(index[i,2]*L)]=1

temp=B1-B2

E=t(temp)%*%(temp)+E

EE=cbind(EE,t(temp))
print(i)
}
Del=bdiag(Del,Del)
EE=EE[,-1]
EE=bdiag(EE,EE)
E=bdiag(E,E)

with_intercept="True"
penalty_type="scad"
track="True"
Emp="True"

ini_type="lasso"
beta_hat=readMat('D:/beta_hat_SCC.mat')
beta_h=beta_hat$beta.hat
a_h=c(beta_h[,1],beta_h[,2])


lambda1_S_seq=c(2,1,2)/n
lambda1_I_seq=c(2,2,1)/n

theta=1

#calculate the BIC values over the simplex
BICC=BIC_ADMM_SCC(lambda1_S_seq,lambda1_I_seq,B,BB,E,"True",500)


#This is to implement the Nelder-Mead algorithm if the BIC values of simplex have been obtained; 

if(FALSE)
{
track="False"
num_lab=2

nla=length(lambda1_S_seq)

lambda_grid=matrix(0,nrow=num_lab,ncol=nla)

for(i in 1:nla)
{
lambda_grid[,i]=c(lambda1_S_seq[i],lambda1_I_seq[i])
}

BIC_value=BICC$BIC

it=1
#the grid points of simplex and corresponding BIC value is usually quite close wthin 30 iterations, 
#if not, we can set larger maxit.
maxit=30


#turn to the invalid negative values of lambda into very small positive values.
if(Emp=="True")
{
thred=0
dd=0.01/n
positive=c(dd,dd)
}else{
thred=0.005
positive=c(0.005,0.005)}

while(it<maxit)
{
print(c("teration time",it))
BIC_order=order(BIC_value)
BIC_value=BIC_value[BIC_order]
lambda_grid=lambda_grid[,BIC_order]

print(rbind(lambda_grid,BIC_value))

temp_grid=apply(lambda_grid[,1:(nla-1)],1,mean)

#compute reflected point
alpha=1
reflect_grid=temp_grid+alpha*(temp_grid-lambda_grid[,nla])
reflect_grid[reflect_grid<=thred]=positive[reflect_grid<=thred]

if(Emp=="True")
{
if(reflect_grid[2]<dd|reflect_grid[1]<dd)
{
reflect_grid[1]=dd
reflect_grid[2]=dd
}
}

BICC_reflect=BIC_ADMM_SCC(reflect_grid[1],reflect_grid[2],B,BB,E,Emp,500)$BIC

if(BICC_reflect<BIC_value[(nla-1)]&BICC_reflect>=BIC_value[1])
{
lambda_grid[,nla]=reflect_grid
BIC_value[nla]=BICC_reflect
}

#expansion
Gamma=2
if(BICC_reflect<BIC_value[1])
{
expand_grid=temp_grid+Gamma*(reflect_grid-temp_grid)
expand_grid[expand_grid<=thred]=positive[expand_grid<=thred]

if(Emp=="True")
{
if(expand_grid[2]<dd|expand_grid[1]<dd)
{
expand_grid[1]=dd
expand_grid[2]=dd
}
}

BICC_expand=BIC_ADMM_SCC(expand_grid[1],expand_grid[2],B,BB,E,Emp,500)$BIC

if(BICC_expand<BICC_reflect)
{
lambda_grid[,nla]=expand_grid
BIC_value[nla]=BICC_expand
}else{
lambda_grid[,nla]=reflect_grid
BIC_value[nla]=BICC_reflect
}
}

#contraction
if(BICC_reflect>=BIC_value[(nla-1)])
{
Rho=0.5
contract_grid=temp_grid+Rho*(lambda_grid[,nla]-temp_grid)
contract_grid[contract_grid<=thred]=positive[contract_grid<=thred]

if(Emp=="True")
{
if(contract_grid[2]<dd|contract_grid[1]<dd)
{
contract_grid[1]=dd
contract_grid[2]=dd
}
}

BICC_contract=BIC_ADMM_SCC(contract_grid[1],contract_grid[2],B,BB,E,Emp,500)$BIC
if(BICC_contract<BIC_value[nla])
{
lambda_grid[,nla]=contract_grid
BIC_value[nla]=BICC_contract
}else{
Sigma=0.5
#shrink
for(s in 2:nla)
{

lambda_grid[,s]=lambda_grid[,1]+Sigma*(lambda_grid[,s]-lambda_grid[,1])
BIC_value[s]=BIC_ADMM_SCC(lambda_grid[1,s],lambda_grid[2,s],B,BB,E,Emp,500)$BIC
}

}

}
it=it+1
}

lambda1_S= lambda_grid[1]
lambda1_I= lambda_grid[2]
}



#turn the above if(FALSE) to if(TRUE) if running the Nelder-Mead algorithm.
#under if(FALSE), we also provide the final selected lambda values by the above Nelder-Mead algorithm.
lambda1_S=6.813272e-04 
lambda1_I=3.172311e-04

temp_matrix=(1/n)*BB+theta*E

svc_results=SCC_ADMM(temp_matrix,B,theta=theta,lambda1_S=lambda1_S,lambda1_I=lambda1_I,tol=10^{-20},maxit=500)

rm(temp_matrix)
gc()

spline_coef=svc_results$regression_estimate
beta_svc_S=I_B%*%spline_coef[1:nL]
beta_svc_I=I_B%*%spline_coef[(nL+1):npL]

group_indicator=svc_results$group_indicator
group_indicator_S=group_indicator[1:(ne*L)]
group_indicator_I=group_indicator[(ne*L+1):(ne*L*2)]

groups_S=t(array(group_indicator_S,c(L,ne)))
groups_S=apply(abs(groups_S),1,sum)

groups_I=t(array(group_indicator_I,c(L,ne)))
groups_I=apply(abs(groups_I),1,sum)

basis_matrix=matrix(0,nrow=length(z),ncol=1)
group_div_S=group(index,groups_S,basis_matrix)
group_div_I=group(index,groups_I,basis_matrix)

length(group_div_S)
length(group_div_I)

beta_estimate_S=vector(length=n)
beta_estimate_I=vector(length=n)

for(gindex in 1:length(group_div_S))
{
beta_estimate_S[group_div_S[[gindex]]]=gindex
}


for(gindex in 1:length(group_div_I))
{
beta_estimate_I[group_div_I[[gindex]]]=gindex
}

writeMat('D:/scholar/SCC_program/new_data/scc_results_scad2.mat', beta_svc=beta_svc_S,beta_svc_I=beta_svc_I,beta_estimate_I=beta_estimate_I,beta_estimate_S=beta_estimate_S)









