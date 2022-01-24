#This is for empirical analysis.

#PSE method.

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
source("low rank radial basis.r")
source("Empirical_ADMM_algorithm.r")
source("group_SVC.r")
source("randindex.r")
source("BIC_ADMM_empirical.r")
source("Spline_estimation.r")
source("Spline_estimation2.r")

H=data$H
latn=data$latn
depthn=data$depthn
z=data$z
x1=data$x1
sample_location=cbind(latn,depthn)
x=latn
y=depthn

#space-filling design
set.seed(11)
node_number=max(20,min((length(x)/4),40))
design_coordinate=cover.design(R=sample_location,nd=node_number)$design
node_location=design_coordinate
pspline_nopenalty=4
basis_matrix=r_basis(sample_location,node_location)

n=nrow(basis_matrix)
npL=(nrow(basis_matrix)*ncol(basis_matrix))*2
L=ncol(basis_matrix)
pL=L*2
nL=(nrow(basis_matrix)*ncol(basis_matrix))
np=nrow(basis_matrix)*2
ne=(nrow(H)-1)

I_B=sparseMatrix(i=1,j=1,x=0,dims=c(n,nL))
S_B=sparseMatrix(i=1,j=1,x=0,dims=c(n,nL))
for(i in 1:n)
{
S_B[i,(L*(i-1)+1):(L*i)]=x1[i]*basis_matrix[i,]
I_B[i,(L*(i-1)+1):(L*i)]=basis_matrix[i,]
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

diag(Del[((i-1)*L+1):(i*L),((index[i,1]-1)*L+1):(index[i,1]*L)])=1
diag(Del[((i-1)*L+1):(i*L),((index[i,2]-1)*L+1):(index[i,2]*L)])=-1

B1=sparseMatrix(i=1,j=1,x=0,dims=c(L,nL))
B2=sparseMatrix(i=1,j=1,x=0,dims=c(L,nL))

diag(B1[,((index[i,1]-1)*L+1):(index[i,1]*L)])=1
diag(B2[,((index[i,2]-1)*L+1):(index[i,2]*L)])=1

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
track="False"
Emp="True"

lambda1_S_seq=c(30,20,30,30,30)/n
lambda1_I_seq=c(30,30,20,30,30)/n

lambda2_S_seq=c(0.05,0.05,0.05,0.02,0.05)/n
lambda2_I_seq=c(0.05,0.05,0.05,0.05,0.02)/n


theta=1

#calculate the BIC values over the simplex
BICC=BIC_TPS(lambda2_S_seq,lambda2_I_seq,B,BB,E,"True")

#This is to implement the Nelder-Mead algorithm if the BIC values of simplex have been obtained; 

track="False"
num_lab=2

nla=length(lambda1_S_seq)

lambda_grid=matrix(0,nrow=num_lab,ncol=nla)

for(i in 1:nla)
{
lambda_grid[,i]=c(lambda2_S_seq[i],lambda2_I_seq[i])
}

BIC_value=BICC$BIC

it=1
#the grid points of simplex and corresponding BIC value is usually quite close wthin 30 iterations, 
#if not, we can set larger maxit.
maxit=30

#turn to the invalid negative values of lambda into very small positive values.
thred=0
if(Emp=="True")
{
positive=c((0.005/n),(0.005/n))
}else{positive=c(0.000001,0.000001)/n}


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

BICC_reflect=BIC_TPS(reflect_grid[1],reflect_grid[2],B,BB,E,Emp)$BIC

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


BICC_expand=BIC_TPS(expand_grid[1],expand_grid[2],B,BB,E,Emp)$BIC

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

BICC_contract=BIC_TPS(contract_grid[1],contract_grid[2],B,BB,E,Emp)$BIC
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
BIC_value[s]=BIC_TPS(lambda_grid[1,s],lambda_grid[2,s],B,BB,E,Emp)$BIC
}

}

}
it=it+1
}



lambda2_S= lambda_grid[1,1]
lambda2_I= lambda_grid[2,1]

result=BIC_TPS2(lambda2_S,lambda2_I,B,BB,E,"True")

beta_svc_S=result$beta_svc_S
beta_svc_I=result$beta_svc_I


writeMat('D:/scholar/SCC_program/new_data/tps_results.mat', beta_svc=beta_svc_S,beta_svc_I=beta_svc_I)









