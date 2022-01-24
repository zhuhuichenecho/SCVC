
#Implement PSE method.

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

setwd("D:/Program/R functions")
source("low rank radial basis.r")
source("Empirical_ADMM_algorithm.r")
source("group_SVC.r")
source("randindex.r")
source("BIC_ADMM_empirical.r")
source("Spline_estimation.r")

Matlab$startServer()
matlab <- Matlab()
isOpen <- open(matlab)

#derive the MST
evaluate(matlab, "addpath(genpath('D://Program//Packages for Matlab//glmnet_matlab'))")
evaluate(matlab, "addpath(genpath('D://Program//Packages for Matlab//jplv7'))")
evaluate(matlab, "addpath(genpath('D://Program//Packages for Matlab//SCC_program'))")
evaluate(matlab, "load('D:/scholar/new_data/one_cluster_m.mat');")
evaluate(matlab, "[H]=SCC_spanning_tree(lon,lat,1,0.1);")
evaluate(matlab, "[H]=sparse(H);")

tria_method='n'

data <- getVariable(matlab, c("lat", "lon","x","y","beta","group_index_save","group_index_save_intercept"))
H=getVariable(matlab, "H")
H=H$H
sample_location=cbind(data$lon,data$lat)
set.seed(11)
node_number=max(20,min(((dim(sample_location)[1])/4),40))
design_coordinate=cover.design(R=sample_location,nd=node_number,nruns=5)$design
node_location=design_coordinate

data$H=H
data$node_location=node_location
data$sample_location=sample_location

#########################
#specify the lambda for which simulation data, it takes value of 1:100.
#To see a quick result, people can fix the lambda based on the values of several initial sim data set,e.g., sim_index=2. 
sim_index=2

pspline_nopenalty=4
covariates=data$x
beta_Cof=data$beta
sample_location=data$sample_location
node_location=data$node_location
group_index=data$group.index.save
group_index_intercept=data$group.index.save.intercept

H=data$H
Z=data$y
z=Z[sim_index,]

x1=covariates[sim_index,,1]
z_true=beta_Cof[,1]*x1+beta_Cof[,2]
basis_matrix=r_basis(sample_location,node_location)

n=nrow(basis_matrix)
npL=(nrow(basis_matrix)*ncol(basis_matrix))*2
L=ncol(basis_matrix)
pL=ncol(basis_matrix)*2
nL=(nrow(basis_matrix)*ncol(basis_matrix))
np=nrow(basis_matrix)*2

I_B=sparseMatrix(i=1,j=1,x=0,dims=c(n,nL))
S_B=sparseMatrix(i=1,j=1,x=0,dims=c(n,nL))
for(i in 1:n)
{
S_B[i,(L*(i-1)+1):(L*i)]=x1[i]*basis_matrix[i,]
I_B[i,(L*(i-1)+1):(L*i)]=basis_matrix[i,]
}
B=cbind(S_B,I_B)

E=0
EE=sparseMatrix(i=1,j=1,x=0,dims=c(nL,1))
index=matrix(0,nrow=(nrow(H)-1),ncol=2)

ne=(nrow(H)-1)
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
}
Del=bdiag(Del,Del)
EE=EE[,-1]
EE=bdiag(EE,EE)
E=bdiag(E,E)
BB=t(B)%*%B

#Save the data for drawing Figure 2. If MST-equal, save as "plot4.mat". If MST-unequal, save as "plot5.mat".
writeMat('D:/plot4.mat', beta_slope = data$beta[,1],beta_intercept = data$beta[,2], x = data$lon, y = data$lat,index=index,node_location=node_location)

#initial the simplex, given a proper start point.
lambda1_S_seq=c(20,5,20,20,20)/n
lambda1_I_seq=c(20,20,5,20,20)/n

lambda2_S_seq=c(0.05,0.05,0.05,0.4,0.05)/n
lambda2_I_seq=c(0.05,0.05,0.05,0.05,0.4)/n

penalty_type="scad"
with_intercept="True"
theta=1

#indicate whether it is simulation study or empirical study.
Emp="False"

#Details of iterative process is shown if track="True", otherwise track="False".
track="True"

#compute the BIC values on the simplex.
BICC=BIC_TPS(lambda2_S_seq,lambda2_I_seq,B,BB,E,Emp)

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

#when iterative process ends, take lambda_grid[,1] as the selected lambda value, whose BIC value is smallest.

lambda2_S= lambda_grid[1,1]
lambda2_I= lambda_grid[2,1]

BICC=BIC_TPS(lambda2_S,lambda2_I,B,BB,E,"False")
BICC

#doing simulation 100 times.
lambda2_S= lambda_grid[1,1]
lambda2_I= lambda_grid[2,1]

track="False"

D=sparseMatrix(i=seq(L),j=seq(L),x=1,dims=c(L,L))
diag(D)[1:(pspline_nopenalty-1)]=0
Identity_matrix=sparseMatrix(i=seq(np),j=seq(np),x=1,dims=c(np,np))
diag(Identity_matrix)[1:n]=lambda2_S
diag(Identity_matrix)[(n+1):np]=lambda2_I
D=kronecker(Identity_matrix,D)

MSE_tps_beta=vector(length=0)
MSE_tps_intercept=vector(length=0)

sim_times=100
for(simm in 1:sim_times)
{

z=Z[simm,]
x1=covariates[simm,,1]


BICC=BIC_TPS(lambda2_S,lambda2_I,B,BB,E,"False")
MSE_tps_beta[simm]=BICC$MSE_svc_beta
MSE_tps_intercept[simm]=BICC$MSE_svc_intercept

print(c(MSE_tps_beta[simm],MSE_tps_intercept[simm]))
print(c("simulated times",simm))
}

#compute the results for PSE method
mean(MSE_tps_beta)
aa=MSE_tps_beta
MSE_var=sqrt((mean((aa)^2)-(mean(aa))^2)/length(aa))
MSE_var

mean(MSE_tps_intercept)
aa=MSE_tps_intercept
MSE_var=sqrt((mean((aa)^2)-(mean(aa))^2)/length(aa))
MSE_var	


#compute the results for GWR and SCC.
SCC_mse_I=vector(length=0)
SCC_RI_I=vector(length=0)
SCC_IC_I=vector(length=0)
SCC_mse_S=vector(length=0)
SCC_RI_S=vector(length=0)
SCC_IC_S=vector(length=0)
GWR_mse_I=vector(length=0)
GWR_mse_S=vector(length=0)

#load the results from "simdata-MSTequal.m" or "simdata-MSTunequal.m"
beta_hat_SCC_o=readMat('D:/beta_hat_SCC.mat')
beta_hat_GWR_o=readMat('D:/beta_hat_GWR.mat')

for(i in 1:sim_times)
{

beta_hat_SCC=beta_hat_SCC_o$ beta.hat.SCC[i,,1]
SCC_IC_S[i]=length(unique(beta_hat_SCC))
SCC_mse_S[i]=mean((beta_hat_SCC-data$beta[,1])^2)
SCC_RI_S[i]=RI(beta_hat_SCC,group_index,1)$RI_min

beta_hat_SCC=beta_hat_SCC_o$ beta.hat.SCC[i,,2]
SCC_IC_I[i]=length(unique(beta_hat_SCC))
SCC_mse_I[i]=mean((beta_hat_SCC-data$beta[,2])^2)
SCC_RI_I[i]=RI(beta_hat_SCC,group_index,1)$RI_min

beta_hat_GWR=beta_hat_GWR_o$ beta.hat.GWR[i,,1]
GWR_mse_S[i]=mean((beta_hat_GWR-data$beta[,1])^2)

beta_hat_GWR=beta_hat_GWR_o$ beta.hat.GWR[i,,2]
GWR_mse_I[i]=mean((beta_hat_GWR-data$beta[,2])^2)



}

mean(GWR_mse_S)
aa=GWR_mse_S
sqrt((mean((aa)^2)-(mean(aa))^2)/length(aa))

mean(GWR_mse_I)
aa=GWR_mse_I
sqrt((mean((aa)^2)-(mean(aa))^2)/length(aa))



mean(SCC_mse_S)
aa=SCC_mse_S
sqrt((mean((aa)^2)-(mean(aa))^2)/length(aa))

mean(SCC_mse_I)
aa=SCC_mse_I
sqrt((mean((aa)^2)-(mean(aa))^2)/length(aa))

mean(SCC_RI_S)
aa=SCC_RI_S
sqrt((mean((aa)^2)-(mean(aa))^2)/length(aa))

mean(SCC_RI_I)
aa=SCC_RI_I
sqrt((mean((aa)^2)-(mean(aa))^2)/length(aa))

mean(SCC_IC_S)
aa=SCC_IC_S
sqrt((mean((aa)^2)-(mean(aa))^2)/length(aa))

mean(SCC_IC_I)
aa=SCC_IC_I
sqrt((mean((aa)^2)-(mean(aa))^2)/length(aa))







