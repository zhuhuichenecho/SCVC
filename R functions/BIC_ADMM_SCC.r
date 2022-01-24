#compute the BIC value based on SCC.

BIC_ADMM_SCC=function(lambda1_S_seq,lambda1_I_seq,B,BB,E,Emp,maxit)
{
tt=1
BIC=vector(length=0)
AIC=vector(length=0)
random_index_S=vector(length=0)
random_index_I=vector(length=0)
lambd1_S=vector(length=0)
lambd1_I=vector(length=0)
MSE=vector(length=0)
MSE_var=vector(length=0)
group_num_vector_I=vector(length=0)
group_num_vector_S=vector(length=0)
MSE_scc_beta=vector(length=0)
MSE_scc_intercept=vector(length=0)

for(lam1 in seq_along(lambda1_S_seq))
{

lambda1_S=lambda1_S_seq[lam1]
lambda1_I=lambda1_I_seq[lam1]

temp_matrix_in=(1/n)*BB+theta*E

delta_I=matrix(0,nrow=((nrow(H)-1)*L),ncol=1)
delta_S=matrix(0,nrow=((nrow(H)-1)*L),ncol=1)

scc_results=SCC_ADMM(temp_matrix_in,B,theta=theta,lambda1_S=lambda1_S,lambda1_I=lambda1_I,tol=10^{-20},maxit=maxit)

rm(temp_matrix_in)
gc()

group_indicator=scc_results$group_indicator
group_indicator_S=group_indicator[1:(ne*L)]
group_indicator_I=group_indicator[(ne*L+1):(ne*L*2)]

groups_S=t(array(group_indicator_S,c(L,ne)))
groups_S=apply(abs(groups_S),1,sum)

groups_I=t(array(group_indicator_I,c(L,ne)))
groups_I=apply(abs(groups_I),1,sum)

spline_coef=scc_results$regression_estimate

basis_matrix=matrix(0,nrow=length(z),ncol=1)
group_div_S=group(index,groups_S,basis_matrix)
group_div_I=group(index,groups_I,basis_matrix)

group_num_S=length(group_div_S)
group_num_I=length(group_div_I)

print(c("group number", group_num_S, group_num_I))

spline_num=ncol(basis_matrix)

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

if(Emp=="False")
{
random_index_S[tt]=RI(beta_estimate_S,group_index,1)$RI_min
random_index_I[tt]=RI(beta_estimate_I,group_index_intercept,1)$RI_min
}

Z_fit=S_B%*%spline_coef[1:nL]+I_B%*%spline_coef[(nL+1):npL]


# compute the degree of freedom
SB=sparseMatrix(i=1,j=1,x=0,dims=c(n,(L*group_num_S)))

for(i in 1:group_num_S)
{
SB[group_div_S[[i]],(L*(i-1)+1):(L*i)]=x1[group_div_S[[i]]]
}

IB=sparseMatrix(i=1,j=1,x=0,dims=c(n,(L*group_num_I)))

for(i in 1:group_num_I)
{
IB[group_div_I[[i]],(L*(i-1)+1):(L*i)]=1
}


SB_IB=cbind(SB,IB)

SBD=t(SB_IB)%*%SB_IB

SS_II=t(SB_IB)%*%SB_IB

SH=try(solve(SBD,(SS_II),sparse=TRUE),silent=TRUE)
if('try-error' %in% class(SH))
{
 dof=10^5
 }else{
 dof=sum(diag(SH))
 }
 
 dof=nrow(SBD)


print(c("degree of freedom",dof))


BIC[tt]=length(z)*log(mean((z-Z_fit)^2))+log(length(z))*dof
AIC[tt]=length(z)*log(mean((z-Z_fit)^2))+2*dof


lambd1_I[tt]=lambda1_I
lambd1_S[tt]=lambda1_S
group_num_vector_S[tt]=group_num_S
group_num_vector_I[tt]=group_num_I

if(Emp=="False")
{
MSE[tt]=mean((z_true-Z_fit)^2)
MSE_var[tt]=sqrt((mean((z_true-Z_fit)^4)-(mean((z_true-Z_fit)^2))^2)/length(z_true))
MSE_scc_beta[tt]=mean((data$beta[,1]-I_B%*%spline_coef[1:nL])^2)
MSE_scc_intercept[tt]=mean((data$beta[,2]-I_B%*%spline_coef[(nL+1):npL])^2)
}


print(c("iteration",tt))
tt=tt+1


rm(SB_IB)
rm(SBD)
rm(SS_II)
rm(SH)
gc()

}

if(Emp=="False")
{
return(list(MSE_scc_beta=MSE_scc_beta,MSE_scc_intercept=MSE_scc_intercept,random_index_S=random_index_S,random_index_I=random_index_I,group_num_S=group_num_vector_S,group_num_I=group_num_vector_I,MSE=MSE,MSE_var=MSE_var,AIC=AIC,BIC=BIC,lambda1_S=lambd1_S,lambda1_I=lambd1_I))

}else{
return(list(group_num_S=group_num_vector_S,group_num_I=group_num_vector_I,AIC=AIC,BIC=BIC,lambda1_S=lambd1_S,lambda1_I=lambd1_I))
}

}