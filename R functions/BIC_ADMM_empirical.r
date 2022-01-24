#compute the BIC value based on ADMM algorithm

BIC_ADMM_empirical=function(lambda1_S_seq,lambda2_S_seq,lambda1_I_seq,lambda2_I_seq,B,BB,E,Emp,maxit)
{
tt=1
BIC=vector(length=0)
AIC=vector(length=0)
random_index_S=vector(length=0)
random_index_I=vector(length=0)
lambd1_S=vector(length=0)
lambd2_S=vector(length=0)
lambd1_I=vector(length=0)
lambd2_I=vector(length=0)
MSE=vector(length=0)
MSE_var=vector(length=0)
group_num_vector_I=vector(length=0)
group_num_vector_S=vector(length=0)
MSE_svc_beta=vector(length=0)
MSE_svc_intercept=vector(length=0)

for(lam1 in seq_along(lambda1_S_seq))
{

lambda1_S=lambda1_S_seq[lam1]
lambda2_S=lambda2_S_seq[lam1]

lambda1_I=lambda1_I_seq[lam1]
lambda2_I=lambda2_I_seq[lam1]

D=sparseMatrix(i=seq(L),j=seq(L),x=1,dims=c(L,L))
diag(D)[1:(pspline_nopenalty-1)]=0
Identity_matrix=sparseMatrix(i=seq(np),j=seq(np),x=1,dims=c(np,np))
diag(Identity_matrix)[1:n]=lambda2_S
diag(Identity_matrix)[(n+1):np]=lambda2_I
D=kronecker(Identity_matrix,D)

if(with_intercept=="False")
{diag(D)[(nL+1):npL]=lambda2_I
print("no intercept considered")
}

temp_matrix_in=(1/n)*BB+theta*E+D

delta_I=matrix(0,nrow=((nrow(H)-1)*L),ncol=1)
delta_S=matrix(0,nrow=((nrow(H)-1)*L),ncol=1)


svc_results=Empirical_ADMM(temp_matrix_in,B,D,theta=theta,lambda1_S=lambda1_S,lambda1_I=lambda1_I,lambda2_S=lambda2_S,lambda2_I=lambda2_I,tol=10^{-8},maxit=maxit)

rm(temp_matrix_in)
gc()

group_indicator=svc_results$group_indicator
group_indicator_S=group_indicator[1:(ne*L)]
group_indicator_I=group_indicator[(ne*L+1):(ne*L*2)]

groups_S=t(array(group_indicator_S,c(L,ne)))
groups_S=apply(abs(groups_S),1,sum)

groups_I=t(array(group_indicator_I,c(L,ne)))
groups_I=apply(abs(groups_I),1,sum)

spline_coef=svc_results$regression_estimate

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
SD=sparseMatrix(i=1,j=1,x=0,dims=c((L*group_num_S),(L*group_num_S)))

for(i in 1:group_num_S)
{
SB[group_div_S[[i]],(L*(i-1)+1):(L*i)]=sweep(as.matrix(basis_matrix[group_div_S[[i]],]),1,x1[group_div_S[[i]]],'*')
diag(SD)[(L*(i-1)+pspline_nopenalty):(L*i)]=length(group_div_S[[i]])
}

if(Emp=="False")
{
diag(SD)[diag(SD)==0]=10^{-11}
}else{
diag(SD)[diag(SD)==0]=10^{-6}
}

IB=sparseMatrix(i=1,j=1,x=0,dims=c(n,(L*group_num_I)))
ID=sparseMatrix(i=1,j=1,x=0,dims=c((L*group_num_I),(L*group_num_I)))

for(i in 1:group_num_I)
{
IB[group_div_I[[i]],(L*(i-1)+1):(L*i)]=as.matrix(basis_matrix[group_div_I[[i]],])
diag(ID)[(L*(i-1)+pspline_nopenalty):(L*i)]=length(group_div_I[[i]])
}

if(Emp=="False")
{
diag(ID)[diag(ID)==0]=10^{-11}
}else{
diag(ID)[diag(ID)==0]=10^{-6}
}


if(with_intercept=="False")
{diag(ID)=n
print("no intercept considered")
}



SB_IB=cbind(SB,IB)
SD_ID=bdiag((lambda2_S*SD),(lambda2_I*ID))

SBD=(t(SB_IB)%*%SB_IB+(n*SD_ID))

SS_II=t(SB_IB)%*%SB_IB

SH=try(solve(SBD,(SS_II),sparse=TRUE),silent=TRUE)
if('try-error' %in% class(SH))
{
 dof=10^5
 }else{
 dof=sum(diag(SH))
 }


print(c("degree of freedom",dof))


BIC[tt]=length(z)*log(mean((z-Z_fit)^2))+log(length(z))*dof
AIC[tt]=length(z)*log(mean((z-Z_fit)^2))+2*dof


lambd1_I[tt]=lambda1_I
lambd2_I[tt]=lambda2_I
lambd1_S[tt]=lambda1_S
lambd2_S[tt]=lambda2_S
group_num_vector_S[tt]=group_num_S
group_num_vector_I[tt]=group_num_I

if(Emp=="False")
{
MSE[tt]=mean((z_true-Z_fit)^2)
MSE_var[tt]=sqrt((mean((z_true-Z_fit)^4)-(mean((z_true-Z_fit)^2))^2)/length(z_true))
MSE_svc_beta[tt]=mean((data$beta[,1]-I_B%*%spline_coef[1:nL])^2)
MSE_svc_intercept[tt]=mean((data$beta[,2]-I_B%*%spline_coef[(nL+1):npL])^2)
}


print(c("iteration",tt))
tt=tt+1


rm(SB_IB)
rm(SD_ID)
rm(SBD)
rm(SS_II)
rm(SH)
gc()

}

if(Emp=="False")
{
return(list(MSE_svc_beta=MSE_svc_beta,MSE_svc_intercept=MSE_svc_intercept,random_index_S=random_index_S,random_index_I=random_index_I,group_num_S=group_num_vector_S,group_num_I=group_num_vector_I,MSE=MSE,MSE_var=MSE_var,AIC=AIC,BIC=BIC,lambda1_S=lambd1_S,lambda1_I=lambd1_I,lambda2_S=lambd2_S,lambda2_I=lambd2_I))

}else{
return(list(group_num_S=group_num_vector_S,group_num_I=group_num_vector_I,AIC=AIC,BIC=BIC,lambda1_S=lambd1_S,lambda1_I=lambd1_I,lambda2_S=lambd2_S,lambda2_I=lambd2_I))
}

}