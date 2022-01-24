#ADMM algorithm

lasso_threh=function(z,lambda1)
{
max(0,(1-lambda1/(theta*sqrt(sum(z^2)))))
}


lasso_objective=function(z,lambda1)
{
norm_2=sqrt(sum(z^2))
value=lambda1*norm_2
return(value)
}


SCAD_objective=function(z,lambda1,gamma_scad)
{
norm_2=sqrt(sum(z^2))
threh_1=lambda1
threh_2=gamma_scad*lambda1
if(norm_2<=threh_1)
{
value=lambda1*norm_2
}
if(norm_2>threh_1&norm_2<=threh_2)
{
value=(2*gamma_scad*lambda1*norm_2-norm_2^2-lambda1^2)/(2*(gamma_scad-1))
}
if(norm_2>threh_2)
{
value=((lambda1^2)*(gamma_scad+1))/2
}
return(value)
}



SCAD_threh=function(z,lambda1,gamma_scad)
{
norm_2=sqrt(sum(z^2))
threh_1=(2*lambda1)/theta
threh_2=(gamma_scad*lambda1)/theta
if(norm_2<=threh_1)
{
value=max(0,(1-lambda1/(theta*norm_2)))
}
if(norm_2>threh_1&norm_2<=threh_2)
{
temp=1-(lambda1*gamma_scad)/((gamma_scad-1)*theta*norm_2)
value=((gamma_scad-1)/(gamma_scad-2))*max(0,temp)
}
if(norm_2>threh_2)
{
value=1
}
return(value)
}

Empirical_ADMM=function(temp_matrix,B,D,theta,lambda1_S,lambda1_I,lambda2_S,lambda2_I,tol=10^{-8},maxit=500)
{

ini_eta=rep(0,(2*(L*ne)))
ini_v=rep(0,(2*(L*ne)))

eta=ini_eta
v=ini_v

it=1
criteria=1000

while((criteria)>tol&it<maxit)
{

eta_before=eta
#迭代a
temp_vector=(1/n)*(t(B)%*%z)+EE%*%(theta*eta-v)   
a_vector=solve(temp_matrix,temp_vector,sparse=TRUE)

#迭代eta
a_dif=Del%*%a_vector
delta=a_dif+(theta^-1)*v
delta_mat=array(delta,c(L,(2*ne)))

delta_threh=vector(length=(2*ne))
if(penalty_type=="lasso")
{
delta_threh[1:ne]=apply(delta_mat[,1:ne],2,lasso_threh,lambda1=lambda1_S)
delta_threh[(ne+1):(2*ne)]=apply(delta_mat[,(ne+1):(2*ne)],2,lasso_threh,lambda1=lambda1_I)


a_dif_threh=vector(length=(2*ne))
a_dif_mat=array(a_dif,c(L,(2*ne)))
a_dif_threh[1:ne]=apply(a_dif_mat[,1:ne],2,lasso_objective,lambda1=lambda1_S)
a_dif_threh[(ne+1):(2*ne)]=apply(a_dif_mat[,(ne+1):(2*ne)],2,lasso_objective,lambda1=lambda1_I)

res=z-B%*%a_vector
objective_value=(1/n)*0.5*sum(res^2)+sum(a_dif_threh)+as.numeric(0.5*(c(a_vector)%*%D%*%a_vector))

if(track=="True")
{
print(1000*c(((1/n)*0.5*sum(res^2)),sum(a_dif_threh),as.numeric(0.5*(c(a_vector)%*%D%*%a_vector))))
}
}
if(penalty_type=="scad")
{
delta_threh[1:ne]=apply(delta_mat[,1:ne],2,SCAD_threh,lambda1=lambda1_S,gamma_scad=3.7)
delta_threh[(ne+1):(2*ne)]=apply(delta_mat[,(ne+1):(2*ne)],2,SCAD_threh,lambda1=lambda1_I,gamma_scad=3.7)


a_dif_threh=vector(length=(2*ne))
a_dif_mat=array(a_dif,c(L,(2*ne)))
a_dif_threh[1:ne]=apply(a_dif_mat[,1:ne],2,SCAD_objective,lambda1=lambda1_S,gamma_scad=3.7)
a_dif_threh[(ne+1):(2*ne)]=apply(a_dif_mat[,(ne+1):(2*ne)],2,SCAD_objective,lambda1=lambda1_I,gamma_scad=3.7)

res=z-B%*%a_vector
objective_value=(1/n)*0.5*sum(res^2)+sum(a_dif_threh)+as.numeric(0.5*(c(a_vector)%*%D%*%a_vector))

if(track=="True")
{
print(1000*c(((1/n)*0.5*sum(res^2)),sum(a_dif_threh),as.numeric(0.5*(c(a_vector)%*%D%*%a_vector))))
}
}

eta=sweep(delta_mat,2,delta_threh,'*')
eta=as.vector(eta)

#迭代v
v=v+theta*(a_dif-eta)

criteria=sum((a_dif-eta)^2)

criteria1=sum((EE%*%(theta*(eta-eta_before)))^2)

it=it+1

if(track=="True")
{
print(c(it,criteria,criteria1,(1000*objective_value)))
}
}

print(c(it,criteria,criteria1,(1000*objective_value)))

return(list(regression_estimate=a_vector,group_indicator=eta,lagrange=v))
}