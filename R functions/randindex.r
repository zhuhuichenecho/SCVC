#compute the RI criterion

RI=function(spline_coef,group_index_save,num_basis)
{

Ri=vector(length=0)
aaa=group_index_save%*%rep(1,length(group_index_save));
aaa=aaa-t(aaa)
aaa[abs(aaa)!=0]=1;

mask=lower.tri(diag(length(group_index_save)))
aaa=aaa[mask]

spline_Cof=t(array(spline_coef,c(num_basis,(length(spline_coef)/num_basis))))
for(i in 1:ncol(spline_Cof))
{
bbb=cbind(spline_Cof[,i])%*%rep(1,length(group_index_save));
bbb=bbb-t(bbb)
bbb[abs(bbb)>1e-6]=1
bbb[abs(bbb)<1e-6]=0
bbb=bbb[mask]
st=sum(aaa-bbb==0);
Ri[i]=st/length(aaa);
} 
bbb=cbind(apply(abs(spline_Cof),1,mean))%*%rep(1,length(group_index_save));
bbb=bbb-t(bbb)
bbb[abs(bbb)>1e-6]=1
bbb[abs(bbb)<1e-6]=0
bbb=bbb[mask]
st=sum(aaa-bbb==0);
Ri_mean=st/length(aaa);

return(list(RI_min=min(Ri),RI_max=max(Ri),RI_mean=Ri_mean))   
}