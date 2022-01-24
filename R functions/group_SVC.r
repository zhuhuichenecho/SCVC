#compute the group membership from ADMM.

group=function(index,groups,basis_matrix)
{

connect=which(groups==0)
connect_pairs=index[connect,]

group_value=matrix(0,nrow=(nrow(basis_matrix)),ncol=1)
tt=1

for(i in 1:nrow(connect_pairs))
{

pair_location=connect_pairs[i,]

if((group_value[pair_location[1]]!=0)&(group_value[pair_location[2]]==0))
{
group_value[pair_location[2]]=group_value[pair_location[1]]
}
if((group_value[pair_location[1]]==0)&(group_value[pair_location[2]]!=0))
{
group_value[pair_location[1]]=group_value[pair_location[2]]
}
if((group_value[pair_location[1]]==0)&(group_value[pair_location[2]]==0))
{
group_value[pair_location[1]]=tt
group_value[pair_location[2]]=tt
tt=tt+1
}
if((group_value[pair_location[1]]!=0)&(group_value[pair_location[2]]!=0)&(group_value[pair_location[2]]==group_value[pair_location[1]]))
{
next
}
if((group_value[pair_location[1]]!=0)&(group_value[pair_location[2]]!=0)&(group_value[pair_location[2]]!=group_value[pair_location[1]]))
{
index1=which(group_value==group_value[pair_location[1]])
index2=which(group_value==group_value[pair_location[2]])
group_value[c(index1,index2)]=min(group_value[pair_location[1]],group_value[pair_location[2]])
}
}

group_divide=list()

iso_group=which(group_value==0)
if(length(iso_group)==0)
{
u_group=unique(group_value)
group_number=length(u_group)
}else{
u_group=unique(group_value[-iso_group])
group_number=length(u_group)+length(iso_group)
}
for(i in 1:group_number)
{
if(i<=length(u_group))
{
group_divide[[i]]=which(group_value==u_group[i])
}else{
group_divide[[i]]=iso_group[i-length(u_group)]
}

}


return(group_divide)







}