
r_basis=function(sample_location,node_location)
# Description: This function generates the low-rank radial basis functions

# INPUT ARGUMENTS
# sample_location:  [n,2] matrix storing the sample locations
# node_location:    [k,2] matrix storing the node locations

# OUTPUT ARGUMENTS
# the basis matrix" [n, k+3]
{
sample_size=dim(sample_location)[1]
node_size=dim(node_location)[1]
expand_s=sample_location[rep(1:sample_size, each = node_size), ]
expand_n=matrix( rep( t(node_location) , sample_size ) , ncol = ncol(node_location) , byrow = TRUE )

expand_difference=expand_s-expand_n

expand_distance=sqrt(expand_difference[,1]^2+expand_difference[,2]^2)

tempor_matrix=cbind(1,sample_location,t(array((expand_distance)*log(expand_distance^expand_distance),c(node_size,sample_size))))

mean_tempor=apply(abs(tempor_matrix),2,mean)


return(t(t(tempor_matrix)/mean_tempor))

}



