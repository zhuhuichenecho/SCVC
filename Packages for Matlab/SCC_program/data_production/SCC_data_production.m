function[x,y]=SCC_simulation_spatial_data_prescribe(beta,epsilon_mag,sim_num,lon,lat,phi,seed0)

% This funcution is used to generate simulation data for the case with
% clustered or smoothly varying coefficients

% Input:
% beta: [n,p] matrix storing the true regression coefficients with the last
%       column corresponds to the spatially varying intercept
% epsilon_mag: standard deviation (sigma) for iid random noise
% sim_num: the number of simulations
% lon: [n,1] vector containing the longitudes of locations
% lat: [n,1] vector containing the latitudes of locations
% phi: range parameter for the covariates
% x: [sim_num,n,p] arrary storing the covariates
% y: [sim_num,n] matrix storing the response variable
% seed0: the seed used to generate random numbers

[n,p]=size(beta);
y=nan(sim_num,n);
x=nan(sim_num,n,p);
[n,p]=size(beta);



 dx=nan(n,n);
 dy=nan(n,n);
 for i=1:n
     for j=1:n
         dx(i,j)=abs(lon(i)-lon(j));
         dy(i,j)=abs(lat(i)-lat(j));
     end
 end

 
    rng(seed0,'twister');
    if p>1
       for k=1:p-1
           if phi==0
               rng(k*sim_num+seed0,'twister');
               x(1,:,k)=randn(n,1);
           else          
               sigma2_mat=exp(-sqrt(dx.^2/phi^2+dy.^2/phi^2));
               rng(k*sim_num+seed0,'twister');
               x(1,:,k)=mvnrnd(zeros(n,1),sigma2_mat);
           end
       end
	   x(1,:,2)=1;
	else
	   x(1,:,p)=ones(1,n); 
    end	
    

 

for t=1:sim_num
    rng(seed0+t-1,'twister');
    epsilon=epsilon_mag*randn(n,1);
    
    if p>1
       for k=1:p-1
           if phi==0
               rng(k*sim_num+seed0+t-1,'twister');
               x(t,:,k)=x(1,:,k);
           else          
               sigma2_mat=exp(-sqrt(dx.^2/phi^2+dy.^2/phi^2));
               rng(k*sim_num+seed0+t-1,'twister');
               x(t,:,k)=x(1,:,k);
           end
       end
       x(t,:,2)=1;
	   y(t,:)=sum(beta.*squeeze(x(t,:,:)),2)'+epsilon';
	else
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	   x(t,:,p)=ones(1,n); 
	   y(t,:)=beta'+epsilon';
    end	
    
end
