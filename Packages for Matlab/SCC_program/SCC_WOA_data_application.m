function[TS_SCC,TS_GWR,BIC_SCC]=SCC_WOA_data_application(latn,depthn,temp,salt,only_intercept)
% Description: Real Data Application -- Estimate the T-S relationship 
% using the SCC and GWR methods

% INPUT ARGUMENTS:
% temp:   temperatue 
% salt:   salinity
% latn:   nondimensional coordinates in the latitudal direction
% depthn: nondimensional coordinates in the vertical direction

% OUTPUT ARGUMENTS
% TS_SCC: estimates of T-S relationship using SCC method
% TS_GWR: estimates of T-S relationship using SCC method

x=temp;
x(:,2)=1;
y=salt;
[n,p]=size(x);

if(only_intercept==0)

% SCC method
[aaa,BIC]=SCC(x,y,latn,depthn,[]);
[~,index]=min(BIC);
BIC_SCC=BIC;
beta_hat_SCC=reshape(aaa(:,index),[n,p]);
TS_SCC=beta_hat_SCC(:,1);

% GWR method
info.dtype='gaussian';
info.bmin=0.01^2;
info.bmax=2^2;
result=gwr(y,x,latn,depthn,info);
beta_hat_GWR=result.beta;
TS_GWR=beta_hat_GWR(:,1);
else
% SCC method
[aaa,BIC]=SCC(ones(n,1),y,latn,depthn,[]);
[~,index]=min(BIC);
beta_hat_SCC=reshape(aaa(:,index),[n,1]);
TS_SCC=beta_hat_SCC;

% GWR method
info.dtype='gaussian';
info.bmin=0.01^2;
info.bmax=2^2;
result=gwr(y,ones(n,1),latn,depthn,info);
beta_hat_GWR=result.beta;
TS_GWR=beta_hat_GWR;
end