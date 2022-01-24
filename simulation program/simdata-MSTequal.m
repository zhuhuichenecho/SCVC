%This is for generating simulation data of MST-equal.  

%Give the results for SCC and GWR, using the program provided by Li and Sang (2019).

%The codes is run in Matlab R2016a, windows system.

clear all 

%load the program for SCC and GWR, they are provided in the file "Packages for Matlab".
%you need to specify the right path.
addpath(genpath('C:\Users\huichenzhu\Dropbox\CUHK\Research\Paper\SCVC\JCGS\Program-SCVC\Program\Packages for Matlab\glmnet_matlab'))
addpath(genpath('C:\Users\huichenzhu\Dropbox\CUHK\Research\Paper\SCVC\JCGS\Program-SCVC\Program\Packages for Matlab\jplv7'))
addpath(genpath('C:\Users\huichenzhu\Dropbox\CUHK\Research\Paper\SCVC\JCGS\Program-SCVC\Program\Packages for Matlab\SCC_program'))

% simulated data
sim_num=100; % simulation numbers
epsilon_mag=0.1; % standard deviation (sigma) for iid random noise
phi=1;%range parameter for the covariates
which_pattern=1;

% load the locations of MST-equal, already saved in "location generation".
load('C:\Users\huichenzhu\Dropbox\CUHK\Research\Paper\SCVC\JCGS\Program-SCVC\Program\location generation\svc_pattern.mat');

% If beta(s) is smooth within each subregion, using "smooth_cluster". If constant, using "constant_cluster".
scenario='smooth_cluster';

% give values of beta on the locations
if strcmp(scenario,'constant_cluster')
   pro=1;
   group_index_save=beta(:,which_pattern);
   group_index_save_intercept=beta(:,which_pattern);
   beta=beta(:,which_pattern).*pro;  
   z_true_smooth=beta;
   unique_z_true=unique(beta);
   for zz=1:length(unique_z_true)
     temp_index=find(beta==unique_z_true(zz));
     z_true_smooth(temp_index)=unique_z_true(zz);
     if (unique_z_true(zz)==1)
         group_index_save_intercept(temp_index)=-0.5;
		 Intercept_smooth(temp_index)=-0.5;
     end
	 if (unique_z_true(zz)==-1)
	    lon_indicator=lon(temp_index);
		lat_indicator=lat(temp_index);
		cut=find(lon_indicator<0.5&lat_indicator<0.7);
		add=find(lon_indicator<0.46&lon_indicator>0.4&lat_indicator<0.73&lat_indicator>0.7);
	    cut=[cut;add];
		temp_index_1=temp_index(cut);
		temp_index(cut)=[];
		temp_index_2=temp_index;
		Intercept_smooth(temp_index_1)=1;
		Intercept_smooth(temp_index_2)=-0.5;
		group_index_save_intercept(temp_index_1)=1;
		group_index_save_intercept(temp_index_2)=-0.5;
	 end
	 if (unique_z_true(zz)==-0.5)
	     group_index_save_intercept(temp_index)=0.5;
         Intercept_smooth(temp_index)=0.5;
     end
	 if (unique_z_true(zz)==0.5)
	     group_index_save_intercept(temp_index)=-1;
         Intercept_smooth(temp_index)=-1;
     end
   end
   beta=z_true_smooth;
   Intercept=Intercept_smooth;
else
   pro=1;
   group_index_save=beta(:,which_pattern);
   group_index_save_intercept=beta(:,which_pattern);
   beta=beta(:,which_pattern).*pro;  
   z_true_smooth=beta;
   unique_z_true=unique(beta);
   for zz=1:length(unique_z_true)
     temp_index=find(beta==unique_z_true(zz));
     z_true_smooth(temp_index)=unique_z_true(zz)+((lon(temp_index))+(lat(temp_index))).^2;
     if (unique_z_true(zz)==1)
         group_index_save_intercept(temp_index)=1;
		 Intercept_smooth(temp_index)=-0.5+((lon(temp_index))+(lat(temp_index))).^1.5;
     end
	 if (unique_z_true(zz)==-1)
	    lon_indicator=lon(temp_index);
		lat_indicator=lat(temp_index);
		cut=find(lon_indicator<0.5&lat_indicator<0.7);
		add=find(lon_indicator<0.46&lon_indicator>0.4&lat_indicator<0.73&lat_indicator>0.7);
	    cut=[cut;add];
		temp_index_1=temp_index(cut);
		temp_index(cut)=[];
		temp_index_2=temp_index;
		Intercept_smooth(temp_index_1)=1+((lon(temp_index_1))+(lat(temp_index_1))).^2;
		Intercept_smooth(temp_index_2)=-0.5+((lon(temp_index_2))+(lat(temp_index_2))).^1.5;
		group_index_save_intercept(temp_index_1)=0;
		group_index_save_intercept(temp_index_2)=1;
	 end
	 if (unique_z_true(zz)==-0.5)
	     group_index_save_intercept(temp_index)=-1;
         Intercept_smooth(temp_index)=0.5+((lon(temp_index))+(lat(temp_index))).^1.5;
     end
	 if (unique_z_true(zz)==0.5)
	     group_index_save_intercept(temp_index)=-2;
         Intercept_smooth(temp_index)=-1+((lon(temp_index))+(lat(temp_index))).^1.7;
     end
   end
   beta=z_true_smooth;
   Intercept=Intercept_smooth;
end

beta=[beta Intercept'];

[x,y]=SCC_data_production(beta,epsilon_mag,sim_num,lon,lat,phi,sim_num*1000);

% save the simulation data
save('D:\scholar\new_data\one_cluster_m.mat','lon','lat','x','y','beta','group_index_save','group_index_save_intercept');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%run SCC and GWR on the simulation data
load('D:\scholar\new_data\one_cluster_m.mat');
sim_index=1:100

% Estimate the regression coefficients using the SCC method 
[beta_hat_SCC,MSE_SCC]=SCC_spatial_regression(x(sim_index,:,:),y(sim_index,:),lon,lat,beta,length(sim_index),[]);
	
% Estimate the regression coefficients using the GWR method  	
info.dtype='gaussian';
info.bmin=0.01^2;
info.bmax=2^2;
[beta_hat_GWR,MSE_GWR]=GWR_spatial_regression(x(sim_index,:,:),y(sim_index,:),lon,lat,beta,length(sim_index),info);
    
MSE_SCC(:,1)
MSE_SCC(:,2)

MSE_GWR(:,1)
MSE_GWR(:,2)
	
mean(MSE_SCC(:,1))
aa=MSE_SCC(:,1);
MSE_var=sqrt((mean((aa).^2)-(mean(aa)).^2)/length(aa))
	
mean(MSE_SCC(:,2))
aa=MSE_SCC(:,2);
MSE_var=sqrt((mean((aa).^2)-(mean(aa)).^2)/length(aa))
	
mean(MSE_GWR(:,1))
aa=MSE_GWR(:,1);
MSE_var=sqrt((mean((aa).^2)-(mean(aa)).^2)/length(aa))
	
mean(MSE_GWR(:,2))
aa=MSE_GWR(:,2);
MSE_var=sqrt((mean((aa).^2)-(mean(aa)).^2)/length(aa))


%save the results of SCC and GWR	
save('D:\beta_hat_SCC.mat','beta_hat_SCC')	
save('D:\beta_hat_GWR.mat','beta_hat_GWR')
	
	