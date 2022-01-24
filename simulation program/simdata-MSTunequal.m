%This is for generating simulation data of MST-unequal.  

%Give the results for SCC and GWR, using the program provided by Li and Sang (2019).

%The codes is run in Matlab R2016a, windows system.

clear all 

%load the program for SCC and GWR, they are provided in the file "Packages for Matlab".
%you need to specify the right path.
addpath(genpath('D:\Program\Packages for Matlab\glmnet_matlab'))
addpath(genpath('D:\Program\Packages for Matlab\jplv7'))
addpath(genpath('D:\Program\Packages for Matlab\SCC_program'))


% simulated data
sim_num=100; % simulation numbers
epsilon_mag=0.3; % standard deviation (sigma) for iid random noise
phi=1;%range parameter for the covariates
which_pattern=1;

% load the locations of MST-unequal, already saved in the file "location generation".
load('D:\Program\simulation program\try_nonperfect_svc_pattern.mat');
 
% If beta(s) is smooth within each subregion, using "smooth_cluster". If constant, using "constant_cluster".
scenario='smooth_cluster';

% give values of beta on the locations
   pro=1;
   group_index_save=beta(:,which_pattern);
   group_index_save_intercept=beta(:,which_pattern);
   beta=beta(:,which_pattern).*pro;  
   z_true_smooth=beta;
   unique_z_true=unique(beta);
   for zz=1:length(unique_z_true)
     temp_index=find(beta==unique_z_true(zz));
     
     if (unique_z_true(zz)==1)
	    lon_indicator=lon(temp_index);
		lat_indicator=lat(temp_index);
		
		Intercept_smooth(temp_index)=1;
		group_index_save_intercept(temp_index)=1;
		
		cut=find(lon_indicator>0.34);
		temp_index_1=temp_index(cut);
		temp_index(cut)=[];
		temp_index_2=temp_index;
		
		group_index_save(temp_index_1)=-1;
		group_index_save(temp_index_2)=1;
 
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
		
		Intercept_smooth(temp_index_1)=3;
		Intercept_smooth(temp_index_2)=1;
		group_index_save_intercept(temp_index_1)=3;
		group_index_save_intercept(temp_index_2)=1;
		
		lon_indicator=lon(temp_index_1);
		lat_indicator=lat(temp_index_1);
		cut=find(lon_indicator<0.26&lon_indicator>0.17&lat_indicator>0.22&lat_indicator<0.27);
		temp_index_3=temp_index_1(cut);
	    
		Intercept_smooth(temp_index_3)=-0.5;
		group_index_save_intercept(temp_index_3)=-0.5;
	    group_index_save(temp_index_3)=0.5;
		
		lon_indicator=lon(temp_index_2);
		lat_indicator=lat(temp_index_2);
		cut=find(lat_indicator<0.6);
		temp_index_4=temp_index_2(cut);
		
		Intercept_smooth(temp_index_4)=-0.5;
		group_index_save_intercept(temp_index_4)=-0.5;
	    group_index_save(temp_index_4)=0.5;
		
		cut=find(lon_indicator<0.46&lat_indicator<0.82);
		cut1=find(lon_indicator<0.67&lon_indicator>0.5&lat_indicator<0.71);
		cut2=find(lon_indicator<0.6&lon_indicator>0.54&lat_indicator<0.8&lat_indicator>0.709);
		cut3=find(lon_indicator<0.64&lon_indicator>0.59&lat_indicator<0.74&lat_indicator>0.7);
		cut4=find(lon_indicator<0.6&lon_indicator>0.56&lat_indicator<0.82&lat_indicator>0.78);
		
		cut=[cut;cut1;cut2;cut3;cut4];
		temp_index_5=temp_index_2(cut);
		
		Intercept_smooth(temp_index_5)=3;
		group_index_save_intercept(temp_index_5)=3;
		
	 end
	 
	 if (unique_z_true(zz)==0.5)
	    lon_indicator=lon(temp_index);
		lat_indicator=lat(temp_index);
		cut=find(lat_indicator>0.83);
		temp_index_1=temp_index(cut);
		temp_index(cut)=[];
		temp_index_2=temp_index;
	    
	     group_index_save_intercept(temp_index_2)=-0.5;
         Intercept_smooth(temp_index_2)=-0.5;
	     group_index_save_intercept(temp_index_1)=1;
         Intercept_smooth(temp_index_1)=1;
		 group_index_save(temp_index_1)=-1;
     end
	 
	 if (unique_z_true(zz)==-0.5)
	    lon_indicator=lon(temp_index);
		lat_indicator=lat(temp_index);
		cut=find(lon_indicator<0.9&lon_indicator>0.8&lat_indicator<0.36&lat_indicator>0.32);
		temp_index_1=temp_index(cut);
		temp_index(cut)=[];
		temp_index_2=temp_index;
	    
		group_index_save_intercept(temp_index_1)=-0.5;
		group_index_save(temp_index_1)=0.5;
        Intercept_smooth(temp_index_1)=-0.5;
	    group_index_save_intercept(temp_index_2)=1.5;
        Intercept_smooth(temp_index_2)=1.5;
     end
	 
   end
   beta=group_index_save;
   Intercept=Intercept_smooth;

if strcmp(scenario,'constant_cluster')
   beta=group_index_save;
   pro=1;
   unique_z_true=unique(group_index_save_intercept);
   for zz=1:length(unique_z_true)
     temp_index=find(group_index_save_intercept==unique_z_true(zz));
     if (unique_z_true(zz)==1)
		 Intercept_smooth(temp_index)=-0.5;
     end
	 if (unique_z_true(zz)==3)
		Intercept_smooth(temp_index)=1;
	 end
	 if (unique_z_true(zz)==1.5)
         Intercept_smooth(temp_index)=0.5;
     end
	 if (unique_z_true(zz)==-0.5)
         Intercept_smooth(temp_index)=-1;
     end
   end
   Intercept=Intercept_smooth;
   
else
   beta=group_index_save+((lon)+(lat)).^2;
   pro=1;
   unique_z_true=unique(group_index_save_intercept);
   for zz=1:length(unique_z_true)
     temp_index=find(group_index_save_intercept==unique_z_true(zz));
     if (unique_z_true(zz)==1)
		 Intercept_smooth(temp_index)=-0.5+((lon(temp_index))+(lat(temp_index))).^1.5;
     end
	 if (unique_z_true(zz)==3)
		Intercept_smooth(temp_index)=1+((lon(temp_index))+(lat(temp_index))).^2;
	 end
	 if (unique_z_true(zz)==1.5)
         Intercept_smooth(temp_index)=0.5+((lon(temp_index))+(lat(temp_index))).^1.5;
     end
	 if (unique_z_true(zz)==-0.5)
         Intercept_smooth(temp_index)=-1+((lon(temp_index))+(lat(temp_index))).^1.7;
     end
   end
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
	
	