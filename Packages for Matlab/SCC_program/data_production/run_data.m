% This script is used to create the real and simulated data files in this package.
% These data files have been already included in the subfolder "data" in
% this package. So it is not necessary to run this script to make the data
% files available.

% Please first download the WOA datasets woa13_decav_t00_04v2.nc and
% woa13_decav_s00_04v2.nc and save the files in this subfolder.
% The WOA datasets can be downloaded from the following URLs.
%https://data.nodc.noaa.gov/thredds/catalog/woa/WOA13/DATAv2/temperature/netcdf/decav/0.25/catalog.html?dataset=woa/WOA13/DATAv2/temperature/netcdf/decav/0.25/woa13_decav_t00_04v2.nc
%https://data.nodc.noaa.gov/thredds/catalog/woa/WOA13/DATAv2/salinity/netcdf/decav/0.25/catalog.html?dataset=woa/WOA13/DATAv2/salinity/netcdf/decav/0.25/woa13_decav_s00_04v2.nc
%用这个
%---------------------------------------------------
clear all 
% 调用函数
addpath(genpath('D:\scholar\glmnet_matlab'))
addpath(genpath('D:\scholar\jplv7'))
addpath(genpath('D:\scholar\SCC_program'))

% simulated data
coef_linear=0.75; % correlation coefficient between the first two covariates
sim_num=100; % simulation numbers
%scenario='smooth_cluster'
scenario='smooth_cluster';
pattern='nonperfect';

epsilon_mag=0.1; % standard deviation (sigma) for iid random noise


if strcmp(pattern,'perfect')
   beta_threshold=0.483
   lon=repmat([linspace(0.01,beta_threshold,15),linspace(0.517,0.99,15)],1,30);
   lat=repmat(linspace(0.01,0.99,30),30,1);
   lat=lat(:)';
   lon=lon';
   lat=lat';
   beta=ones(length(lon),1);
   beta(lon<=beta_threshold)=1;
   beta(lon>beta_threshold)=0.5;
   which_pattern=1;
else
   %load('cluster_pre.mat');
   load('C:\Users\me\Desktop\可能要做的新论文\program for this new paper\svc_pattern.mat')
   which_pattern=1;
end




if strcmp(scenario,'constant_cluster')

   pro=1;
   group_index_save=beta(:,which_pattern);
   beta=beta(:,which_pattern).*pro;
else
   pro=1;
   group_index_save=beta(:,which_pattern);
   beta=beta(:,which_pattern).*pro;  
   z_true_smooth=beta;
   unique_z_true=unique(beta);
   for zz=1:length(unique_z_true)
     temp_index=(beta==unique_z_true(zz));
     z_true_smooth(temp_index)=unique_z_true(zz)+((lon(temp_index))+(lat(temp_index))).^2;
     if (unique_z_true(zz)==1)|(unique_z_true(zz)==-0.5)
         Intercept_smooth(temp_index)=1+((lon(temp_index))+(lat(temp_index))).^2;
     else
         Intercept_smooth(temp_index)=-0.5+((lon(temp_index))+(lat(temp_index))).^2;
     end
   end
   beta=z_true_smooth;
   Intercept=Intercept_smooth;
end

beta=[beta Intercept'];

[x,y]=SCC_data_production(beta,epsilon_mag,sim_num,coef_linear,lon,lat,0.1,sim_num*10);
save('D:\scholar\new_data\one_cluster_m.mat','lon','lat','x','y','beta','group_index_save');

    load('D:\scholar\new_data\one_cluster_m.mat');
    sim_index=1:100
    % Estimate the regression coefficients using the SCC method 
    %one-dimensional[beta_hat_SCC,MSE_SCC]=SCC_spatial_regression(x(sim_index,:,:)',y(sim_index,:),lon,lat,beta,1,[]);
	[beta_hat_SCC,MSE_SCC]=SCC_spatial_regression(x(sim_index,:,:),y(sim_index,:),lon,lat,beta,100,[]);
	
    save('C:\Users\me\Desktop\可能要做的新论文\program for this new paper\beta_hat_SCC.mat','beta_hat_SCC')
	
	SCC_RI_loop(beta,beta_hat_SCC)

	info.dtype='gaussian';
    info.bmin=0.01^2;
    info.bmax=2^2;
    %one-dimensional[beta_hat_GWR,MSE_GWR]=GWR_spatial_regression(x(sim_index,:,:)',y(sim_index,:),lon,lat,beta,1,info);
    [beta_hat_GWR,MSE_GWR]=GWR_spatial_regression(x(sim_index,:,:),y(sim_index,:),lon,lat,beta,100,info);
    
	
	SCC_RI_loop(beta,beta_hat_GWR)

	
	MSE_SCC
	MSE_GWR
	
	mean(MSE_SCC(:,1))
	mean(MSE_SCC(:,2))
	
	mean(MSE_GWR(:,1))
	mean(MSE_GWR(:,2))
	
	
	length(unique(beta_hat_SCC))
	
	mean((y(sim_index,:)-beta').^2)
	
	plot(sort(beta_hat_SCC),'*:')
	
    beta_hat_SCC=squeeze(beta_hat_SCC);
    beta_hat_GWR=squeeze(beta_hat_GWR);
	
	
    [n,p]=size(beta);
	x=ones(n,1);
	y=y(sim_index,:);
    options=[]
    

	load('colormap.mat');
	figure(1)	
    width=400;%宽度，像素数
    height=100;%高度
    left=200;%距屏幕左下角水平距离
    bottem=100;%距屏幕左下角垂直距离
    set(gcf,'position',[left,bottem,width,height])
    h1=subplot(1,3,1);
    scatter(lon,lat,20,beta,'fill');
	axis square;
	title('\beta')
    caxis([-5,5])
	originalSize1 = get(gca, 'Position')
    h2=subplot(1,3,2);
    scatter(lon,lat,20,beta_hat_SCC,'fill');
	title('\beta_{SCC}')
    caxis([-5,5])
	axis square;
	originalSize2 = get(gca, 'Position')
    h3=subplot(1,3,3);
    scatter(lon,lat,20,beta_hat_GWR,'fill');
	title('\beta_{GWR}')
    caxis([-5,5])
	axis square;
	originalSize3 = get(gca, 'Position')
    %colormap(cmap_scatterplot)
	colormap(jet(200))
	%colormap(cmap_scatterplot)
    colorbar
	set(h1, 'Position', originalSize1);
	set(h2, 'Position', originalSize2);
	set(h3, 'Position', originalSize3);
	
	
	load('D:/scholar/new_data/R.mat','Z_fit','z_true')
	load('colormap.mat');
	figure(1)	
    width=400;%宽度，像素数
    height=100;%高度
    left=200;%距屏幕左下角水平距离
    bottem=100;%距屏幕左下角垂直距离
    set(gcf,'position',[left,bottem,width,height])
    h1=subplot(2,2,1);
    scatter(lon,lat,20,beta,'fill');
	axis square;
	title('\beta')
    caxis([-5,5])
	originalSize1 = get(gca, 'Position')
    h2=subplot(2,2,2);
    scatter(lon,lat,20,beta_hat_SCC,'fill');
	title('\beta_{SCC}')
    caxis([-5,5])
	axis square;
	originalSize2 = get(gca, 'Position')
    h3=subplot(2,2,3);
    scatter(lon,lat,20,beta_hat_GWR,'fill');
	title('\beta_{GWR}')
    caxis([-5,5])
	axis square;
	originalSize3 = get(gca, 'Position')
	h4=subplot(2,2,4);
    scatter(lon,lat,20,Z_fit,'fill');
	title('\beta_{SVC}')
    caxis([-5,5])
	axis square;
	originalSize4 = get(gca, 'Position')
    %colormap(cmap_scatterplot)
	%colormap(jet(200))
	colormap(cmap_scatterplot)
    colorbar
	set(h1, 'Position', originalSize1);
	set(h2, 'Position', originalSize2);
	set(h3, 'Position', originalSize3);
    set(h4, 'Position', originalSize4);


	
%calculate MSE

    % Estimate the regression coefficients using the SCC method     
    load('D:\scholar\new_data\one_cluster_m.mat');
	sim_number=100;
    [beta_hat_SCC,MSE_SCC]=SCC_spatial_regression(x,y,lon,lat,beta,sim_number,[]);
	
	
	SCC_RI_loop(group_index_save,beta_hat_SCC(:,:,1))
	mean(MSE_SCC(:,1))
	
	group_number=ones(sim_number,1);
    for i=1:sim_number
	    group_number(i)=length(unique(beta_hat_SCC(i,:,1)));
	end
	mean(group_number)
	    
	
    load('D:\scholar\new_data\one_cluster_m.mat');
	sim_number=100;
    % Estimate the regression coefficients using the GWR method
    info.dtype='gaussian';
    info.bmin=0.01^2;
    info.bmax=2^2;
    
    [beta_hat_GWR,MSE_GWR]=GWR_spatial_regression(x,y,lon,lat,beta,sim_number,info);

	mean(MSE_GWR)
	mean(MSE_GWR(:,1))
 
 
 
 %for Plot examples
    load('colormap.mat');
	figure(1)	
    width=400;%宽度，像素数
    height=100;%高度
    left=200;%距屏幕左下角水平距离
    bottem=100;%距屏幕左下角垂直距离
    set(gcf,'position',[left,bottem,width,height])
    h1=subplot(1,1,1);
    scatter(lon,lat,20,beta(:,1),'fill');
	axis square;
	title('\beta')
    caxis([min(beta(:,1)),max(beta(:,1))])
	originalSize1 = get(gca, 'Position')
	colormap(jet(200))
	%colormap(cmap_scatterplot)
    colorbar
	set(h1, 'Position', originalSize1);

 