%This is for empirical analysis.

%SCC and GWR methods.

%The codes is run in Matlab R2016a, windows system.

addpath(genpath('D:\scholar\glmnet_matlab'))
addpath(genpath('D:\scholar\jplv7'))
addpath(genpath('D:\scholar\SCC_program'))

    %Generate Figure 3
    load('D:\data_for_sever_south.mat');
	
    subplot(2,1,1)
    scatter(latn,-depthn,20,temp,'fill');
    caxis([min(temp),10]);
	set(gca,'xtick',[-0.5 -0.375 -0.25 -0.125 0]);
    set(gca,'xticklabel',{'60S','45S','30S','15S','Eq'});
	set(gca,'ytick',[-50/55 -40/55 -30/55 -20/55 -10/55 0]);
    set(gca,'yticklabel',{'-5000','-4000','-3000','-2000','-1000','0'});
	title('(a)temperature') 
    ylabel('depth(m)') 
	colormap('jet');
	d = colorbar;
	set(d,'xtick',[0 5 10])
	set(d,'xticklabel',{'0', '5', '10~27'})
    subplot(2,1,2)
	scatter(latn,-depthn,20,salt,'fill');
    caxis([min(salt),35.5]);
	set(gca,'xtick',[-0.5 -0.375 -0.25 -0.125 0]);
    set(gca,'xticklabel',{'60S','45S','30S','15S','Eq'});
	set(gca,'ytick',[-50/55 -40/55 -30/55 -20/55 -10/55 0]);
    set(gca,'yticklabel',{'-5000','-4000','-3000','-2000','-1000','0'});
	title('(b)salinity')
    ylabel('depth(m)') 
	colormap('jet');
	c = colorbar;
	set(c,'xtick',[34 34.5 35 35.5])
	set(c,'xticklabel',{'34', '34.5', '35', '35.5~37.2'})
	
    %SCC and GWR method.
    [beta_hat_SCC,TS_SCC,TS_GWR,BIC_SCC]=SCC_WOA_data_application(latn,depthn,temp,salt,0);
    save('D:\beta_hat_SCC.mat','beta_hat_SCC')
	
	%Generate Figure4(ii)
	TS_SCC_group=TS_SCC;
	TS_unique=unique(TS_SCC);
	rng(2);
	rr=randperm(length(TS_unique));
	for i=1:length(TS_unique)
	index=(TS_SCC==TS_unique(i));
	TS_SCC_group(index)=rr(i);
	end
	subplot(2,1,1)
    scatter(latn,-depthn,20,TS_SCC_group,'fill');
    caxis([min(TS_SCC_group),max(TS_SCC_group)]);
	colormap('jet');
	set(gca,'xtick',[-0.5 -0.375 -0.25 -0.125 0]);
    set(gca,'xticklabel',{'60S','45S','30S','15S','Eq'});
	set(gca,'ytick',[-50/55 -40/55 -30/55 -20/55 -10/55 0]);
    set(gca,'yticklabel',{'-5000','-4000','-3000','-2000','-1000','0'});
	title('(b) SCC') 
    ylabel('depth(m)')
    subplot(2,1,2)
	%Here we load the results from SCC*,i.e.,scc_results_scad2.mat, use svc_results_scad2.mat to load the results from SCVC. 
	load('D:/scholar/SCC_program/new_data/scc_results_scad2.mat')
    beta_S_group=beta_estimate_S;
	beta_unique=unique(beta_estimate_S);
	rng(5);
	rr=randperm(length(beta_unique));
	for i=1:length(beta_unique)
	index=(beta_estimate_S==beta_unique(i));
	beta_S_group(index)=rr(i);
	end
	beta_svc=beta_S_group;
	scatter(latn,-depthn,20,beta_svc,'fill');
    caxis([min(beta_svc),max(beta_svc)]);
	colormap('jet');
	set(gca,'xtick',[-0.5 -0.375 -0.25 -0.125 0]);
    set(gca,'xticklabel',{'60S','45S','30S','15S','Eq'});
	set(gca,'ytick',[-50/55 -40/55 -30/55 -20/55 -10/55 0]);
    set(gca,'yticklabel',{'-5000','-4000','-3000','-2000','-1000','0'});
	title('(c) SCC*') 
    ylabel('depth(m)')
    
	
	%Generate Figure4(i), also needs to load different .mat files to create the Figure of different methods.
	load('D:/scholar/SCC_program/new_data/svc_results_scad2.mat')
	upper_col=max(max(TS_SCC),max(beta_svc));
	low_col=min(min(TS_SCC),min(beta_svc));
	subplot(2,1,1)
	load('D:/scholar/SCC_program/new_data/tps_results.mat')
	beta_svc=beta_svc;
	scatter(latn,-depthn,20,beta_svc,'fill');
    caxis([low_col,upper_col]);
    colormap('jet');
	set(gca,'xtick',[-0.5 -0.375 -0.25 -0.125 0]);
    set(gca,'xticklabel',{'60S','45S','30S','15S','Eq'});
	set(gca,'ytick',[-50/55 -40/55 -30/55 -20/55 -10/55 0]);
    set(gca,'yticklabel',{'-5000','-4000','-3000','-2000','-1000','0'});
	title('(e) PSE') 
    ylabel('depth(m)') 
	colormap('jet');
	colorbar;
	subplot(2,1,2)
    scatter(latn,-depthn,20,TS_GWR,'fill');
    caxis([low_col,upper_col]);
	set(gca,'xtick',[-0.5 -0.375 -0.25 -0.125 0]);
    set(gca,'xticklabel',{'60S','45S','30S','15S','Eq'});
	set(gca,'ytick',[-50/55 -40/55 -30/55 -20/55 -10/55 0]);
    set(gca,'yticklabel',{'-5000','-4000','-3000','-2000','-1000','0'});
	title('(d) GWR') 
    ylabel('depth(m)') 
	colormap('jet');
	colorbar;
    

	
	
	
	