%This is matlab program for drawing Figure 2.

%The codes is run in Matlab R2016a, windows system.

clear all 
    %plot4.mat/plot5.mat is obtained in SCVC.r.  
	data1=load('D:/plot4.mat');
	data2=load('D:/plot5.mat');
	
	min_lim=min([min(data1.beta_slope),min(data1.beta_intercept),min(data2.beta_slope),min(data2.beta_intercept)]);
	max_lim=max([max(data1.beta_slope),max(data1.beta_intercept),max(data2.beta_slope),max(data2.beta_intercept)]);
	
	%load the data of MST-equal.
	data=load('D:/plot4.mat');
	
	psize=20;
	tsize=20;
	lsize=20;
	cc=15;
	lon=data.x;
	lat=data.y;
	beta=data.beta_slope;
	index=data.index;
	node_location=data.node_location;
	figure(1)	
    width=1600;
    height=400;
    left=200;
    bottem=100;
    set(gcf,'position',[left,bottem,width,height]);
    h1=subplot(1,4,1);
    scatter(lon,lat,psize,beta,'fill');
	hold on
	
	for i=1:size(index,1)
       plot(lon(index(i,:)),lat(index(i,:)),'r')	
	end
	hold off	
	axis square;
	xlabel('s_h','fontsize',lsize)
	ylabel('s_v','fontsize',lsize)
	title('\beta_2(s_i): MST-equal','fontsize',tsize)
    caxis([min_lim,max_lim])
	set(gca,'linewidth',2,'fontsize',cc,'fontname','Times');
	originalSize1 = get(gca, 'Position')

	
    h2=subplot(1,4,2);
	lon=data.x;
	lat=data.y;
	beta=data.beta_intercept;
	index=data.index;
	node_location=data.node_location;
    scatter(lon,lat,psize,beta,'fill');
	hold on
	for i=1:size(index,1)
       plot(lon(index(i,:)),lat(index(i,:)),'r')	
	end
	hold off	
	axis square;
	xlabel('s_h','fontsize',lsize)
	title('\beta_1(s_i): MST-equal','fontsize',tsize)
    caxis([min_lim,max_lim])
	set(gca,'linewidth',2,'fontsize',cc,'fontname','Times');
	originalSize2 = get(gca, 'Position')
	
	%load the data of MST-unequal.
	data=load('D:/plot5.mat');

    h3=subplot(1,4,3);
	lon=data.x;
	lat=data.y;
	beta=data.beta_slope;
	index=data.index;
	node_location=data.node_location;
    scatter(lon,lat,psize,beta,'fill');
	hold on
	for i=1:size(index,1)
       plot(lon(index(i,:)),lat(index(i,:)),'r')	
	end
	hold off	
	axis square;
	xlabel('s_h','fontsize',lsize)
	title('\beta_2(s_i): MST-unequal','fontsize',tsize)
    caxis([min_lim,max_lim])
	set(gca,'linewidth',2,'fontsize',cc,'fontname','Times');
	originalSize3 = get(gca, 'Position')
	
	
    h4=subplot(1,4,4);
	lon=data.x;
	lat=data.y;
	beta=data.beta_intercept;
	index=data.index;
	node_location=data.node_location;
    scatter(lon,lat,psize,beta,'fill');
	hold on
	for i=1:size(index,1)
       plot(lon(index(i,:)),lat(index(i,:)),'r')	
	end
	hold off	
	axis square;
	xlabel('s_h','fontsize',lsize)
	title('\beta_1(s_i): MST-unequal','fontsize',tsize)
    caxis([min_lim,max_lim])
	set(gca,'linewidth',2,'fontsize',cc,'fontname','Times');
	originalSize4 = get(gca, 'Position')
	
	
	colormap(jet(64))
    colorbar('fontsize',15)
	set(h1, 'Position', originalSize1);
	set(h2, 'Position', originalSize2);
	set(h3, 'Position', originalSize3);
	set(h4, 'Position', originalSize4);
	