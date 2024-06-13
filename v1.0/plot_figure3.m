
%% plot figure 3 in ZHANG & Yu JPO

FontSize = 18;
FontSize_a = 24;
FontSizel = 14;
u10_arr_XY = 1:1: 71;
load('../data/Cd_storage_XY_1580500.mat');


Cd_storage_XY_1580500(Cd_storage_XY_1580500<0.0001) = NaN;
Cd_storage_breaking(Cd_storage_breaking<0.0001) = NaN;

RGB = othercolor('Spectral4');
close all
color_dif =0.25* [1 , 1 , 1];
RGB_masure=0*((RGB -0.25)<=0) + (RGB -0.25).*((RGB -0.25)>0);
linewidth = 1.5;
MarkerSize =5;
CapSize = 5 ;
figure_length = 500;
figure(1)
subx = 0.1;
suby= 0.1;
sub_dy = suby +0.5;
subplot(311)
%subplot('Position',[0.05,0.2,subx,suby])
load('../data/zhao_2015.mat')
err_zhao_2015 = zhao_2015(:,3) - zhao_2015(:,2);
errorbar(zhao_2015(:,1) ,zhao_2015(:,2) ,err_zhao_2015...
   ,'linestyle' ,'--','Marker','square','LineWidth',linewidth,'MarkerSize' , MarkerSize,'CapSize',CapSize,...
    'color',RGB_masure(0.75/4*256,:),"MarkerFaceColor",RGB_masure(0.75/4*256,:)); hold on 
load('../data/chen_2022.mat');
err_chen_2022 = chen_2022(:,3) - chen_2022(:,2);
errorbar(chen_2022(:,1) ,chen_2022(:,2) ,err_chen_2022...
    ,'linestyle' ,'--','Marker','square','LineWidth',linewidth,'MarkerSize' , MarkerSize,'CapSize',CapSize...
    ,'color',RGB_masure(1.25/4*256,:),"MarkerFaceColor",RGB_masure(1.25/4*256,:))
i= 1;
H_num=3;
 n= fix((log10(min(max(15,10),500))-1)/ (log10(500)-1)*255)+1;
plot(u10_arr_XY,Cd_storage_XY_1580500(i,:),'--','linewidth',3,'Color',RGB(n,:)) ; hold on
plot(u10_arr,Cd_storage_breaking(i,:),'-','linewidth',3,'Color',RGB(n,:)) ; hold on
grid on
hl=legend({'Zhao (2015) {\itd}=14m','Chen (2022) {\itd}=16m','Xu and Yu (2021) {\itd}=15m','e-AWBLM {\itd}=15m'},'FontSize',FontSizel);  %,'NumColumns',2 'Orientation','horizontal'
%set(hl,'box','off')
ylim([0,5]*10^-3);xlim([0,70])
ylabel('{\itC_d}'); xlabel('{\itU}_1_0 (m/s)')
set(gca,"FontName","Times New Roman","FontSize",FontSize,"LineWidth",1);
text(1,0.0046,'\bf(a)','FontSize',FontSize_a,"FontName","Times New Roman");

subplot(312)
%subplot('Position',[0.05,0.2,subx,suby+sub_dy])
load('../data/jarosz_2007.mat')
x = [ jarosz_2007(:,3)'  jarosz_2007(end:-1:1,5)'];
y = [ jarosz_2007(:,4)' jarosz_2007(end:-1:1,6)'];
p=fill(x,y,'r','HandleVisibility','off');hold on
p.FaceColor = RGB_masure(2.5/4*256,:);
set(p,'edgealpha',0.5,'facealpha',0.5)
p.EdgeColor = 'none';
plot( jarosz_2007(:,1) , jarosz_2007(:,2),'--','LineWidth',linewidth,'MarkerSize' , MarkerSize,'color',RGB_masure(2.5/4*256,:),"MarkerFaceColor",RGB_masure(3/4*256,:)); hold on 
i= 2;
 n= fix((log10(min(max(80,10),500))-1)/ (log10(500)-1)*255)+1;


plot(u10_arr_XY,Cd_storage_XY_1580500(i,:),'--','linewidth',3,'Color',RGB(n,:)) ; hold on
plot(u10_arr,Cd_storage_breaking(i,:),'-','linewidth',3,'Color',RGB(n,:)) ; hold on
grid on
hl=legend({'Jarosz (2007) {\itd}=69–89m','Xu and Yu (2021) {\itd}=80m','e-AWBLM {\itd}=80m'},'FontSize',FontSizel);  %,'NumColumns',2 'Orientation','horizontal'
%set(hl,'box','off')
ylim([0,5]*10^-3);xlim([0,70])
ylabel('{\itC_d}'); xlabel('{\itU}_1_0 (m/s)')
set(gca,"FontName","Times New Roman","FontSize",FontSize,"LineWidth",1);
text(1,0.0046,'\bf(b)','FontSize',FontSize_a,"FontName","Times New Roman");

subplot(313)
% subplot('Position',[0.05,0.2,subx,suby + 2 * sub_dy])
i= 3;
load('../data/powell_2003.mat')
err_powell_2003pos = powell_2003(:,3) - powell_2003(:,2);
err_powell_2003neg = powell_2003(:,2) - powell_2003(:,4);
errorbar(powell_2003(:,1) ,powell_2003(:,2) ,err_powell_2003neg,err_powell_2003pos...
    ,'linestyle' ,'--','Marker','square','LineWidth',linewidth,'MarkerSize' , MarkerSize,'CapSize',CapSize,'color',RGB_masure(4/4*256,:),"MarkerFaceColor",RGB_masure(4/4*256,:)); hold on 
grid on
ylabel('Cd'); xlabel('U_1_0 m/s')
ylim([0,5]*10^-3);xlim([0,70])
 n= fix((log10(min(max(500,10),500))-1)/ (log10(500)-1)*255)+1;
plot(u10_arr_XY,Cd_storage_XY_1580500(end,:),'--','linewidth',3,'Color',RGB(n,:)) ; hold on
plot(u10_arr,Cd_storage_breaking(end,:),'-','linewidth',3,'Color',RGB(n,:)) ; hold on
hl=legend({'Powell (2003) deep water','Xu and Yu (2021) deep water','e-AWBLM deep water'},'FontSize',FontSizel);  %,'NumColumns',2 'Orientation','horizontal'
%set(hl,'box','off')
set(figure(1),'Position',[0,0,1.5,4.3]*figure_length)
ylabel('{\itC_d}'); xlabel('{\itU}_1_0 (m/s)')
set(gca,"FontName","Times New Roman","FontSize",FontSize,"LineWidth",1);
text(1,0.0046,'\bf(c)','FontSize',FontSize_a,"FontName","Times New Roman");

 img=gcf;
 print(img,'-dtiff','-r1200','./figure3.tif')
saveas(figure(1),'fig3.fig')

