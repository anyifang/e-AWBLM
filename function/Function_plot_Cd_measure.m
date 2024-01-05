function  [p1,p2,p3,p4]=Function_plot_Cd_measure(RGB_masure)
linewidth = 1.5;
MarkerSize =4;
CapSize = 5 ;
figure_length = 400;
load('D:\MATLAB2021a\bin\m\aWBLM\Cd_measurement\zhao_2015.mat')
err_zhao_2015 = zhao_2015(:,3) - zhao_2015(:,2);
p1=errorbar(zhao_2015(:,1) ,zhao_2015(:,2) ,err_zhao_2015...
    ,'linestyle' ,'--','Marker','square','LineWidth',linewidth,'MarkerSize' , MarkerSize,'CapSize',CapSize,...
    'color',RGB_masure(0.75/4*256,:),"MarkerFaceColor",RGB_masure(0.75/4*256,:)); hold on
load('D:\MATLAB2021a\bin\m\aWBLM\Cd_measurement\chen_2022.mat');
err_chen_2022 = chen_2022(:,3) - chen_2022(:,2);
p2=errorbar(chen_2022(:,1) ,chen_2022(:,2) ,err_chen_2022...
    ,'linestyle' ,'--','Marker','square','LineWidth',linewidth,'MarkerSize' , MarkerSize,'CapSize',CapSize...
    ,'color',RGB_masure(1.25/4*256,:),"MarkerFaceColor",RGB_masure(1.25/4*256,:));
load('D:\MATLAB2021a\bin\m\aWBLM\Cd_measurement\jarosz_2007.mat')
x = [ jarosz_2007(:,3)'  jarosz_2007(end:-1:1,5)'];
y = [ jarosz_2007(:,4)' jarosz_2007(end:-1:1,6)'];
p=fill(x,y,'r','HandleVisibility','off');hold on
p.FaceColor = RGB_masure(2.5/4*256,:);
set(p,'edgealpha',0.5,'facealpha',0.5)
p.EdgeColor = 'none';
p3=plot( jarosz_2007(:,1) , jarosz_2007(:,2),'--','LineWidth',linewidth,'MarkerSize' , MarkerSize,'color',RGB_masure(2.5/4*256,:),"MarkerFaceColor",RGB_masure(3/4*256,:)); hold on 
load('D:\MATLAB2021a\bin\m\aWBLM\Cd_measurement\powell_2003.mat')
err_powell_2003pos = powell_2003(:,3) - powell_2003(:,2);
err_powell_2003neg = powell_2003(:,2) - powell_2003(:,4);
p4=errorbar(powell_2003(:,1) ,powell_2003(:,2) ,err_powell_2003neg,err_powell_2003pos...
    ,'linestyle' ,'--','Marker','square','LineWidth',linewidth,'MarkerSize' , ...
    MarkerSize,'CapSize',CapSize,'color',RGB_masure(4/4*256,:),"MarkerFaceColor",RGB_masure(4/4*256,:));
grid on
ylim([0 0.005]);xlim([0 70]);
end

