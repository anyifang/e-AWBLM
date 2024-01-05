function kp = k_calcu_store(fp,H,k_store)
% 根据频散关系反算波数
%k = logspace(log10(0.0001), log10(1000), 100000);
%load('D:\matlab\bin\m\aWBLM\logk\k_calcu_log.mat')
% k = [0:0.0001:0.1 0.01:0.0001:1  1.01:0.001:10  10.1:0.1:1000 1001:0.5:2500 ] ; %波数范围
% f = sqrt(9.81/4/pi^2 .* k .* tanh(k*H));
% dif=abs(fp - f);
% [~,I ] = min(dif);
% kp= k (I) ;
fp_n=round((log10(fp)+4)/0.0002)+1;
H_n= round(H);
fp_n = min(fp_n,30001);
H_n = min(H_n,500);
kp = k_store(fp_n,H_n);

end
