% write by Zhang and Yu JPO(2024)
function kp = k_calcu(fp,H)
% calculate k According to wave the law of dispersion
k = logspace(log10(0.0001), log10(1000), 100000);
%load('D:\matlab\bin\m\aWBLM\logk\k_calcu_log.mat')
%k = [0:0.0001:0.1 0.01:0.0001:1  1.01:0.001:10  10.1:0.1:1000 1001:0.5:2500 ] ; %²¨Êý·¶Î§
f = sqrt(9.81/4/pi^2 .* k .* tanh(k*H));
dif=abs(fp - f);
[~,I ] = min(dif);
kp= k (I) ;
% load ('D:\matlab\bin\m\aWBLM\k_store.mat')
% fp_n=round(fp/0.001)+1;
% H_n= round(H);
% kp = k_store(fp_n,H_n);

end

% for i = 1:length(fp)
%     dif=abs(fp(i) - f);
%     [~,I ] = min(dif);
%     kp(i) = k (I) ;
% end

