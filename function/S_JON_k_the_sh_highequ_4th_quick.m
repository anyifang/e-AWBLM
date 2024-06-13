function S_k_the_JON_sh_highequ = S_JON_k_the_sh_highequ_4th_quick(k,the,U10,H,x_fetch_,C_beta,type,k_store)
% 浅水波浪谱，JON谱直接加上浅水函数psai，考虑4次波-波非线性
% 无量纲谱峰频率
fp_ = 3.5 * x_fetch_^(-0.33);
fp = max(0.13,fp_) * 9.81 / U10; % 谱峰频率
kp = k_calcu_store(fp,H,k_store); % 计算得到水深对应的谱峰波数
kp_ = kp * U10^2 / 9.81;

if strcmp(type,'XY2021' )
    %XY2021
    alpha = 0.0078 * kp_^0.49;
    gama = 1;
    
    % 2016Takagaki
elseif strcmp(type,'TK2016' )
%     Cp =  sqrt(tanh(kp * H) * 9.81 ./ kp);
%     
%     U10LIMIT =40;
%     u_staturated = min(40*tanh(u10/35),u10);
%     u_staturated = 40*tanh(u10/35);
%     fp = max(0.13,fp_) * 9.81 / u_staturated; % 谱峰频率
%     kp = k_calcu(fp,H); 
%     kp_ = kp * u_staturated^2 / 9.81;
%     kp_ = kp * u10^2 / 9.81;
    alpha = 0.0091 * kp_^0.49;
  %     alpha = 0.00136 * kp_;
%     gama_normal = 3.851 * (u10./Cp).^0.4047;
%     weight = (40 - u_staturated)/40;
    gama =1;
    %gama =1 *weight+(1-weight)*( 1.242 *gama_normal .* exp ( -0.7266*u10/U10LIMIT ));
end

siga = 0.07;
sigb = 0.09;

% 求出k截断频率
% kc =2 * kp;%3 * kp; 0.5; % 暂时先认为截断波数是谱峰波数的3倍关系
% fc = sqrt(9.81/4/pi^2 * kc * tanh(kc*H));  % 对应输入k的频率
%fc = 3*fp; % 暂时先认为截断波数是谱峰波数的5倍关系
fc = 10*fp; % 暂时先认为截断波数是谱峰波数的10倍关系
% fc =2.5* (log(H)+1)*fp;
kc = k_calcu_store(fc,H,k_store);% 暂时先认为截断波数是谱峰波数的10倍关系

f = sqrt(9.81/4/pi^2 * k * tanh(k*H)); % 对应输入k的频率
if f <= fp
    delta1 = exp(-(f-fp)^2/(2*siga^2*fp^2));
    S_f_JON = alpha * 9.81^2 * (2*pi)^(-4) * f^(-5) * exp(-5/4*(fp/f)^4) * gama ^ delta1;
    S_k_JON_sh_highequ = S_f_JON * (tanh(k*H))^3 * 9.81/8/pi^2/f / k; % 在原来转换基础上再除以k，波数谱psai(k)kdk=psai(f)df
    if abs(the) <= pi/2
        S_k_the_JON_sh_highequ = S_k_JON_sh_highequ * 2/pi * (cos(the))^2;
    else
        S_k_the_JON_sh_highequ = 0;
    end
elseif f > fp && f <= fc
    delta2 = exp(-(f-fp)^2/(2*sigb^2*fp^2));
    S_f_JON = alpha * 9.81^2 * (2*pi)^(-4) * f^(-5) * exp(-5/4*(fp/f)^4) * gama ^ delta2;
    S_k_JON_sh_highequ = S_f_JON * (tanh(k*H))^3 * 9.81/8/pi^2/f / k;
    if abs(the) <= pi/2
        S_k_the_JON_sh_highequ = S_k_JON_sh_highequ * 2/pi * (cos(the))^2;
    else
        S_k_the_JON_sh_highequ = 0;
    end
elseif f > fc
    % 求出截断频率下的波数谱值
    delta2 = exp(-(fc-fp)^2/(2*sigb^2*fp^2));
    psai_low_fc = alpha * 9.81^2 * (2*pi)^(-4) * fc^(-5) * exp(-5/4*(fp/fc)^4) * gama^delta2;
    psai_low_kc = psai_low_fc * (tanh(kc*H))^3 * 9.81/8/pi^2/fc / kc; % 方向积分下的截断频率处的波谱值
    
    % 接下来求解Δ复杂表达式部分
    C_theta_D = 3/8*pi;
    h_theta_1 = 2;
    deltafun = 4*kc^(-3)/(C_beta*C_theta_D* (kc*tanh(kc*H))^0.5) * h_theta_1 / psai_low_kc;  % 2.63前面
    % 接下来计算积分kc-k范围内的
    fun_k = @(x) (x.*(tanh(x*H))).^(-0.5);
    fun_jifen = integral(fun_k, kc, k);
    % 接下来求出高频谱值
    S_k_JON_sh_highequ = 8/(C_beta*C_theta_D)*k^(-3) / ((deltafun+2*fun_jifen)*(k*tanh(k*H))^0.5); %
    %disp( S_k_JON_sh_highequ)
    %S_k_JON_sh_highequ = 4 *  h_theta_1/(C_beta*C_theta_D)*k^(-3) / ((deltafun+2*fun_jifen)*(k*tanh(k*H))^0.5);
    %disp( S_k_JON_sh_highequ)
    if abs(the) <= pi/2
        S_k_the_JON_sh_highequ = S_k_JON_sh_highequ * 2/pi * (cos(the))^2;  % %1/2 * cos(the) 这里的方向函数是2/pi*cos(the)^2
    else
        S_k_the_JON_sh_highequ = 0;
    end
end
