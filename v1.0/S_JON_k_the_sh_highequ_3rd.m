function S_k_the_JON_sh_highequ_3rd = S_JON_k_the_sh_highequ_3rd(k,the,u10,H,x_fetch_,C_beta,type)
% 制作平衡域方法下的浅水高频谱，只考虑3次波相互作用和频散关系
% 低频谱部分采用：浅水波浪谱，Young and Verhagen 1996
% 无量纲谱峰频率
fp_ = 3.5 * x_fetch_^(-0.33);
fp = max(0.13,fp_) * 9.81 / u10; % 谱峰频率
kp = k_calcu(fp,H); % 计算得到水深对应的谱峰波数 通过频散关系反求波束
kp_ = kp * u10^2 / 9.81;

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
%       alpha = 0.00136 * kp_;
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
 kc = k_calcu(fc,H); % 暂时先认为截断波数是谱峰波数的10倍关系
 
f = sqrt(9.81/4/pi^2 * k * tanh(k*H)); % 对应输入k的频率
%% 低频谱
if f <= fp
    delta1 = exp(-(f-fp)^2/(2*siga^2*fp^2));
    S_f_JON = alpha * 9.81^2 * (2*pi)^(-4) * f^(-5) * exp(-5/4*(fp/f)^4) * gama^delta1;
    S_k_JON_sh_highequ_3rd = S_f_JON * (tanh(k*H))^3 * 9.81/8/pi^2/f / k; % 在原来转换基础上再除以k，波数谱psai(k)kdk=psai(f)df
    if abs(the) <= pi/2 
        S_k_the_JON_sh_highequ_3rd = S_k_JON_sh_highequ_3rd * 2/pi * (cos(the))^2;
    else
        S_k_the_JON_sh_highequ_3rd = 0;
    end
elseif f > fp && f <= fc
    delta2 = exp(-(f-fp)^2/(2*sigb^2*fp^2));
    S_f_JON = alpha * 9.81^2 * (2*pi)^(-4) * f^(-5) * exp(-5/4*(fp/f)^4) * gama^delta2;
    S_k_JON_sh_highequ_3rd = S_f_JON * (tanh(k*H))^3 * 9.81/8/pi^2/f / k;
    if abs(the) <= pi/2 
        S_k_the_JON_sh_highequ_3rd = S_k_JON_sh_highequ_3rd * 2/pi * (cos(the))^2;
    else
        S_k_the_JON_sh_highequ_3rd = 0;
    end
elseif f > fc
    % 求出截断频率下的波数谱值
    delta2 = exp(-(fc-fp)^2/(2*sigb^2*fp^2));
    psai_low_fc = alpha * 9.81^2 * (2*pi)^(-4) * fc^(-5) * exp(-5/4*(fp/fc)^4) * gama^delta2;
    psai_low_kc = psai_low_fc * (tanh(kc*H))^3 * 9.81/8/pi^2/fc / kc; % 方向积分下的截断频率处的波谱值

    % 接下来求解fs(kΔ)部分
    C_theta_S = 16/15;
    h_theta_2 = pi/2;
    f_kc = log(sinh(kc * H)); 
    if isinf(f_kc)
        f_kc = kc * H - log(2);
    end
    %f_ks = 1/(C_beta * C_theta_S) * kc^(-2)*H / (kc*tanh(kc*H)) * h_theta_2 / psai_low_kc - f_kc;   %2.64
    f_ks = h_theta_2 * H /(C_beta * C_theta_S * kc^3 * tanh(kc*H) * psai_low_kc) - f_kc;            %2.64

    f_k = log(sinh(k*H));
    if isinf(f_k)
        f_k = k * H - log(2);
    end
   
    %% 接下来求出高频谱值
    S_k_JON_sh_highequ_3rd = 1/(C_beta * C_theta_S) * k^(-2)*H / ((f_ks + f_k)*k*tanh(k*H)) * h_theta_2;  % 2.54  有量纲的形式
    %   disp( S_k_JON_sh_highequ_3rd  )
     %S_k_JON_sh_highequ_3rd = k * H * h_theta_2 / (C_beta * C_theta_S * (f_ks + f_k) * tanh(k*H)) * k^(-4);  % 2.54 
    %  disp( S_k_JON_sh_highequ_3rd  )
    if abs(the) <= pi/2 
        S_k_the_JON_sh_highequ_3rd = S_k_JON_sh_highequ_3rd * 2/pi * (cos(the))^2; % 这里的方向函数是2/pi*cos(the)^2
    else
        S_k_the_JON_sh_highequ_3rd = 0;
    end    
end
