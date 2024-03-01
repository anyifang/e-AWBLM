function S_k_the_JON_sh_highequ_3rd = TMA_k_theta_shallow(k,the,u10,H,C_beta,gama,alpha,fp,k_store)
kp = k_calcu_store(fp,H,k_store); 
siga = 0.07;
sigb = 0.09;
fc = 10*fp;
 kc = k_calcu_store(fc,H,k_store);
 
f = sqrt(9.81/4/pi^2 * k * tanh(k*H)); 

if f <= fp
    delta1 = exp(-(f-fp)^2/(2*siga^2*fp^2));
    S_f_JON = alpha * 9.81^2 * (2*pi)^(-4) * f^(-5) * exp(-5/4*(fp/f)^4) * gama^delta1;
    S_k_JON_sh_highequ_3rd = S_f_JON * (tanh(k*H))^3 * 9.81/8/pi^2/f / k; 
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
    delta2 = exp(-(fc-fp)^2/(2*sigb^2*fp^2));
    psai_low_fc = alpha * 9.81^2 * (2*pi)^(-4) * fc^(-5) * exp(-5/4*(fp/fc)^4) * gama^delta2;
    psai_low_kc = psai_low_fc * (tanh(kc*H))^3 * 9.81/8/pi^2/fc / kc; 
    C_theta_S = 16/15;
    h_theta_2 = pi/2;
    f_kc = log(sinh(kc * H)); 
    if isinf(f_kc)
        f_kc = kc * H - log(2);
    end
    f_ks = h_theta_2 * H /(C_beta * C_theta_S * kc^3 * tanh(kc*H) * psai_low_kc) - f_kc; 
    f_k = log(sinh(k*H));
    if isinf(f_k)
        f_k = k * H - log(2);
    end
    S_k_JON_sh_highequ_3rd = 1/(C_beta * C_theta_S) * k^(-2)*H / ((f_ks + f_k)*k*tanh(k*H)) * h_theta_2; 
    if abs(the) <= pi/2 
        S_k_the_JON_sh_highequ_3rd = S_k_JON_sh_highequ_3rd * 2/pi * (cos(the))^2; 
    else
        S_k_the_JON_sh_highequ_3rd = 0;
    end    
end
