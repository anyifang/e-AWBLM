% e_AWBLM to estimate wind stress

% writed by Xu and Yu (2021)
% include the water depth effect on the Cd through the TMA wave-spectrum
% refer to Moon (2004) to solve the taov iteratively
% Zhang and Yu enhanced by adding breaking-induce stress 2023 / 4 /20
% The default unit is /kg /m /s

%  fix in 2023/10/1 by Zhang
%  to reduce the compute cost. the us has been changed with U10 in Eq(17)
%  the sea spray is fixed with kapa * U10 ^ 3 * $\int$ H（k）LAMDA (k)
% in this equation H(k)= B(k)^0.5 * k
% thc changed results also have a good agreement with measurement


% Contact with zayf21@mails.tsinghua.edu.cn if you have any problem with
% code or you have some methods to modify this code.
close all
clear
addpath('../data')
addpath('../function')
load ('k_store.mat')
%% Cd measurement
RGB = othercolor('YlGnBu9');
figure(3)
Function_plot_Cd_measure(RGB)
%% e_AWLBM
u10_arr = 1:5:71;
H_arr = [15 80 500];
fetch_ = 10^5; % wind fetch
Cd = zeros(length(u10_arr),1);
feixianxing = zeros(length(u10_arr),1);
for  n_H_arr = 1 : length(H_arr) % different water depth
    H = H_arr( n_H_arr);
    for n_u10_arr = 1 : length(u10_arr)% different wind velocity at 10-m level
        U10 = u10_arr(n_u10_arr) ;
        disp(['H=',num2str(H),' U10=',num2str(U10)])
        fp_ = 3.5 * fetch_^(-0.33);
        fp = max(0.13,fp_) * 9.81 / U10;   % peak frequency (Alves et al., 2003)
        kp = k_calcu_store(fp,H,k_store);
        kp_ = kp * U10^2 / 9.81;
        XI = 1;  % enhancement factor
        alpha = 0.0078 * kp_^0.49; %  Phillips parameter
        taow_max_Z = max(0.02/kp, 0.05/400) ;
        tao_vx = 0.001;   % initial viscous stress
        tao_vy = 0;
        
        u_star = U10/28; % initial u*
        k_min = 0.01; %0.0049 * 9.81 / (u10/28)^2;
        k_max = 400;
        k_arr = logspace(log10(k_min), log10(k_max), 100);
        theta_arr = -pi : 2*pi/36 : pi; % wave direction
        w_arr = sqrt(k_arr* 9.81 .* tanh(k_arr * H) ); % wave angular frequency
        c_arr = 9.81 ./ w_arr .* tanh(k_arr * H);  % wave velocity
        
        b_kmin = 0.008; % Melville, Veron, and White (2002) 0.007
        b_kmax= 0.02;%
        kapa = 8 * 10^-6; % empirical constant for sea spray
        C_beta = 25; % C_beta = 25 Chen and Yu 2016
        C_b = 0.35; % pressure drop coefficient, Kudryavtsev（2014）
        
        delta = 0.05; % decay length scale of the breaking-wave induced stress
        rou_w = 10^3; % water density
        rou_a = 1.293; % air density
        nu_a = 14.8*10^(-6); % coefficient of kinematic viscosity
        g = 9.81;
        
        
        fai_k_the_arr_3rd = zeros(length(theta_arr), length(k_arr));
        fai_k_the_arr_4th = zeros(length(theta_arr), length(k_arr));
        for i = 1 : length(theta_arr)
            for j = 1 : length(k_arr)
                if abs(theta_arr(i)) <= pi/2
                    fai_k_the_arr_3rd(i, j) =  TMA_k_theta_shallow(k_arr(j),theta_arr(i),U10,H,C_beta,XI,alpha,fp,k_store);
                    fai_k_the_arr_4th (i, j) =  TMA_k_theta_deep(k_arr(j),theta_arr(i),U10,H,C_beta,XI,alpha,fp,k_store);
                else
                    fai_k_the_arr_3rd(i, j) = 0;
                    fai_k_the_arr_4th(i, j) = 0;
                end
            end
        end
        
        fai_k_the_arr = (fai_k_the_arr_4th+fai_k_the_arr_3rd)/2;
        W_quanzhong_before = 0.5;
        %% to calculate TMA wave spectrum in arbitrary water depth
        % this part is calculted based on XY(2021)
        for iteration_quanzhong = 1 : 200
            integral_1 = zeros(length(k_arr),1); % to integral in theta
            for i = 1 : length(k_arr)
                fai_arr = fai_k_the_arr(:, i)' * k_arr(i);
                integral_1(i) = sum(fai_arr(2:end-1))*2*pi/36 + (fai_arr(1)+fai_arr(end))*2*pi/36/2;
            end
            
            del_k_arr = zeros(length(k_arr),1);
            if length(k_arr) > 2
                for kk = 2 : length(k_arr)-1
                    del_k_arr(kk) = (k_arr(kk+1) - k_arr(kk-1))/2;
                end
            end
            del_k_arr(1) = (k_arr(2) - k_arr(1)) / 2;
            del_k_arr(end) = (k_arr(end) - k_arr(end-1)) / 2; % to integral in k
            integral_2 = sum(integral_1 .* del_k_arr);
            Hs_o = 4*sqrt(integral_2);
     
            f_arr = sqrt(9.81/4/pi^2 * k_arr .* tanh(k_arr*H));  
            integral_1 = zeros(length(f_arr),1); 
            for i = 1 : length(f_arr)
                fai_arr_f = fai_k_the_arr(:, i)' * k_arr(i) * 8*pi^2*f_arr(i) / (9.81*tanh(k_arr(i)*H) + 9.81*k_arr(i)*H*(sech(k_arr(i)*H))^2) * f_arr(i); % 对应k_integral波数的所有方向上的波谱值
                integral_1(i) = sum(fai_arr_f(2:end-1))*2*pi/36 + (fai_arr_f(1)+fai_arr_f(end))*2*pi/36/2;
            end
            del_f_arr = zeros(length(f_arr),1);
            if length(f_arr) > 2
                for ff = 2 : length(f_arr)-1
                    del_f_arr(ff) = (f_arr(ff+1) - f_arr(ff-1))/2;
                end
            end
            del_f_arr(1) = (f_arr(2) - f_arr(1)) / 2;
            del_f_arr(end) = (f_arr(end) - f_arr(end-1)) / 2;
            integral_2 = sum(integral_1 .* del_f_arr);
            Tm01 = (Hs_o/4)^2 / integral_2;
            Ur = 9.81/8/sqrt(2)/pi^2 * Hs_o * Tm01^2 / H^2;
            if Ur <= 0.25
                W_quanzhong = 1;
            elseif Ur >= 1
                W_quanzhong = 0;
            else
                W_quanzhong = ((1-Ur)/0.75);
            end
            W_quanzhong_after = W_quanzhong;
            if abs(W_quanzhong_after - W_quanzhong_before) <0.01
                break
            end
            W_quanzhong_before = W_quanzhong_after;
            fai_k_the_arr = fai_k_the_arr_4th * W_quanzhong_after + fai_k_the_arr_3rd * (1-W_quanzhong_after);
        end
        %% calculate gama_b
        %gama_b_mean =   [0.0625 0.0571 0.0501 0.0448 0.0339];% 10m 20m 50m 100m 500m
        gama_b_mean =   [ 0.0580  0.0456   0.0339];  % 15m 80m 500m
        % gama_b = gama_b_mean(n_H_arr)/gama_b_mean(end)*0.3;
        gama_b = gama_b_U_H_quick(U10,H,k_store)/gama_b_U_H_quick(U10,500,k_store)*0.3;
        gama_bmax = -0.19*log10(H)+0.83;
        gama_b = min(gama_bmax,gama_b);
        gama_b = min(0.6,gama_b);
        gama_b = max(0.2,gama_b);
        Weigth =W_quanzhong_after;  % 10m 100 1000m
        
        k_arr_breaking = gama_b /delta * k_arr;
        w_arr_breaking = sqrt(k_arr_breaking* 9.81 .* tanh(k_arr_breaking * H) );
        c_arr_breaking = 9.81 ./ w_arr_breaking .* tanh(k_arr_breaking * H);
        fai_k_the_arr_3rd_breaking = zeros(length(theta_arr), length(k_arr));
        fai_k_the_arr_4th_breaking = zeros(length(theta_arr), length(k_arr));
        for i = 1 : length(theta_arr)
            for j = 1 : length(k_arr)
                if abs(theta_arr(i)) <= pi/2
                    
                    fai_k_the_arr_3rd_breaking (i, j) =  TMA_k_theta_shallow(k_arr_breaking(j),theta_arr(i),U10,H,C_beta,XI,alpha,fp,k_store);
                    fai_k_the_arr_4th_breaking (i, j) =  TMA_k_theta_deep(k_arr_breaking(j),theta_arr(i),U10,H,C_beta,XI,alpha,fp,k_store);
                else
                    fai_k_the_arr_3rd_breaking(i, j) = 0;
                    fai_k_the_arr_4th_breaking(i, j) = 0;
                end
            end
        end
        fai_k_the_arr_breaking = fai_k_the_arr_4th_breaking * W_quanzhong_after + fai_k_the_arr_3rd_breaking * (1-W_quanzhong_after);   %Wave spectrum at arbitrary water depth
        
        %%  The main iterative part of the procedure such that U(10) = U10.
        for iteration_taov = 1 : 200
            % wave growth rate
            beta_g = zeros(length(k_arr), length(theta_arr));
            beta_g_breaking = zeros(length(k_arr), length(theta_arr));
            tao_tx_z = zeros(length(k_arr),1);
            tao_ty_z = zeros(length(k_arr),1);
            
            
            tao_v = sqrt(tao_vx^2 + tao_vy^2);
            if tao_v <= 10^-8
                cos_taov = 1; sin_taov = 0;
            else
                cos_taov = tao_vx/tao_v; sin_taov = tao_vy/tao_v;
            end
            h_the_arr = (cos(theta_arr)*cos_taov + sin(theta_arr)*sin_taov).^2;
            h_the_arr(abs(theta_arr) > pi/2) = 0;
            beta_g(end, :) = C_beta * w_arr(end) / rou_w / c_arr(end)^2 * tao_v * h_the_arr;
            beta_g_breaking(end, :) = C_beta * w_arr_breaking(end) / rou_w / c_arr_breaking(end)^2 * tao_v * h_the_arr;
            lamda_breaking_the_arr(:,length(k_arr)) =zeros(length(theta_arr),1) ;
            integral_1_lamda_k_breaking(length(k_arr)) = 0;
            tao_tx_z(end) = tao_vx;
            tao_ty_z(end) = tao_vy;
            zmin = delta / k_max;
            zmax = delta / k_min;
            z_arr = delta ./ k_arr;
            num1 = 20;
            num3 = 20;
            z_arr3 = zmax : (10-zmax)/20 : 10;
            zv = 0.1 * nu_a / sqrt(tao_v/rou_a);
            zv = min(zv , zmin /2);
            z_arr1 = zv : (zmin-zv)/num3 : zmin;
            z_arr_total = [z_arr1(1:end-1), z_arr(end:-1:2), z_arr3];
            iteration_n_arr = zeros(length(k_arr)-1,1); %  Store the number of iterations to see if there is any non-convergence
            for iteration_u = 1 : 100    %iterative part to calculate tao_b
                % Next iteratively solve for beta_g(k) and taot on each k.
                % Since tao and beta are functions of each other, solving for both requires iteration
                for i = length(k_arr)-1 : -1 : 1
                    for iteration = 1 : 200
                        iteration_before = beta_g(i,:);
                        % calculate lamda（k）based on  Sin = Sdb
                        if k_arr_breaking(i)>20*pi % Kudryavtsev（2001） tao_b = 0 with  k >20pi
                            lamda_breaking_the_arr(:,i) =zeros(length(theta_arr),1) ;
                            integral_1_lamda_k_breaking(i) = 0;
                        else
                            fai_arr_slope_breaking =(fai_k_the_arr_breaking(:, i)'  *  k_arr_breaking(i) ^4).^0.5;
                            integral_1_slope_breaking(i) =0.25*( sum(fai_arr_slope_breaking(2:end-1))*2*pi/36 + (fai_arr_slope_breaking(1)+fai_arr_slope_breaking(end))*2*pi/36/2)^0.5;
                            
                            %  Muller(2009)
                            del_k_arr_breaking = zeros(length(k_arr_breaking),1);
                            if length(k_arr_breaking) > 2
                                for kk = 2 : length(k_arr_breaking)-1
                                    del_k_arr_breaking(kk) = (k_arr_breaking(kk+1) - k_arr_breaking(kk-1))/2;
                                end
                            end
                            del_k_arr_breaking(1) = (k_arr_breaking(2) - k_arr_breaking(1)) / 2;
                            del_k_arr_breaking(end) = (k_arr_breaking(end) - k_arr_breaking(end-1)) / 2;
                            % It's all calculated and then integrated, because the iterations for beta and tau are backwards from kmax -> kmin.
                            for num_s = 1:length(k_arr)
                                S_Muller_k_the_breaking =fai_k_the_arr_breaking(:, num_s)'  *  k_arr_breaking(num_s) ^3;
                                S_Muller_k_breaking(num_s) =( sum(S_Muller_k_the_breaking(2:end-1))*2*pi/36 + (S_Muller_k_the_breaking(1)+S_Muller_k_the_breaking(end))*2*pi/36/2) ;
                            end
                            S_Muller_k_breaking_integral = S_Muller_k_breaking(1:i+1);
                            del_k_arr_s_breaking = del_k_arr_breaking(1:i+1);
                            
                            W_quanzhong_after = 1;
                            b_k_breaking_s(i)=min(max( 0.01*sum( S_Muller_k_breaking_integral' .* del_k_arr_s_breaking )^0.5, b_kmin),b_kmax);   % spilling
                            b_k_breaking_d(i)=min(max( 0.25*sum( 2*S_Muller_k_breaking_integral'.^0.5 .* del_k_arr_s_breaking )^2.5, b_kmin),b_kmax);   % plunging
                            b_k_breaking_d_fix(i)=min(max( 0.25*(2*sum( S_Muller_k_breaking_integral' .* del_k_arr_s_breaking )^0.5)^2.5, b_kmin),b_kmax);   % plunging
                            b_k_breaking(i)= W_quanzhong_after * b_k_breaking_d_fix(i) + (1-W_quanzhong_after) *b_k_breaking_d(i);
                            
                            lamda_breaking_the_arr(:,i) =beta_g_breaking(i, :).*fai_k_the_arr_breaking(:, i)' * g ^3 /  c_arr_breaking(i)^6 / w_arr_breaking(i) / b_k_breaking(i) ;
                            integral_1_lamda_k_breaking(i) = sum(lamda_breaking_the_arr(2:end-1,i))*2*pi/36 + (lamda_breaking_the_arr(1,i)+lamda_breaking_the_arr(end,i))*2*pi/36/2;
                        end
                        
                        
                        beta_g_integral = beta_g(i:end,:);
                        k_integral_arr = k_arr(i:end);
                        w_integral_arr = w_arr(i:end);
                        k_integral_arr_breaking =  k_arr_breaking (i:end);
                        w_integral_arr_breaking = w_arr_breaking(i:end);
                        
                        integralx_1 = zeros(length(k_integral_arr),1);
                        integraly_1 = zeros(length(k_integral_arr),1);
                        integralx_1_breaking = zeros(length(k_integral_arr),1);
                        integraly_1_breaking = zeros(length(k_integral_arr),1);
                        for ii = 1 : length(k_integral_arr)
                            k_integral = k_integral_arr(ii);
                            w_integral = w_integral_arr(ii);
                            k_integral_breaking = k_integral_arr_breaking(ii);
                            w_integral_breaking = w_integral_arr_breaking(ii);
                            c_integral_breaking  = w_integral_breaking/k_integral_breaking;
                            the_integral_arr = -pi : 2*pi/36 : pi;
                            
                            fai_arr = fai_k_the_arr(:, i+ii-1)';
                            integralx_arr = beta_g_integral(ii,:) .* fai_arr * rou_w * w_integral * k_integral .* cos(the_integral_arr);
                            integralx_1(ii) = sum(integralx_arr(2:end-1))*2*pi/36 + (integralx_arr(1)+integralx_arr(end))*2*pi/36/2;
                            integraly_arr = beta_g_integral(ii,:) .* fai_arr * rou_w * w_integral * k_integral .* sin(the_integral_arr);
                            integraly_1(ii) = sum(integraly_arr(2:end-1))*2*pi/36 + (integraly_arr(1)+integraly_arr(end))*2*pi/36/2;
                            lamda_arr =  lamda_breaking_the_arr(:, i+ii-1)';
                            
                            z_breaking = gama_b / k_integral_breaking;  % ak = gama
                            z0 = 10/(exp(1)-1);
                            if iteration_u == 1
                                ub = U10 * log((z_breaking+z0)/z0) * h_the_arr.^0.5 - c_integral_breaking;
                                ub = ub .* (ub > 0) ;
                            else
                                [MIN, I ]=min(abs (z_arr_total-z_breaking));
                                ub = (u_arr(I))* h_the_arr.^0.5  - c_integral_breaking;
                                ub = ub .* (ub > 0) ;
                            end
                            integralx_arr_breaking = rou_a * 2 * gama_b * C_b * ub.^2 .* lamda_arr .* cos(the_integral_arr) ;   %Kudryavtsev（2014）
                            integraly_arr_breaking = rou_a * 2 * gama_b * C_b * ub.^2 .* lamda_arr .* sin(the_integral_arr) ;   %Kudryavtsev（2014）
                            integralx_1_breaking(ii) = sum(integralx_arr_breaking(2:end-1))*2*pi/36 + (integralx_arr_breaking(1)+integralx_arr_breaking(end))*2*pi/36/2;
                            integraly_1_breaking(ii) = sum(integraly_arr_breaking(2:end-1))*2*pi/36 + (integraly_arr_breaking(1)+integraly_arr_breaking(end))*2*pi/36/2;
                        end
                        
                        del_k_arr = zeros(length(k_integral_arr),1);
                        if length(k_integral_arr) > 2
                            for kk = 2 : length(k_integral_arr)-1
                                del_k_arr(kk) = (k_integral_arr(kk+1) - k_integral_arr(kk-1))/2;
                            end
                        end
                        del_k_arr(1) = (k_integral_arr(2) - k_integral_arr(1)) / 2;
                        del_k_arr(end) = (k_integral_arr(end) - k_integral_arr(end-1)) / 2;
                        integralx_2 = sum(integralx_1 .* del_k_arr);
                        integraly_2 = sum(integraly_1 .* del_k_arr);
                        
                        del_k_arr_breaking = zeros(length(k_integral_arr_breaking),1);
                        if length(k_integral_arr_breaking) > 2
                            for kk = 2 : length(k_integral_arr_breaking)-1
                                del_k_arr_breaking(kk) = (k_integral_arr_breaking(kk+1) - k_integral_arr_breaking(kk-1))/2;
                            end
                        end
                        del_k_arr_breaking(1) = (k_integral_arr_breaking(2) - k_integral_arr_breaking(1)) / 2;
                        del_k_arr_breaking(end) = (k_integral_arr_breaking(end) - k_integral_arr_breaking(end-1)) / 2;
                        integralx_2_breaking = sum(integralx_1_breaking .* del_k_arr_breaking);
                        integraly_2_breaking = sum(integraly_1_breaking .* del_k_arr_breaking);
                        % include breaking stress
                        tao_tx =  tao_vx + integralx_2 + integralx_2_breaking;
                        tao_ty =  tao_vy + integraly_2 + integraly_2_breaking;
                        
                        
                        tao_t = sqrt(tao_tx^2 + tao_ty^2);
                        if tao_t <= 10^-8
                            cos_taot = 1; sin_taot = 0;
                        else
                            cos_taot = tao_tx/tao_t; sin_taot = tao_ty/tao_t;
                        end
                        h_the_arr = (cos(theta_arr)*cos_taot + sin(theta_arr)*sin_taot).^2;
                        h_the_arr(abs(theta_arr) > pi/2) = 0;
                        iteration_after = C_beta * w_arr(i) / rou_w / c_arr(i)^2 * tao_t * h_the_arr;
                        beta_g_breaking(i, :) =  C_beta * w_arr_breaking(i) / rou_w / c_arr_breaking(i)^2 * tao_t * h_the_arr;
                        iteration_after(iteration_after <= 10^(-8)) = 0;
                        if max(abs(iteration_after-iteration_before)) < 10^(-8)
                            break
                        end
                        %beta_g(i, :) = (iteration_after + iteration_before) / 2;
                        beta_g(i, :) = iteration_after;  %  This one converges faster
                    end
                    tao_tx_z(i) = tao_tx;
                    tao_ty_z(i) = tao_ty;
                    iteration_n_arr(i) = iteration;
                end
                
                %% calculate the sea spray
                integral_1_lamda_k = zeros(length(k_arr),1);
                integral_1_slope= zeros(length(k_arr),1);
                integral_1_W= zeros(length(k_arr),1);
                del_k_arr = zeros(length(k_arr),1);
                if length(k_arr) > 2
                    for kk = 2 : length(k_arr)-1
                        del_k_arr(kk) = (k_arr(kk+1) - k_arr(kk-1))/2;
                    end
                end
                del_k_arr(1) = (k_arr(2) - k_arr(1)) / 2;
                del_k_arr(end) = (k_arr(end) - k_arr(end-1)) / 2;
                for i = 1 : length(k_arr)
                    fai_arr_slope =(fai_k_the_arr(:, i)'  *  k_arr(i) ^4).^0.5;
                    integral_1_slope(i) =0.25*(2* (sum(fai_arr_slope(2:end-1))*2*pi/36 + (fai_arr_slope(1)+fai_arr_slope(end))*2*pi/36/2)^0.5) ^2.5;
                    %%  Muller(2009)
                    S_Muller_k_the =fai_k_the_arr(:, i)'  *  k_arr(i) ^3;
                    S_Muller_k(i) =( sum(S_Muller_k_the(2:end-1))*2*pi/36 + (S_Muller_k_the(1)+S_Muller_k_the(end))*2*pi/36/2) ;
                    S_Muller_integral = S_Muller_k(1:i);
                    del_k_arr_s = del_k_arr(1:i);
                    b_k_s(i)= min(max( 0.01*sum( S_Muller_integral' .*  del_k_arr_s )^0.5, b_kmin),b_kmax);
                    b_k_d(i)= min(max( 0.25*sum( 2*(S_Muller_integral').^0.5 .* del_k_arr_s )^2.5, b_kmin),b_kmax);
                    b_k_d_fix(i)= min(max( 0.25*(2*sum( S_Muller_integral' .* del_k_arr_s )^0.5) ^2.5, b_kmin),b_kmax);
                    b_k(i)= W_quanzhong_after * b_k_d_fix(i) + (1-W_quanzhong_after) *b_k_d(i);
                    fai_arr_H =fai_k_the_arr(:, i)' .^0.5 *  k_arr(i) ;   %  2ak = hk = B（k）^0.5
                    integral_1_H = (sum(fai_arr_H(2:end-1))*2*pi/36 + (fai_arr_H(1)+fai_arr_H(end))*2*pi/36/2); %          %integral_1_H =2* gama_b/k_arr(i);;
                    fai_arr_lamda =beta_g(i, :).*fai_k_the_arr(:, i)' * g ^3 /  c_arr(i)^6 / w_arr(i) / b_k(i) *  k_arr(i) .* fai_arr_H ;
                    integral_1_lamda_k(i) = sum(fai_arr_lamda(2:end-1))*2*pi/36 + (fai_arr_lamda(1)+fai_arr_lamda(end))*2*pi/36/2;
                end
                del_k_arr(k_arr>10) = 0; % Kudryavtsev (2006)
                integral_2_s0 = sum( integral_1_lamda_k .* del_k_arr )  * kapa;%   fai_arr_H =fai_k_the_arr(:, i)' .^0.5 *  k_arr(i) ;
                
                %s0_Kudry = 2.5*10^(-5)*u_star^5;     % Kudryavrev 2011
                %delta_H = 9.81 * H / U10^2;
                %s0 = 2.5*10^(-5)*u_star^5;    % Kudryavrev 2011
                %s0_chen = 7.5*10^(-11)*U10^4*u_star;   %chen 2018(XY 2021)
                %s0_xy = 3.25*10^( -9 ) * 9.23 * 10^(-4) * U10 ^ 5 ;  %XY 2021
                % if delta_H >= 0.3
                %  elseif delta_H < 0.3
                %  s0_xy =min((delta_H)^(-1/2)-0.83 ,10) * s0_xy;
                %  s0 = 7.5*10^(-11)*U10^4*u_star;   % Kudryavrev 2011
                % end
                zmin = delta / k_max;
                zmax = delta / k_min;
                %% wave-induced stress
                z_arr = delta ./ k_arr;
                tao_wx_z = zeros(length(z_arr),1);
                tao_wy_z = zeros(length(z_arr),1);
                for i = 2 : length(z_arr)
                    beta_g_integral = beta_g(1:i,:);
                    k_integral_arr = k_arr(1:i);
                    w_integral_arr = w_arr(1:i);
                    integralx_1 = zeros(length(k_integral_arr),1);
                    integraly_1 = zeros(length(k_integral_arr),1);
                    for ii = 1 : length(k_integral_arr)
                        k_integral = k_integral_arr(ii);
                        w_integral = w_integral_arr(ii);
                        the_integral_arr = -pi : 2*pi/36 : pi;
                        fai_arr = fai_k_the_arr(:, ii)';
                        integralx_arr = beta_g_integral(ii,:) .* fai_arr * rou_w * w_integral * k_integral .* cos(the_integral_arr);
                        integralx_1(ii) = sum(integralx_arr(2:end-1))*2*pi/36 + (integralx_arr(1)+integralx_arr(end))*2*pi/36/2;
                        integraly_arr = beta_g_integral(ii,:) .* fai_arr * rou_w * w_integral * k_integral .* sin(the_integral_arr);
                        integraly_1(ii) = sum(integraly_arr(2:end-1))*2*pi/36 + (integraly_arr(1)+integraly_arr(end))*2*pi/36/2;
                    end
                    
                    del_k_arr = zeros(length(k_integral_arr),1);
                    if length(k_integral_arr) > 2
                        for kk = 2 : length(k_integral_arr)-1
                            del_k_arr(kk) = (k_integral_arr(kk+1) - k_integral_arr(kk-1))/2;
                        end
                    end
                    del_k_arr(1) = (k_integral_arr(2) - k_integral_arr(1)) / 2;
                    del_k_arr(end) = (k_integral_arr(end) - k_integral_arr(end-1)) / 2;
                    integralx_2 = sum(integralx_1 .* del_k_arr);
                    integraly_2 = sum(integraly_1 .* del_k_arr);
                    tao_wx_z(i) = integralx_2;
                    tao_wy_z(i) = integraly_2;
                end
                tao_w_z = (tao_wy_z.^2 + tao_wx_z.^2 ).^ 0.5;
                %% breaking stress
                for i = 2 : length(z_arr)
                    k_integral_arr_breaking = k_arr_breaking(1:i);
                    w_integral_arr_breaking = w_arr_breaking(1:i);
                    integralx_1_breaking = zeros(length(k_integral_arr_breaking),1);
                    integraly_1_breaking = zeros(length(k_integral_arr_breaking),1);
                    for ii = 1 : length(k_integral_arr_breaking)
                        k_integral_breaking = k_integral_arr_breaking(ii);
                        w_integral_breaking = w_integral_arr_breaking(ii);
                        c_integral_breaking  = w_integral_breaking/k_integral_breaking;
                        the_integral_arr = -pi : 2*pi/36 : pi;
                        lamda_arr =  lamda_breaking_the_arr(:, ii)';
                        
                        z_breaking = gama_b / k_integral_breaking;
                        z0 = 10/(exp(1)-1);
                        if iteration_u ==1
                            ub = U10 * log((z_breaking+z0)/z0) * h_the_arr.^0.5 - c_integral_breaking;
                            ub = ub .* (ub > 0) ;
                        else
                            [MIN, I ]=min(abs (z_arr_total-z_breaking));
                            ub = (u_arr(I))* h_the_arr.^0.5  - c_integral_breaking;
                            ub = ub .* (ub > 0) ;
                        end
                        integralx_arr_breaking = rou_a * 2 * gama_b * C_b * ub.^2 .* lamda_arr.* cos(the_integral_arr) ;   %Kudryavtsev（2014）
                        integraly_arr_breaking = rou_a * 2 * gama_b * C_b * ub.^2 .* lamda_arr .* sin(the_integral_arr) ;   %Kudryavtsev（2014）
                        integralx_1_breaking(ii) = sum(integralx_arr_breaking(2:end-1))*2*pi/36 + (integralx_arr_breaking(1)+integralx_arr_breaking(end))*2*pi/36/2;
                        integraly_1_breaking(ii) = sum(integraly_arr_breaking(2:end-1))*2*pi/36 + (integraly_arr_breaking(1)+integraly_arr_breaking(end))*2*pi/36/2;
                    end
                    
                    del_k_arr_breaking = zeros(length(k_integral_arr_breaking),1);
                    if length(k_integral_arr_breaking) > 2
                        for kk = 2 : length(k_integral_arr_breaking)-1
                            del_k_arr_breaking(kk) = (k_integral_arr_breaking(kk+1) - k_integral_arr_breaking(kk-1))/2;
                        end
                    end
                    del_k_arr_breaking(1) = (k_integral_arr_breaking(2) - k_integral_arr_breaking(1)) / 2;
                    del_k_arr_breaking(end) = (k_integral_arr_breaking(end) - k_integral_arr_breaking(end-1)) / 2;
                    integralx_2_breaking = sum(integralx_1_breaking .* del_k_arr_breaking);
                    integraly_2_breaking = sum(integraly_1_breaking .* del_k_arr_breaking);
                    tao_sx_z(i) = integralx_2_breaking;
                    tao_sy_z(i) = integraly_2_breaking;
                end
                tao_s_z =( tao_sy_z.^2+tao_sx_z.^2).^0.5;
                
                %% turbulent stress
                tao_tx_z2 =tao_tx_z;
                tao_ty_z2 =tao_ty_z;
                tao_t_z2 = sqrt(tao_tx_z2.^2 + tao_ty_z2.^2);
                
                taox = tao_tx(1);% Total stress equal to the turbulent stress at the highest altitude
                taoy = tao_ty(1);
                tao = sqrt(taox^2 + taoy^2);
                u_star = sqrt(tao/rou_a);% resistive flow rate
                
                % Characteristic wave height solution
                integral_1 = zeros(length(k_arr),1);
                for i = 1 : length(k_arr)
                    fai_arr = fai_k_the_arr(:, i)' * k_arr(i);
                    integral_1(i) = sum(fai_arr(2:end-1))*2*pi/36 + (fai_arr(1)+fai_arr(end))*2*pi/36/2;
                end
                
                del_k_arr = zeros(length(k_arr),1);
                if length(k_arr) > 2
                    for kk = 2 : length(k_arr)-1
                        del_k_arr(kk) = (k_arr(kk+1) - k_arr(kk-1))/2;
                    end
                end
                del_k_arr(1) = (k_arr(2) - k_arr(1)) / 2;
                del_k_arr(end) = (k_arr(end) - k_arr(end-1)) / 2;
                integral_2 = sum(integral_1 .* del_k_arr);
                Hs = 4*sqrt(integral_2);
                b = 0.01;
                
                %% to calculate sea spray use U10
                %             Hs_d (n_u10_arr) = Hs ;
                %clear u_s_before
                for iteration_s0 = 1 : 1
                    
                    u_s_before = U10;
                    s0 = (u_s_before)^3 * integral_2_s0 ;
                    
                    
                    num1 = 20;
                    num3 = 20;
                    z_arr3 = zmax : (10-zmax)/20 : 10;
                    zv = 0.1 * nu_a / sqrt(tao_v/rou_a);
                    zv = min(zv , zmin /2);
                    z_arr1 = zv : (zmin-zv)/num3 : zmin;
                    z_arr_total = [z_arr1(1:end-1), z_arr(end:-1:2), z_arr3];
                    z_plus = z_arr_total * u_star  / nu_a;
                    weight =(z_plus - 5) / 65;
                    linear_z = 2.5 * nu_a/u_star;
                    power = 0.3;
                    down_ = 5;
                    up_ = 70 ;
                    z_turbulence = (z_plus < down_).* linear_z +...
                        (z_plus <= up_ & z_plus >= down_) .* (linear_z*(1-weight.^power) + weight.^power.*z_arr_total) +  ...
                        (z_plus > up_) .* z_arr_total;
                    
                    z_arr2_tubulence = z_turbulence (num3+100:-1:num3+1);
                    z_arr1_tubulence  =  z_turbulence (1:num3+1);
                    % z_tubulence is calculated to avoid zv> delta/kmax
                    % with very low wind speed U10 < 1m/s
                    % the fix method can be refered to <Boundary-Layer Theory> 2017 by Schlichting
                    %% uplayer wind profile
                    
                    duxdz3 = (rou_a/0.4./z_arr3.*(tao/rou_a)^1.5 ...
                        - 9.81*(rou_w-rou_a)*2.0*0.4.*z_arr3*u_star*s0*log(b)/Hs.*exp(z_arr3*log(b)/Hs))...
                        * taox / tao^2;
                    duydz3 = (rou_a/0.4./z_arr3.*(tao/rou_a)^1.5 ...
                        - 9.81*(rou_w-rou_a)*2.0*0.4.*z_arr3*u_star*s0*log(b)/Hs.*exp(z_arr3*log(b)/Hs))...
                        * taoy / tao^2;
                    
                    %% midlayer wind profile
                    
                    % Fw
                    Fw_arr = zeros(1,length(k_arr));
                    for i = 1 : length(k_arr)
                        k_integral = k_arr(i);
                        fai_arr = fai_k_the_arr(:, i)';
                        integral_arr = beta_g(i,:) .* fai_arr * rou_w * 9.81 * k_integral;
                        Fw_arr(i) = sum(integral_arr(2:end-1))*2*pi/36 + (integral_arr(1)+integral_arr(end))*2*pi/36/2;
                    end
                    %Fws
                    for i = 1 : length(k_arr_breaking)
                        k_integral_breaking = k_arr_breaking(i);
                        w_integral_breaking = w_arr_breaking(i);
                        c_integral_breaking =c_arr_breaking(i);
                        z_breaking = gama_b / k_integral_breaking;
                        if iteration_u ==1
                            ub = U10 * log((z_breaking+z0)/z0) * h_the_arr.^0.5 - c_integral_breaking;
                            ub = ub .* (ub > 0) ;
                        else
                            [MIN, I ]=min(abs (z_arr_total-z_breaking));
                            ub = (u_arr(I))* h_the_arr.^0.5  - c_integral_breaking;
                            ub = ub .* (ub > 0) ;
                        end
                        lamda_arr =  lamda_breaking_the_arr(:, i)';
                        integral_arr_breaking = 2 * rou_a * gama_b * C_b * ub.^2 .* lamda_arr * c_integral_breaking  * gama_b / delta;   %Kudryavtsev（2014）
                        Fws_arr(i) = sum(integral_arr_breaking(2:end-1))*2*pi/36 + (integral_arr_breaking(1)+integral_arr_breaking(end))*2*pi/36/2;
                    end
                    
                    Fwa_arr = Fw_arr+Fws_arr;
                    
                    duxdz2 = (delta./(z_arr).^2.*Fwa_arr ...
                        + rou_a/0.4./z_arr2_tubulence.*(tao_t_z2'/rou_a).^1.5 ...
                        - 9.81*(rou_w-rou_a)*2.0*0.4.*z_arr.*sqrt((tao_t_z2'/rou_a))*s0*log(b)/Hs.*exp(z_arr*log(b)/Hs)) ...
                        .* tao_tx_z2'./(taox*tao_tx_z2'+taoy*tao_ty_z2');
                    duydz2 = (delta./(z_arr).^2.*Fwa_arr ...
                        + rou_a/0.4./z_arr2_tubulence.*(tao_t_z2'/rou_a).^1.5 ...
                        - 9.81*(rou_w-rou_a)*2.0*0.4.*z_arr.*sqrt((tao_t_z2'/rou_a))*s0*log(b)/Hs.*exp(z_arr*log(b)/Hs)) ...
                        .* tao_ty_z2'./(taox*tao_tx_z2'+taoy*tao_ty_z2');
                    
                    %% downlayer wind profile
                    duxdz1 = (rou_a/0.4./z_arr1_tubulence * (tao_v/rou_a)^1.5 ...
                        - 9.81*(rou_w-rou_a)*2.0*0.4.*z_arr1*sqrt((tao_v/rou_a))*s0*log(b)/Hs.*exp(z_arr1*log(b)/Hs)) ...
                        * tao_vx / (taox*tao_vx + taoy*tao_vy);
                    duydz1 = (rou_a/0.4./z_arr1_tubulence * (tao_v/rou_a)^1.5 ...
                        - 9.81*(rou_w-rou_a)*2.0*0.4.*z_arr1*sqrt((tao_v/rou_a))*s0*log(b)/Hs.*exp(z_arr1*log(b)/Hs)) ...
                        * tao_vy / (taox*tao_vx + taoy*tao_vy);
                    
                    %% Integration of wind speed profiles
                    duxdz_total = [duxdz1(1:end-1), duxdz2(end:-1:2), duxdz3];
                    duydz_total = [duydz1(1:end-1), duydz2(end:-1:2), duydz3];
                    ux_arr = zeros(length(z_arr_total),1);
                    uy_arr = zeros(length(z_arr_total),1);
                    for i = 2 : length(z_arr_total)
                        ux_arr(i) = ux_arr(i-1) + duxdz_total(i-1)*(z_arr_total(i)-z_arr_total(i-1));
                        uy_arr(i) = uy_arr(i-1) + duydz_total(i-1)*(z_arr_total(i)-z_arr_total(i-1));
                    end
                    %% Iterate wind speed profiles for breaking stresses
                    [m,I_Hs]=min(abs(z_arr_total -  min([Hs/2,10]))); %,
                    u_s_after = min(sqrt(ux_arr(I_Hs)^2 + uy_arr(I_Hs)^2),U10/ 2^0.33);
                    if abs(u_s_after - u_s_before) <= 10^(-2)  || Hs > 20
                        break
                    end
                    u_s_before = u_s_after;
                end
                
                u_arr = sqrt(ux_arr.^2 + uy_arr.^2);
                %% Iterate wind speed profiles for breaking stresses
                
                if iteration_u == 1
                    u_iteration_before =  U10 * log((z_arr_total+z0)/z0) ;
                end
                u_iteration_after = u_arr ;
                if sum(abs(u_iteration_before - u_iteration_after)) <= max(10^(-2),U10/100)
                    break
                end
                u_iteration_before= u_iteration_after;
                if isnan(u_iteration_before(end))
                    break
                end
            end
            delta_iteration = (U10 - u_arr(end))/u_arr(end);
            tao_vx_new = tao_vx * (1+delta_iteration);
            if abs(u_arr(end) - U10) <= max( U10/1000 ,10^(-2)) %Satisfaction of iterative accuracy to calculate wind stress
                break
            end
            if tao_vx_new <= 0
                tao_vx_new = tao_vx / 10;
            end
            tao_vx = tao_vx_new;
            tao_vy = 0;
            if isnan(tao_vx)
                break
            end
        end
        if iteration_taov == 200
            disp('Non Convergence of Iteration')
            %             quit cancel
        end
        
        u_star.^2 / U10^2
        
        Cd(n_u10_arr) = u_star.^2 / U10^2;
    end
    
    %% Cd-U10
    figure(3);
    Cd_storage_breaking(n_H_arr,:) = Cd;
    plot(u10_arr,Cd,'-','linewidth',3,'Color',RGB(fix(n_H_arr/length(H_arr)*256),:)) ; hold on
    ylabel('Cd'); xlabel('U_1_0 m/s')
    ylim([0,5]*10^-3);xlim([0,70])
    
end
RGB = othercolor('Spectral4');
Function_plot_Cd_measure(RGB)
img=gcf;
print(img,'-dpng','-r1200',['./Cd.png'])


