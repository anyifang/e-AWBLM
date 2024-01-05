% e_AWBLM to estimate wind stress

% writed by Xu and Yu (2021)
% include the water depth effect on the Cd through the TMA wave-spectrum
% refer to Moon (2004) to solve the taov iteratively
% Zhang and Yu enhanced by adding breaking-induce stress 2023 / 4 /20
% The default unit is /kg /m /s
close all
clear
addpath('../data')
addpath('../function')
load ('k_store.mat')
%% Cd measurement
%RGB = [249, 87, 56; 238, 150, 75; 244, 211, 94; 13, 59, 102]/256;
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
    for n_u10_arr = 1 : length(u10_arr)  % different wind velocity at 10-m level
        U10 = u10_arr(n_u10_arr) ;
        disp(['H=',num2str(H),' U10=',num2str(U10)])
        %     H = H_arr(n_u10_arr);
        fp_ = 3.5 * fetch_^(-0.33);
        fp = max(0.13,fp_) * 9.81 / U10; % peak frequency (Alves et al., 2003)
        kp = k_calcu(fp,H);
        kp1 = k_calcu(0.7*fp,H);
        kp2 = k_calcu(1.3*fp,H);
        kp_ = kp * U10^2 / 9.81;
        taow_max_Z = max(0.02/kp, 0.05/400) ;
       
        tao_vx = 0.001;   % initial viscous stress
        tao_vy = 0;           
        
        u_star = U10/28; % initial u*
        k_min = 0.01; %0.0049 * 9.81 / (u10/28)^2; 
        k_max = 400;
        k_arr = logspace(log10(k_min), log10(k_max), 100);
        the_arr = -pi : 2*pi/36 : pi; % wave direction
        w_arr = sqrt(k_arr* 9.81 .* tanh(k_arr * H) ); 
        c_arr = 9.81 ./ w_arr .* tanh(k_arr * H);    
        b_kmin = 0.001;   % Melville, Veron, and White (2002) 0.007
        b_kmax= 0.1;
        C_beta = 25; % empirical constant for sea spray
        C_db = 0.35; % pressure drop coefficient, Kudryavtsev（2014）
        
        delta = 0.05; % decay length scale of the breaking-wave induced stress
        rou_w = 10^3; % water density
        rou_a = 1.293;  % air density
        nu_a = 14.8*10^(-6);  % coefficient of kinematic viscosity
        g = 9.81;
      
        fai_k_the_arr_3rd = zeros(length(the_arr), length(k_arr));
        fai_k_the_arr_4th = zeros(length(the_arr), length(k_arr));
        U10up =U10;
        for i = 1 : length(the_arr)
            for j = 1 : length(k_arr)
                if abs(the_arr(i)) <= pi/2
                    fai_k_the_arr_3rd(i, j) = S_JON_k_the_sh_highequ_3rd(k_arr(j),the_arr(i),U10up,H,fetch_,C_beta,'XY2021');  
                    fai_k_the_arr_4th(i, j) = S_JON_k_the_sh_highequ_4th(k_arr(j),the_arr(i),U10up,H,fetch_,C_beta,'XY2021');   
                else
                    fai_k_the_arr_3rd(i, j) = 0;
                    fai_k_the_arr_4th(i, j) = 0;
                end
            end
        end
        
        fai_k_the_arr = (fai_k_the_arr_4th+fai_k_the_arr_3rd)/2;   
        W_quanzhong_before = 0.5;
        for diedai_quanzhong = 1 : 200
           %% to calculate TMA wave spectrum in arbitrary water depth
            % this part is calculted based on XY(2021)
            jifen_1 = zeros(length(k_arr),1); % to integral in theta
            for i = 1 : length(k_arr)
                fai_arr = fai_k_the_arr(:, i)' * k_arr(i); 
                jifen_1(i) = sum(fai_arr(2:end-1))*2*pi/36 + (fai_arr(1)+fai_arr(end))*2*pi/36/2;
            end
            del_k_arr = zeros(length(k_arr),1);
            if length(k_arr) > 2
                for kk = 2 : length(k_arr)-1
                    del_k_arr(kk) = (k_arr(kk+1) - k_arr(kk-1))/2;
                end
            end
            del_k_arr(1) = (k_arr(2) - k_arr(1)) / 2;
            del_k_arr(end) = (k_arr(end) - k_arr(end-1)) / 2; % to integral in k
            jifen_2 = sum(jifen_1 .* del_k_arr);
            Hs_o = 4*sqrt(jifen_2); 
 
            f_arr = sqrt(9.81/4/pi^2 * k_arr .* tanh(k_arr*H));  
            jifen_1 = zeros(length(f_arr),1); 
            for i = 1 : length(f_arr)
                fai_arr_f = fai_k_the_arr(:, i)' * k_arr(i) * 8*pi^2*f_arr(i) / (9.81*tanh(k_arr(i)*H) + 9.81*k_arr(i)*H*(sech(k_arr(i)*H))^2) * f_arr(i); 
                jifen_1(i) = sum(fai_arr_f(2:end-1))*2*pi/36 + (fai_arr_f(1)+fai_arr_f(end))*2*pi/36/2;
            end
            del_f_arr = zeros(length(f_arr),1);
            if length(f_arr) > 2
                for ff = 2 : length(f_arr)-1
                    del_f_arr(ff) = (f_arr(ff+1) - f_arr(ff-1))/2;
                end
            end
            del_f_arr(1) = (f_arr(2) - f_arr(1)) / 2;
            del_f_arr(end) = (f_arr(end) - f_arr(end-1)) / 2;
            jifen_2 = sum(jifen_1 .* del_f_arr);
            Tm01 = (Hs_o/4)^2 / jifen_2;
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
       % Epslion_b_mean =   [ 0.0625    0.0571    0.0501    0.0448    0.0339];  % 10 20 50 100 500
         Epslion_b_mean =   [ 0.0580    0.0456   0.0339];  % 15 80 500
        %Epslion_b = Epslion_b_mean(n_H_arr)/Epslion_b_mean(end)*0.3;
        Epslion_b = gama_b_U_H_quick(U10,H,k_store)/gama_b_U_H_quick(U10,500,k_store)*0.3;
        %Epslion_b = 16*Epslion_b_U_H(U10,H)
        Epslion_b = min(0.6,Epslion_b);
        Epslion_b = max(0.2,Epslion_b);
        %Epslion_b = 0.35
        Weigth =W_quanzhong_after;  % 10m 100 1000m
        k_arr_AFS = Epslion_b /delta * k_arr;
        w_arr_AFS = sqrt(k_arr_AFS* 9.81 .* tanh(k_arr_AFS * H) ); 
        c_arr_AFS = 9.81 ./ w_arr_AFS .* tanh(k_arr_AFS * H); 
        fai_k_the_arr_3rd_AFS = zeros(length(the_arr), length(k_arr));
        fai_k_the_arr_4th_AFS = zeros(length(the_arr), length(k_arr));
        for i = 1 : length(the_arr)
            for j = 1 : length(k_arr)
                if abs(the_arr(i)) <= pi/2
                    fai_k_the_arr_3rd_AFS (i, j) = S_JON_k_the_sh_highequ_3rd(k_arr_AFS (j),the_arr(i),U10up,H,fetch_,C_beta,'XY2021');
                    fai_k_the_arr_4th_AFS (i, j) = S_JON_k_the_sh_highequ_4th(k_arr_AFS (j),the_arr(i),U10up,H,fetch_,C_beta,'XY2021');
                else
                    fai_k_the_arr_3rd_AFS(i, j) = 0;
                    fai_k_the_arr_4th_AFS(i, j) = 0;
                end
            end
        end
        fai_k_the_arr_AFS = fai_k_the_arr_4th_AFS * W_quanzhong_after + fai_k_the_arr_3rd_AFS * (1-W_quanzhong_after);   %Wave spectrum at arbitrary water depth
        
       %%  The main iterative part of the procedure such that U(10) = U10.
        for diedai_taov = 1 : 200
                 % wave growth rate
            beta_g = zeros(length(k_arr), length(the_arr));
            beta_g_AFS = zeros(length(k_arr), length(the_arr)); 
            tao_tx_z = zeros(length(k_arr),1); 
            tao_ty_z = zeros(length(k_arr),1);
           
            tao_v = sqrt(tao_vx^2 + tao_vy^2); 

            if tao_v <= 10^-8
                cos_taov = 1; sin_taov = 0;
            else
                cos_taov = tao_vx/tao_v; sin_taov = tao_vy/tao_v;
            end
            h_the_arr = (cos(the_arr)*cos_taov + sin(the_arr)*sin_taov).^2; 
            h_the_arr(abs(the_arr) > pi/2) = 0;
            beta_g(end, :) = C_beta * w_arr(end) / rou_w / c_arr(end)^2 * tao_v * h_the_arr;                    
            beta_g_AFS(end, :) = C_beta * w_arr_AFS(end) / rou_w / c_arr_AFS(end)^2 * tao_v * h_the_arr;  
            lamda_AFS_the_arr(:,length(k_arr)) =zeros(length(the_arr),1) ;
            jifen_1_lamda_k_AFS(length(k_arr)) = 0;
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
            diedai_n_arr = zeros(length(k_arr)-1,1); %  Store the number of iterations to see if there is any non-convergence
            for diedai_u = 1 : 100    %iterative part to calculate tao_b
                  % Next iteratively solve for beta_g(k) and taot on each k.
                % Since tao and beta are functions of each other, solving for both requires iteration
                for i = length(k_arr)-1 : -1 : 1
                    for diedai = 1 : 200
                        diedai_before = beta_g(i,:);
                        % calculate lamda（k）based on  Sin = Sdb
                        if k_arr_AFS(i)>20*pi  % Kudryavtsev（2001） tao_b = 0 with  k >20pi
                            lamda_AFS_the_arr(:,i) =zeros(length(the_arr),1) ;
                            jifen_1_lamda_k_AFS(i) = 0;
                        else
                            fai_arr_slope_AFS =fai_k_the_arr_AFS(:, i)'  *  k_arr_AFS(i) ^4;
                            jifen_1_slope_AFS(i) =2*( sum(fai_arr_slope_AFS(2:end-1))*2*pi/36 + (fai_arr_slope_AFS(1)+fai_arr_slope_AFS(end))*2*pi/36/2) ^0.5;
                            %  Muller(2009)
                            del_k_arr_AFS = zeros(length(k_arr_AFS),1);
                            if length(k_arr_AFS) > 2
                                for kk = 2 : length(k_arr_AFS)-1
                                    del_k_arr_AFS(kk) = (k_arr_AFS(kk+1) - k_arr_AFS(kk-1))/2;
                                end
                            end
                            del_k_arr_AFS(1) = (k_arr_AFS(2) - k_arr_AFS(1)) / 2;
                            del_k_arr_AFS(end) = (k_arr_AFS(end) - k_arr_AFS(end-1)) / 2;
                            for num_s = 1:length(k_arr)  
                                S_Muller_k_the_AFS =fai_k_the_arr_AFS(:, num_s)'  *  k_arr_AFS(num_s) ^3;
                                S_Muller_k_AFS(num_s) =2*( sum(S_Muller_k_the_AFS(2:end-1))*2*pi/36 + (S_Muller_k_the_AFS(1)+S_Muller_k_the_AFS(end))*2*pi/36/2) ^0.5;
                            end
                            S_Muller_k_AFS_jifen = S_Muller_k_AFS(1:i+1);  
                           del_k_arr_s_AFS = del_k_arr_AFS(1:i+1);
                            W_quanzhong_after = 1;
                            b_k_AFS_s(i)=min(max( 0.01*sum( S_Muller_k_AFS_jifen' .* del_k_arr_s_AFS )^0.5, b_kmin),b_kmax);   % spilling
                            b_k_AFS_d(i)=min(max( 0.25*sum( S_Muller_k_AFS_jifen' .* del_k_arr_s_AFS )^2.5, b_kmin),b_kmax);   % plunging
                            b_k_AFS(i)= W_quanzhong_after * b_k_AFS_d(i) + (1-W_quanzhong_after) *b_k_AFS_s(i);
  
                            lamda_AFS_the_arr(:,i) =beta_g_AFS(i, :).*fai_k_the_arr_AFS(:, i)' * g ^3 /  c_arr_AFS(i)^6 / w_arr_AFS(i) / b_k_AFS(i) ; 
                            jifen_1_lamda_k_AFS(i) = sum(lamda_AFS_the_arr(2:end-1,i))*2*pi/36 + (lamda_AFS_the_arr(1,i)+lamda_AFS_the_arr(end,i))*2*pi/36/2;
                        end
                        
                        beta_g_jifen = beta_g(i:end,:); 
                        k_jifen_arr = k_arr(i:end); 
                        w_jifen_arr = w_arr(i:end); 
                        k_jifen_arr_AFS =  k_arr_AFS (i:end);  
                        w_jifen_arr_AFS = w_arr_AFS(i:end); 

                        jifenx_1 = zeros(length(k_jifen_arr),1); 
                        jifeny_1 = zeros(length(k_jifen_arr),1); 
                        jifenx_1_AFS = zeros(length(k_jifen_arr),1); 
                        jifeny_1_AFS = zeros(length(k_jifen_arr),1); 
                        for ii = 1 : length(k_jifen_arr)
                            k_jifen = k_jifen_arr(ii);
                            w_jifen = w_jifen_arr(ii);
                            k_jifen_AFS = k_jifen_arr_AFS(ii);
                            w_jifen_AFS = w_jifen_arr_AFS(ii);
                            c_jifen_AFS  = w_jifen_AFS/k_jifen_AFS;
                            the_jifen_arr = -pi : 2*pi/36 : pi; 
                            fai_arr = fai_k_the_arr(:, i+ii-1)'; 
                            jifenx_arr = beta_g_jifen(ii,:) .* fai_arr * rou_w * w_jifen * k_jifen .* cos(the_jifen_arr);
                            jifenx_1(ii) = sum(jifenx_arr(2:end-1))*2*pi/36 + (jifenx_arr(1)+jifenx_arr(end))*2*pi/36/2;
                            jifeny_arr = beta_g_jifen(ii,:) .* fai_arr * rou_w * w_jifen * k_jifen .* sin(the_jifen_arr);
                            jifeny_1(ii) = sum(jifeny_arr(2:end-1))*2*pi/36 + (jifeny_arr(1)+jifeny_arr(end))*2*pi/36/2;
                            lamda_arr =  lamda_AFS_the_arr(:, i+ii-1)';
                          
                            z_AFS = Epslion_b / k_jifen_AFS;
                            z0 = 10/(exp(1)-1);
                            if diedai_u == 1
                                ub = U10 * log((z_AFS+z0)/z0) * h_the_arr.^0.5 - c_jifen_AFS;
                                ub = ub .* (ub > 0) ;
                            else
                                [MIN, I ]=min(abs (z_arr_total-z_AFS));
                                ub = (u_arr(I))* h_the_arr.^0.5  - c_jifen_AFS;
                                ub = ub .* (ub > 0) ;
                            end
                            jifenx_arr_AFS = rou_a * 2 * Epslion_b * C_db * ub.^2 .* lamda_arr .* cos(the_jifen_arr) ;   %Kudryavtsev（2014）
                            jifeny_arr_AFS = rou_a * 2 * Epslion_b * C_db * ub.^2 .* lamda_arr .* sin(the_jifen_arr) ;   %Kudryavtsev（2014）
                            jifenx_1_AFS(ii) = sum(jifenx_arr_AFS(2:end-1))*2*pi/36 + (jifenx_arr_AFS(1)+jifenx_arr_AFS(end))*2*pi/36/2;
                            jifeny_1_AFS(ii) = sum(jifeny_arr_AFS(2:end-1))*2*pi/36 + (jifeny_arr_AFS(1)+jifeny_arr_AFS(end))*2*pi/36/2;
                        end
                        del_k_arr = zeros(length(k_jifen_arr),1);
                        if length(k_jifen_arr) > 2
                            for kk = 2 : length(k_jifen_arr)-1
                                del_k_arr(kk) = (k_jifen_arr(kk+1) - k_jifen_arr(kk-1))/2;
                            end
                        end
                        del_k_arr(1) = (k_jifen_arr(2) - k_jifen_arr(1)) / 2;
                        del_k_arr(end) = (k_jifen_arr(end) - k_jifen_arr(end-1)) / 2;
                        jifenx_2 = sum(jifenx_1 .* del_k_arr);
                        jifeny_2 = sum(jifeny_1 .* del_k_arr);
                        
                        del_k_arr_AFS = zeros(length(k_jifen_arr_AFS),1);
                        if length(k_jifen_arr_AFS) > 2
                            for kk = 2 : length(k_jifen_arr_AFS)-1
                                del_k_arr_AFS(kk) = (k_jifen_arr_AFS(kk+1) - k_jifen_arr_AFS(kk-1))/2;
                            end
                        end
                        del_k_arr_AFS(1) = (k_jifen_arr_AFS(2) - k_jifen_arr_AFS(1)) / 2;
                        del_k_arr_AFS(end) = (k_jifen_arr_AFS(end) - k_jifen_arr_AFS(end-1)) / 2;
                        jifenx_2_AFS = sum(jifenx_1_AFS .* del_k_arr_AFS);
                        jifeny_2_AFS = sum(jifeny_1_AFS .* del_k_arr_AFS);
                         % include breaking stress
                        tao_tx =  tao_vx + jifenx_2 + jifenx_2_AFS;
                        tao_ty =  tao_vy + jifeny_2 + jifeny_2_AFS;
     
                        tao_t = sqrt(tao_tx^2 + tao_ty^2);
                        if tao_t <= 10^-8
                            cos_taot = 1; sin_taot = 0;
                        else
                            cos_taot = tao_tx/tao_t; sin_taot = tao_ty/tao_t;
                        end
                        h_the_arr = (cos(the_arr)*cos_taot + sin(the_arr)*sin_taot).^2;
                        h_the_arr(abs(the_arr) > pi/2) = 0;
                        diedai_after = C_beta * w_arr(i) / rou_w / c_arr(i)^2 * tao_t * h_the_arr; 
                        beta_g_AFS(i, :) =  C_beta * w_arr_AFS(i) / rou_w / c_arr_AFS(i)^2 * tao_t * h_the_arr;
                        diedai_after(diedai_after <= 10^(-8)) = 0;
                        if max(abs(diedai_after-diedai_before)) < 10^(-8)
                            break
                        end
                        beta_g(i, :) = diedai_after;  %  This one converges faster
                    end
                    tao_tx_z(i) = tao_tx; 
                    tao_ty_z(i) = tao_ty;
                    diedai_n_arr(i) = diedai;
                end
                
               %% calculate the sea spray
                jifen_1_lamda_k = zeros(length(k_arr),1);
                jifen_1_slope= zeros(length(k_arr),1);
                jifen_1_W= zeros(length(k_arr),1);
                del_k_arr = zeros(length(k_arr),1);
                if length(k_arr) > 2
                    for kk = 2 : length(k_arr)-1
                        del_k_arr(kk) = (k_arr(kk+1) - k_arr(kk-1))/2;
                    end
                end
                del_k_arr(1) = (k_arr(2) - k_arr(1)) / 2;
                del_k_arr(end) = (k_arr(end) - k_arr(end-1)) / 2;
                for i = 1 : length(k_arr)
                    fai_arr_slope =fai_k_the_arr(:, i)'  *  k_arr(i) ^4;
                    jifen_1_slope(i) =2*( sum(fai_arr_slope(2:end-1))*2*pi/36 + (fai_arr_slope(1)+fai_arr_slope(end))*2*pi/36/2) ^0.5;
                    %%  Muller(2009)
                    S_Muller_k_the =fai_k_the_arr(:, i)'  *  k_arr(i) ^3;
                    S_Muller_k(i) =2*( sum(S_Muller_k_the(2:end-1))*2*pi/36 + (S_Muller_k_the(1)+S_Muller_k_the(end))*2*pi/36/2) ^0.5;     %这个是算一次积分一次
                    S_Muller_jifen = S_Muller_k(1:i);
                    del_k_arr_s = del_k_arr(1:i);
                    b_k_s(i)= min(max( 0.01*sum( S_Muller_jifen' .*  del_k_arr_s )^0.5, b_kmin),b_kmax);  %浅水下b
                    b_k_d(i)=min(max( 0.25*sum( S_Muller_jifen' .* del_k_arr_s )^2.5, b_kmin),b_kmax);    %深水下b
                    b_k(i)= W_quanzhong_after * b_k_d(i) + (1-W_quanzhong_after) *b_k_s(i);
                    fai_arr_H =fai_k_the_arr(:, i)' .^0.5 *  k_arr(i) ;   % 饱和谱波陡得到波高
                    jifen_1_H = sum(fai_arr_H(2:end-1))*2*pi/36 + (fai_arr_H(1)+fai_arr_H(end))*2*pi/36/2;
                    fai_arr_lamda =beta_g(i, :).*fai_k_the_arr(:, i)' * g ^3 /  c_arr(i)^6 / w_arr(i) / b_k(i) *  k_arr(i) * jifen_1_H ;
                    %fai_arr_lamda =beta_g(i, :).*fai_k_the_arr(:, i)' * g ^3 /  c_arr(i)^6 / w_arr(i) *  k_arr(i)  / b_k(i) ;
                    jifen_1_lamda_k(i) = sum(fai_arr_lamda(2:end-1))*2*pi/36 + (fai_arr_lamda(1)+fai_arr_lamda(end))*2*pi/36/2;
                end
                
                
                kp1 = k_calcu(0.5*fp,H); kp2 = k_calcu(1.5*fp,H);
                del_k_arr_W = del_k_arr;
                del_k_arr_W(c_arr<2) = 0;
                
                del_k_arr(k_arr>10) = 0;
                b_slope(n_u10_arr) =max( 0.01*sum( jifen_1_slope .* del_k_arr )^0.5, 0.01);
                W(n_u10_arr)=sum( jifen_1_W .* del_k_arr_W );
                jifen_2_s0 = sum( jifen_1_lamda_k .* del_k_arr )  * 8 *10 ^(-6);%b=bdeep
                %  s0_Kudry = 2.5*10^(-5)*u_star^5;     % Kudryavrev 2011
     
                
               %% wave-induced stress
                zmin = delta / k_max;
                zmax = delta / k_min;
                z_arr = delta ./ k_arr; 
                tao_wx_z = zeros(length(z_arr),1); 
                tao_wy_z = zeros(length(z_arr),1);
                for i = 2 : length(z_arr) 
                    beta_g_jifen = beta_g(1:i,:); 
                    k_jifen_arr = k_arr(1:i); 
                    w_jifen_arr = w_arr(1:i); 
                    jifenx_1 = zeros(length(k_jifen_arr),1); 
                    jifeny_1 = zeros(length(k_jifen_arr),1); 
                    for ii = 1 : length(k_jifen_arr)
                        k_jifen = k_jifen_arr(ii);
                        w_jifen = w_jifen_arr(ii);
                        the_jifen_arr = -pi : 2*pi/36 : pi; 

                        fai_arr = fai_k_the_arr(:, ii)';
                        jifenx_arr = beta_g_jifen(ii,:) .* fai_arr * rou_w * w_jifen * k_jifen .* cos(the_jifen_arr);
                        jifenx_1(ii) = sum(jifenx_arr(2:end-1))*2*pi/36 + (jifenx_arr(1)+jifenx_arr(end))*2*pi/36/2;
                        jifeny_arr = beta_g_jifen(ii,:) .* fai_arr * rou_w * w_jifen * k_jifen .* sin(the_jifen_arr);
                        jifeny_1(ii) = sum(jifeny_arr(2:end-1))*2*pi/36 + (jifeny_arr(1)+jifeny_arr(end))*2*pi/36/2;
                    end
                    del_k_arr = zeros(length(k_jifen_arr),1);
                    if length(k_jifen_arr) > 2
                        for kk = 2 : length(k_jifen_arr)-1
                            del_k_arr(kk) = (k_jifen_arr(kk+1) - k_jifen_arr(kk-1))/2;
                        end
                    end
                    del_k_arr(1) = (k_jifen_arr(2) - k_jifen_arr(1)) / 2;
                    del_k_arr(end) = (k_jifen_arr(end) - k_jifen_arr(end-1)) / 2;
                    jifenx_2 = sum(jifenx_1 .* del_k_arr);
                    jifeny_2 = sum(jifeny_1 .* del_k_arr);
                    tao_wx_z(i) = jifenx_2;
                    tao_wy_z(i) = jifeny_2;
                end
                tao_w_z = (tao_wy_z.^2 + tao_wx_z.^2 ).^ 0.5;
                 %% breaking stress
                for i = 2 : length(z_arr) 
                    k_jifen_arr_AFS = k_arr_AFS(1:i);
                    w_jifen_arr_AFS = w_arr_AFS(1:i); 
                    jifenx_1_AFS = zeros(length(k_jifen_arr_AFS),1); % theta积分结束后，k积分前
                    jifeny_1_AFS = zeros(length(k_jifen_arr_AFS),1); 
                    for ii = 1 : length(k_jifen_arr_AFS)
                        k_jifen_AFS = k_jifen_arr_AFS(ii);
                        w_jifen_AFS = w_jifen_arr_AFS(ii);
                        c_jifen_AFS  = w_jifen_AFS/k_jifen_AFS;
                        the_jifen_arr = -pi : 2*pi/36 : pi; 
                        lamda_arr =  lamda_AFS_the_arr(:, ii)';
                        z_AFS = Epslion_b / k_jifen_AFS;
                        z0 = 10/(exp(1)-1);
                        if diedai_u ==1
                            ub = U10 * log((z_AFS+z0)/z0) * h_the_arr.^0.5 - c_jifen_AFS;
                            ub = ub .* (ub > 0) ;
                        else
                            [MIN, I ]=min(abs (z_arr_total-z_AFS));
                            ub = (u_arr(I))* h_the_arr.^0.5  - c_jifen_AFS;
                            ub = ub .* (ub > 0) ;
                        end
                        jifenx_arr_AFS = rou_a * 2 * Epslion_b * C_db * ub.^2 .* lamda_arr.* cos(the_jifen_arr) ;   %Kudryavtsev（2014）
                        jifeny_arr_AFS = rou_a * 2 * Epslion_b * C_db * ub.^2 .* lamda_arr .* sin(the_jifen_arr) ;   %Kudryavtsev（2014）
                        jifenx_1_AFS(ii) = sum(jifenx_arr_AFS(2:end-1))*2*pi/36 + (jifenx_arr_AFS(1)+jifenx_arr_AFS(end))*2*pi/36/2;
                        jifeny_1_AFS(ii) = sum(jifeny_arr_AFS(2:end-1))*2*pi/36 + (jifeny_arr_AFS(1)+jifeny_arr_AFS(end))*2*pi/36/2;
                    end
                    del_k_arr_AFS = zeros(length(k_jifen_arr_AFS),1);
                    if length(k_jifen_arr_AFS) > 2
                        for kk = 2 : length(k_jifen_arr_AFS)-1
                            del_k_arr_AFS(kk) = (k_jifen_arr_AFS(kk+1) - k_jifen_arr_AFS(kk-1))/2;
                        end
                    end
                    del_k_arr_AFS(1) = (k_jifen_arr_AFS(2) - k_jifen_arr_AFS(1)) / 2;
                    del_k_arr_AFS(end) = (k_jifen_arr_AFS(end) - k_jifen_arr_AFS(end-1)) / 2;
                    jifenx_2_AFS = sum(jifenx_1_AFS .* del_k_arr_AFS);
                    jifeny_2_AFS = sum(jifeny_1_AFS .* del_k_arr_AFS);
                    tao_sx_z(i) = jifenx_2_AFS;
                    tao_sy_z(i) = jifeny_2_AFS;
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
                jifen_1 = zeros(length(k_arr),1); 
                for i = 1 : length(k_arr)
                    fai_arr = fai_k_the_arr(:, i)' * k_arr(i); 
                    jifen_1(i) = sum(fai_arr(2:end-1))*2*pi/36 + (fai_arr(1)+fai_arr(end))*2*pi/36/2;
                end

                del_k_arr = zeros(length(k_arr),1);
                if length(k_arr) > 2
                    for kk = 2 : length(k_arr)-1
                        del_k_arr(kk) = (k_arr(kk+1) - k_arr(kk-1))/2;
                    end
                end
                del_k_arr(1) = (k_arr(2) - k_arr(1)) / 2;
                del_k_arr(end) = (k_arr(end) - k_arr(end-1)) / 2;
                jifen_2 = sum(jifen_1 .* del_k_arr);
                Hs = 4*sqrt(jifen_2); 
                b = 0.01; 
                
             %% to calculate sea spray use U10
                for diedai_s0 = 1 : 1
                 u_s_before = U10;
                    s0 = (u_s_before)^3 * jifen_2_s0 ;
                    num1 = 20;
                    num3 = 20;
                    z_arr3 = zmax : (10-zmax)/20 : 10; 
                    zv = 0.1 * nu_a / sqrt(tao_v/rou_a); 
                    zv = min(zv , zmin /2);
                    z_arr1 = zv : (zmin-zv)/num3 : zmin; 
                    %             if (zv>=zmin)
                    %                 z_arr1 = z_arr1 *0;
                    %             end
                    
                    z_arr_total = [z_arr1(1:end-1), z_arr(end:-1:2), z_arr3];
                    z_plus = z_arr_total * u_star  / nu_a;
                    %z_arr_total(z_arr_total<= zv) = 0;
                    weight =(z_plus - 5) / 65;
                    linear_z = 2.5 * nu_a/u_star;%max(z_arr_total((z_plus < 5)));   
                    power = 0.3;
                    down_ = 5;
                    up_ = 70 ;
                    z_turbulence = (z_plus < down_).* linear_z +...
                        (z_plus <= up_ & z_plus >= down_) .* (linear_z*(1-weight.^power) + weight.^power.*z_arr_total) +  ...
                        (z_plus > up_) .* z_arr_total;
                    
                    z_arr2_tubulence = z_turbulence (num3+100:-1:num3+1);
                    z_arr1_tubulence  =  z_turbulence (1:num3+1);
    
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
                        k_jifen = k_arr(i);
        
                        fai_arr = fai_k_the_arr(:, i)'; 
                        jifen_arr = beta_g(i,:) .* fai_arr * rou_w * 9.81 * k_jifen;
                        Fw_arr(i) = sum(jifen_arr(2:end-1))*2*pi/36 + (jifen_arr(1)+jifen_arr(end))*2*pi/36/2;
                    end
                   %Fws
                    for i = 1 : length(k_arr_AFS)
                        k_jifen_AFS = k_arr_AFS(i);
                        w_jifen_AFS = w_arr_AFS(i);
                        c_jifen_AFS =c_arr_AFS(i);
                        z_AFS = Epslion_b / k_jifen_AFS;
                        if diedai_u ==1
                            ub = U10 * log((z_AFS+z0)/z0) * h_the_arr.^0.5 - c_jifen_AFS;
                            ub = ub .* (ub > 0) ;
                        else
                            [MIN, I ]=min(abs (z_arr_total-z_AFS));
                            ub = (u_arr(I))* h_the_arr.^0.5  - c_jifen_AFS;
                            ub = ub .* (ub > 0) ;
                        end
                        lamda_arr =  lamda_AFS_the_arr(:, i)'; 
                        jifen_arr_AFS = 2 * rou_a * Epslion_b * C_db * ub.^2 .* lamda_arr * c_jifen_AFS  * Epslion_b / delta;   %Kudryavtsev（2014）
                        Fws_arr(i) = sum(jifen_arr_AFS(2:end-1))*2*pi/36 + (jifen_arr_AFS(1)+jifen_arr_AFS(end))*2*pi/36/2;
                    end
                    
                    Fw_arr = Fw_arr+Fws_arr; 
    
                    duxdz2 = (delta./(z_arr).^2.*Fw_arr ...
                        + rou_a/0.4./z_arr2_tubulence.*(tao_t_z2'/rou_a).^1.5 ...
                        - 9.81*(rou_w-rou_a)*2.0*0.4.*z_arr.*sqrt((tao_t_z2'/rou_a))*s0*log(b)/Hs.*exp(z_arr*log(b)/Hs)) ...
                        .* tao_tx_z2'./(taox*tao_tx_z2'+taoy*tao_ty_z2');
                    duydz2 = (delta./(z_arr).^2.*Fw_arr ...
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
                    u_s_after = sqrt(ux_arr(I_Hs)^2 + uy_arr(I_Hs)^2);
                    if abs(u_s_after - u_s_before) <= 10^(-2)  || Hs > 20 
                        break
                    end
                    u_s_before = u_s_after;
                end
                
                u_arr = sqrt(ux_arr.^2 + uy_arr.^2); 
                  %% Iterate wind speed profiles for breaking stresses
                if diedai_u == 1
                    u_diedai_before =  U10 * log((z_arr_total+z0)/z0) ;
                end
                u_diedai_after = u_arr ;
                if sum(abs(u_diedai_before - u_diedai_after)) <= max(10^(-2),U10/100)  
                    break
                end
                u_diedai_before= u_diedai_after;
                if isnan(u_diedai_before(end))
                    break
                end
            end
            delta_diedai = (U10 - u_arr(end))/u_arr(end);
          
            tao_vx_new = tao_vx * (1+delta_diedai);
            if abs(u_arr(end) - U10) <= max( U10/1000 ,10^(-2))%Satisfaction of iterative accuracy to calculate wind stress
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
        if diedai_taov == 200
            disp('Non Convergence of Iteration')
            %             quit cancel
        end
        a =4;
        
        delta_H = 9.81 * H / U10^2; 
        % s0 = 2.5*10^(-5)*u_star^5; % Kudryavrev 2011
        %s0_xy = 7.5*10^(-11)*U10^4*u_star;     %chen 2018(XY 2021)
        s0_xy = 3.25*10^( -9 ) * 9.23 * 10^(-4) * U10 ^ 5 ; %XY 2021
        if delta_H >= 0.3
        elseif delta_H < 0.3
            s0_xy =min((delta_H)^(-1/2)-0.83 ,10) * s0_xy;                               
        end

        s0_arr(n_u10_arr) =s0;
        s0_cpmpare(n_u10_arr) =s0_xy ;
        jifen_2_arr(n_u10_arr) = jifen_2 ;
        u_star.^2 / U10^2
        Cd(n_u10_arr) = u_star.^2 / U10^2;

    end


    %% Cd-U10
    figure(3);
    Cd_storage_AFS(n_H_arr,:) = Cd;
    plot(u10_arr,Cd,'-','linewidth',3,'Color',RGB(fix(n_H_arr/length(H_arr)*256),:)) ; hold on
    ylabel('Cd'); xlabel('U_1_0 m/s')
    ylim([0,5]*10^-3);xlim([0,70])

end



