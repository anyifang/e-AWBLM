% Write by zhang 2024 
% To estimate gama in the arbitrary water depth
function gama_b = gama_b_U_H_quick(U,H,k_store)
fetch_ = 10^5;
U10 = U;
fp_ = 3.5 * fetch_^(-0.33);

fp = max(0.13,fp_) * 9.81 / U10;
kp = k_calcu_store(fp,H,k_store);

kp = k_calcu_store(fp,H,k_store);
kp1 =k_calcu_store(0.5* fp,H,k_store);
 kp2 =max(k_calcu_store(10* fp,H,k_store));
 kp2 / kp;
k_min = kp1;
k_max = kp2;

delta =50;
k_arr = linspace(k_min, k_max,delta);


the_arr = -pi : 2*pi/36 : pi;
C_beta = 25; 


fai_k_the_arr_3rd = zeros(length(the_arr), length(k_arr));
fai_k_the_arr_4th = zeros(length(the_arr), length(k_arr));
fai_k_the_arr_deep= zeros(length(the_arr), length(k_arr));

U10up =U10;
for i = 1 : length(the_arr)
    for j = 1 : length(k_arr)
        if abs(the_arr(i)) <= pi/2
            fai_k_the_arr_3rd(i, j) = S_JON_k_the_sh_highequ_3rd_quick(k_arr(j),the_arr(i),U10up,H,fetch_,C_beta,'XY2021',k_store);  %有限水深 TK2016
            fai_k_the_arr_4th(i, j) = S_JON_k_the_sh_highequ_4th_quick(k_arr(j),the_arr(i),U10up,H,fetch_,C_beta,'XY2021',k_store);    %深水下 XY2021     
        else
            fai_k_the_arr_3rd(i, j) = 0;
            fai_k_the_arr_4th(i, j) = 0;
            fai_k_the_arr_deep(i, j) = 0;
        end
    end
end

fai_k_the_arr = (fai_k_the_arr_4th+fai_k_the_arr_3rd)/2;   
W_quanzhong_before = 0.5;
for diedai_quanzhong = 1 : 200

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
    Hs_o = 4*sqrt(jifen_2); 
    f_arr = sqrt(9.81/4/pi^2 * k_arr .* tanh(k_arr*H));  
    jifen_1 = zeros(length(f_arr),1);
    for i = 1 : length(f_arr)
        fai_arr_f = fai_k_the_arr(:, i)' * k_arr(i) * 8*pi^2*f_arr(i) / (9.81*tanh(k_arr(i)*H) + 9.81*k_arr(i)*H*(sech(k_arr(i)*H))^2) * f_arr(i); % 对应k_jifen波数的所有方向上的波谱值
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
del_k_arr = del_k_arr.*((kp1<=k_arr) &( kp2>=k_arr))';
jifen_2 = sum(jifen_1 .* del_k_arr);
gama_b= 2*jifen_2^0.5*kp;
