load('q4_v.mat');
k_max = 1;
for i = 1:201
    v = result_v(:, i);
    v0 = v(1, 1);
    v_max = max(v);
    k = v_max / v0;
    if k > k_max
        k_max = k;
    end
end

V_upper = 2; % 各节点全时刻最大允许速度 (m/s)
V_final = V_upper / k_max; % 龙头最大允许速度 (m/s)
k_max
V_final