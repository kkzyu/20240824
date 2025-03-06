function result_xy = generating_position(start_t, end_t, h, L)
    % 初始化参数
    d1 = 2.86; % 第一个和第二个把手的距离
    d2 = 1.65; % 第i个把手与第i+1个把手的距离
    a = L / (2 * pi); % 螺线参数
    v = 1; % 速度
    theta_0 = 32 * pi; % 初始角度    
    
    if start_t == end_t
        t_lst = [start_t];
        t_max = 1;
    else
        % 生成时间序列，步长为 h
        t_lst = start_t:h:end_t;
        t_max = length(t_lst); % 时间步数
    end

    % 初始化螺线长度函数
    l_theta = @(theta) a / 2 * (theta .* sqrt(1 + theta.^2) + log(theta + sqrt(1 + theta.^2))); 
    
    % 初始螺线长度
    l_theta_0 = l_theta(theta_0);

    % 初始化位置矩阵，224个把手的 x 和 y 坐标
    result_xy = zeros(448, t_max); % 每一列存储某时间步所有把手的位置

    % 循环遍历每一个时间步
    for t_idx = 1:t_max
        current_time = t_lst(t_idx); % 获取当前时间
        
        % 计算 t 时刻螺线长度
        l_t = l_theta_0 - v * current_time;

        % 防止螺线长度出现负数
        if l_t < 0
            l_t = 0;
        end

        % 利用 fzero 解 t 时刻的 theta 值
        theta_t = fzero(@(theta) l_theta(theta) - l_t, theta_0);
        
        % 计算龙头的极坐标
        x_t = [a * theta_t, theta_t]; % 龙头的极坐标
        x_cartesian = [x_t(1) * cos(x_t(2)), x_t(1) * sin(x_t(2))]; % 转换为直角坐标

        % 记录龙头直角坐标到 result_xy
        result_xy(1, t_idx) = x_cartesian(1); % 龙头 x
        result_xy(2, t_idx) = x_cartesian(2); % 龙头 y

        % 计算每个把手的位置并存储到 result_xy
        x_i = x_t; % 当前把手位置
        for i = 2:224
            % 判断使用 d1 还是 d2
            d = (i == 2) * d1 + (i > 2) * d2;

            % 使用 fzero 求解 beta_i
            beta_i = fzero(@(beta) (x_i(1)^2 + (x_i(1) + a * beta)^2 - d^2) - 2 * x_i(1) * (x_i(1) + a * beta) * cos(beta), 0.1);
            
            % 更新第 i 个把手的极坐标
            x_i = [x_i(1) + a * beta_i, x_i(2) + beta_i];

            % 将当前把手的位置转换为直角坐标
            x_cartesian = [x_i(1) * cos(x_i(2)), x_i(1) * sin(x_i(2))];

            % 记录当前把手直角坐标到 result_xy
            result_xy(2*i-1, t_idx) = x_cartesian(1); % 第 i 把手的 x 坐标
            result_xy(2*i, t_idx) = x_cartesian(2);   % 第 i 把手的 y 坐标
        end
    end
end
