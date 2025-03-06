function r = calc_collision_r(L)
    % 输入螺距L，返回停止时龙头到原点的距离
    d1 = 2.86; % 第一个和第二个把手的距离
    d2 = 1.65; % 第i个把手与第i+1个把手的距离
    a = L / (2 * pi); % 螺线参数
    v = 1; % 速度
    theta_0 = 32 * pi; % 初始角度
    % 计算螺线长度函数
    l_theta = @(theta) a/2 * (theta .* sqrt(1 + theta.^2) + log(theta + sqrt(1 + theta.^2))); 
    % 初始螺线长度
    l_theta_0 = l_theta(theta_0);

    % 搜索起点为0，时间步长为10s
    t_start = 0;
    delta = 10;
    % 搜索终点为螺线长度/龙头速度
    t_end = l_theta_0;
    result_xy_10 = generating_position(t_start, t_end, delta, L);
    t_lst = t_start:delta:t_end;
    n = length(t_lst);
    t_start = t_lst(1);
    t_end = t_lst(n);
    flag_lst = zeros(1, n);

    % 循环遍历每个时间步
    for k = 1:n
        flag = 0;  % 初始化当前时间步的标志位
    % 遍历前32个把手
    for i = 1:32
        lst = result_xy_10(:, k)';
        
        % 获取第i个把手的坐标 (x1, y1), (x2, y2)
        x1 = lst(2*i-1);
        y1 = lst(2*i);
        if i < 224
            x2 = lst(2*i+1);
            y2 = lst(2*i+2);
        end
        
        % 遍历后续把手，判断是否可能发生碰撞
        for j = i+2:223
            % 获取第j个把手的坐标 (x3, y3), (x4, y4)
            x3 = lst(2*j-1);
            y3 = lst(2*j);
            if j < 224
                x4 = lst(2*j+1);
                y4 = lst(2*j+2);
            end
            
            % 计算四点之间的最小距离
            d1 = sqrt((x1 - x3)^2 + (y1 - y3)^2);
            d2 = sqrt((x1 - x4)^2 + (y1 - y4)^2);
            d3 = sqrt((x2 - x3)^2 + (y2 - y3)^2);
            d4 = sqrt((x2 - x4)^2 + (y2 - y4)^2);
            d_min = min([d1, d2, d3, d4]);
            
            % 如果最小距离小于1，检查是否碰撞
            if d_min < 1
                flag = if_coordinates(x1, y1, x2, y2, x3, y3, x4, y4);
            end

            % 如果发生碰撞，跳出循环
            if flag
                break
            end
        end

        % 如果发生碰撞，跳出外层循环
        if flag
            break
        end
    end
    
        % 存储当前时间步的碰撞标志位
        flag_lst(k) = flag;
    end

    for i = 1:n
        if flag_lst(i)
            t_end = t_lst(i);
            break
        end
    end
    t_start = t_end - 10;

    % 新的搜索范围，步长为0.1
    delta_fine = 0.1;  % 细致搜索的步长为0.1
    result_xy_fine = generating_position(t_start, t_end, delta_fine, L);
    n_fine = (t_end - t_start) / delta_fine + 1;
    flag_lst_fine = zeros(1, n_fine);

    % 循环遍历每个时间步，进行精细搜索
    for k = 1:n_fine
        flag = 0;  % 初始化当前时间步的标志位
    
        % 遍历前32个把手
        for i = 1:32
            lst = result_xy_fine(:, k)';
        
            % 获取第i个把手的坐标 (x1, y1), (x2, y2)
            x1 = lst(2*i-1);
            y1 = lst(2*i);
            if i < 224
                x2 = lst(2*i+1);
                y2 = lst(2*i+2);
            end
        
            % 遍历后续把手，判断是否可能发生碰撞
            for j = i+2:223
                % 获取第j个把手的坐标 (x3, y3), (x4, y4)
                x3 = lst(2*j-1);
                y3 = lst(2*j);
                if j < 224
                    x4 = lst(2*j+1);
                    y4 = lst(2*j+2);
                end
            
                % 计算四点之间的最小距离
                d1 = sqrt((x1 - x3)^2 + (y1 - y3)^2);
                d2 = sqrt((x1 - x4)^2 + (y1 - y4)^2);
                d3 = sqrt((x2 - x3)^2 + (y2 - y3)^2);
                d4 = sqrt((x2 - x4)^2 + (y2 - y4)^2);
                d_min = min([d1, d2, d3, d4]);
            
                % 如果最小距离小于1，检查是否碰撞
                if d_min < 1
                    flag = if_coordinates(x1, y1, x2, y2, x3, y3, x4, y4);
                end

                % 如果发生碰撞，跳出循环
                if flag
                    break
                end
            end

            % 如果发生碰撞，跳出外层循环
            if flag
                break
            end
        end
    
        % 存储当前时间步的碰撞标志位
        flag_lst_fine(k) = flag;
    end

    % 更新最终时间
    time_lst_fine = linspace(t_start, t_end, n_fine);
    for i = 1:n_fine
        if flag_lst_fine(i)
            t_end_fine = time_lst_fine(i);
            break
        end
    end



    result_xy_final = generating_position(t_end_fine, t_end_fine, 0, L);
    x = result_xy_final(1, 1);
    y = result_xy_final(2, 1);
    r = sqrt(x^2 + y^2);
end