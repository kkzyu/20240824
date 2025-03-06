function position = route2position(l)
    % 参数初始化
    R = 4.5; % 掉头区域半径
    L = 1.7; % 螺距
    a = L / (2 * pi); % 螺线的参数a
    alpha = 2 * pi / 3; % 盘出点的圆心角
    x_radius = 2.805639223409097; % 调头曲线大圆半径
    y_radius = 1.402819611704549; % 调头曲线小圆半径
    angle1 = 2.450545798771806; % 大圆弧的圆心角
    angle2 = 3.497743349968404; % 小圆弧的圆心角
    L_best = 11.782050379839399; % 最优曲线长度
    
    % 螺线盘入点和盘出点的极角
    theta_in = R / a; % 螺线盘入点的极角
    theta_out = theta_in - alpha; % 螺线盘出点的极角
    x_in = [R * cos(theta_in), R * sin(theta_in)]; % 螺线盘入点的坐标
    x_out = [R * cos(theta_out), R * sin(theta_out)]; % 螺线盘出点的坐标
    
    % 圆心坐标
    o1 = x_in * (R - x_radius) / R; % 大圆的圆心
    o2 = x_out * (R - y_radius) / R; % 小圆的圆心

    % 计算大圆和小圆的弧长
    l1 = x_radius * angle1;
    l2 = y_radius * angle2;
    
    % 判断路径在大圆还是小圆
    if l < l1
        % 在大圆上
        theta_l = theta_in - l / x_radius; % 计算对应的角度
        x_temp = o1(1) + x_radius * cos(theta_l);
        y_temp = o1(2) + x_radius * sin(theta_l);
        position = [x_temp, y_temp];
    else
        % 在小圆上
        remaining_length = L_best - l;
        theta_l = theta_out - remaining_length / y_radius; % 计算对应的角度

        x_temp = o2(1) + y_radius * cos(theta_l);
        y_temp = o2(2) + y_radius * sin(theta_l);
        position = [x_temp, y_temp];
    end
end
