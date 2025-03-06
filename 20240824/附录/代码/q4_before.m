l = 1.7; % 螺距
d1 = 2.86; % 第一个和第二个把手的距离
d2 = 1.65; % 第i个把手与第i+1个把手的距离
a = l / (2 * pi); % 螺线参数
v = -1; % 速度
theta_0 = 16.6320; % 初始角度

% 计算螺线长度函数
l_theta = @(theta) a/2 * (theta .* sqrt(1 + theta.^2) + log(theta + sqrt(1 + theta.^2))); 

% 初始螺线长度
l_theta_0 = l_theta(theta_0);

% 时间步长
t_max = 100;
dt = 1;

% 初始化存储数组，448行 (x,y 交替), 100列 (t_max)
result_xy = zeros(448, t_max+1);

% 初始化极坐标图
polar_axes = polaraxes;
hold on;
theta_vals = linspace(0, theta_0, 1000);
r_vals = a * theta_vals;
polarplot(polar_axes, theta_vals, r_vals, 'k--'); % 绘制参考螺线

% 初始化龙头轨迹
x_1_history = []; % 存储第一个把手的轨迹
x_1 = [a * theta_0, theta_0]; % 龙头的极坐标
x_cartesian = [x_1(1) * cos(x_1(2)), x_1(1) * sin(x_1(2))]; % 直角坐标

% 逐步计算时间t时刻的螺线长度、各个把手的位置及其极坐标
for t = 0:dt:t_max
    l_t = l_theta_0 - v * t; % t时刻螺线长度
    theta_t = fzero(@(theta) l_theta(theta) - l_t, theta_0); % 利用螺线长度求解 theta_t
    x_t = [a * theta_t, theta_t]; % 当前时间 t 时的极坐标

    % 记录龙头轨迹
    x_1_history = [x_1_history; x_t(2), x_t(1)];

    % 清除上一帧的把手和连接线
    cla(polar_axes);

    % 绘制参考螺线
    polarplot(polar_axes, theta_vals, r_vals, 'k--');

    % 更新并绘制龙头位置（红色点）和轨迹
    if size(x_1_history, 1) > 1
        polarplot(polar_axes, x_1_history(:, 1), x_1_history(:, 2), 'w', 'LineWidth', 1.5); % 红色轨迹
    end
    polarplot(polar_axes, x_t(2), x_t(1), 'ro', 'MarkerFaceColor', 'r'); % 红色点标记龙头

    % 计算每个把手的位置并连接
    x_i = x_t; % 当前把手位置
    % 逆向存储，第 t_max+1-t 列存储当前时间步的数据
    result_xy(1, t_max+1-t) = x_t(1) * cos(x_t(2)); % x 坐标
    result_xy(2, t_max+1-t) = x_t(1) * sin(x_t(2)); % y 坐标

    for i = 2:224
        d = (i == 2) * d1 + (i > 2) * d2; % 判断使用 d1 还是 d2

        % 计算 beta_i
        beta_i = fzero(@(beta) (x_i(1)^2 + (x_i(1) + a * beta)^2 - d^2) - 2 * x_i(1) * (x_i(1) + a * beta) * cos(beta), 0.1);

        % 计算第 i 个把手的极坐标
        prev_x_i = x_i; % 保存上一个把手的位置
        x_i = [x_i(1) + a * beta_i, x_i(2) + beta_i];
        x_i_cartesian = [x_i(1) * cos(x_i(2)), x_i(1) * sin(x_i(2))];
        prev_x_i_cartesian = [prev_x_i(1) * cos(prev_x_i(2)), prev_x_i(1) * sin(prev_x_i(2))];

        % 绘制把手连接线
        polarplot(polar_axes, [prev_x_i(2), x_i(2)], [prev_x_i(1), x_i(1)], 'g-', 'LineWidth', 1);
        % 用蓝色点表示把手的位置
        polarplot(polar_axes, x_i(2), x_i(1), 'bo', 'MarkerFaceColor', 'b');
        
        % 逆向存储，第 t_max+1-t 列存储当前时间步的数据
        result_xy(2*i-1, t_max+1-t) = x_i_cartesian(1); % x 坐标 (奇数行)
        result_xy(2*i, t_max+1-t) = x_i_cartesian(2);   % y 坐标 (偶数行)
    end

    % 显示时间
    text(1/4*pi, 2* max(r_vals), ['t = ' num2str(t),'s'], 'HorizontalAlignment', 'center', 'FontSize', 12);

    % 暂停以展示动画效果
    pause(0.001);
end

hold off;
result_xy_before = result_xy;
% 存储结果到文件
save('points_data_448x100_reversed.mat', 'result_xy_before');
