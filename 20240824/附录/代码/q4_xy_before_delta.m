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
delta_t = 0.05; % 用于计算瞬时速度的时间增量

% 初始化存储数组，448行 (x,y 交替), 100列 (t_max)
result_xy = zeros(448, t_max+1);
result_xy_deltax = zeros(448, t_max+1); % 用于存储delta_t后的位置
result_v = zeros(224, t_max+1); % 用于存储每个把手的瞬时速度

% 逐步计算时间t时刻的螺线长度、各个把手的位置及其极坐标
for t = 0:dt:t_max
    l_t = l_theta_0 - v * t; % t时刻螺线长度
    theta_t = fzero(@(theta) l_theta(theta) - l_t, theta_0); % 利用螺线长度求解 theta_t
    x_t = [a * theta_t, theta_t]; % 当前时间 t 时的极坐标
    % 计算当前把手的位置并存储
    x_i = x_t; % 当前把手位置
    result_xy(1, t_max+1-t) = x_t(1) * cos(x_t(2)); % x 坐标
    result_xy(2, t_max+1-t) = x_t(1) * sin(x_t(2)); % y 坐标

    % 计算 delta_t 后的螺线长度
    l_t_deltat = l_theta_0 - v * (t+delta_t); % delta_t 时间后的螺线长度
    % 计算 delta_t 后的位置
    theta_t_deltat = fzero(@(theta) l_theta(theta) - l_t_deltat, delta_t + theta_0); % 计算新的 theta
    x_t_deltat = [a * theta_t_deltat, theta_t_deltat];
    x_i_deltat = x_t_deltat;
    result_xy_deltax(1, t_max+1-t) = x_t_deltat(1) * cos(x_t_deltat(2)); % x 坐标
    result_xy_deltax(2, t_max+1-t) = x_t_deltat(1) * sin(x_t_deltat(2)); % y 坐标
    % 计算第一个把手的瞬时速度
    delta_x_1 = sqrt((x_t_deltat(1) * cos(x_t_deltat(2)) - x_t(1) * cos(x_t(2)))^2 + ...
                     (x_t_deltat(1) * sin(x_t_deltat(2)) - x_t(1) * sin(x_t(2)))^2);
    result_v(1, t_max+1-t) = delta_x_1 / delta_t; % 存储第一个把手的瞬时速度

    for i = 2:224
        d = (i == 2) * d1 + (i > 2) * d2; % 判断使用 d1 还是 d2

        % 计算 beta_i
        beta_i = fzero(@(beta) (x_i(1)^2 + (x_i(1) + a * beta)^2 - d^2) - 2 * x_i(1) * (x_i(1) + a * beta) * cos(beta), 0.1);

        % 计算第 i 个把手的极坐标
        x_i = [x_i(1) + a * beta_i, x_i(2) + beta_i];
        x_i_cartesian = [x_i(1) * cos(x_i(2)), x_i(1) * sin(x_i(2))];

        % 逆向存储，第 t_max+1-t 列存储当前时间步的数据
        result_xy(2*i-1, t_max+1-t) = x_i_cartesian(1); % x 坐标 (奇数行)
        result_xy(2*i, t_max+1-t) = x_i_cartesian(2);   % y 坐标 (偶数行)

        % 计算 beta_i_deltat
        beta_i_deltat = fzero(@(beta) (x_i_deltat(1)^2 + (x_i_deltat(1) + a * beta)^2 - d^2) - 2 * x_i_deltat(1) * (x_i_deltat(1) + a * beta) * cos(beta), 0.1);

        % 计算 delta_t 后的位置
        x_i_deltat = [x_i_deltat(1)+a*beta_i_deltat,x_i_deltat(2)+beta_i_deltat];
        x_i_deltat_cartesian = [x_i_deltat(1) * cos(x_i_deltat(2)), x_i_deltat(1) * sin(x_i_deltat(2))];

        % 存储 delta_t 后的位置
        result_xy_deltax(2*i-1, t_max+1-t) = x_i_deltat_cartesian(1); % x 坐标 (奇数行)
        result_xy_deltax(2*i, t_max+1-t) = x_i_deltat_cartesian(2);   % y 坐标 (偶数行)

        % 计算瞬时速度
        delta_x = sqrt((x_i_deltat_cartesian(1) - x_i_cartesian(1))^2 + (x_i_deltat_cartesian(2) - x_i_cartesian(2))^2);
        result_v(i, t_max+1-t) = delta_x / delta_t; % 存储瞬时速度
    end
end
result_xy_before_delta = result_xy_deltax;
% 存储结果到文件
save('q4_xy_before_delta.mat', 'result_xy_before_delta');

