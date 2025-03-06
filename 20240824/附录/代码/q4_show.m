load('q4_xy_before.mat')
load('q4_xy_after.mat')
load('q4_xy_before_delta.mat')
load('q4_xy_after_delta.mat')
delta_time = 0.05;
result_xy = [result_xy_before, result_xy_after];
result_xy_delta = [result_xy_before_delta, result_xy_after_delta];
result_v_xy = (result_xy_delta - result_xy) ./ delta_time;
result_v = zeros(224, 201);

% 设置异常值的阈值 (根据实际情况调整)
v_threshold1 = 1.03;
v_threshold2 = 0.97;

for i = 1:224
    for j = 1:201
        v_x = result_v_xy(2*i-1, j);
        v_y = result_v_xy(2*i, j);
        v = sqrt(v_x^2 + v_y^2); % 计算速度大小
        
        % 判断是否为异常值
        if v > v_threshold1 || isnan(v) || isinf(v) || v < v_threshold2
            v = 1; % 异常值处理，将v设为1
        end
        
        % 存储速度值
        result_v(i,j) = v;
    end
end



% 参数初始化
k = 2; % 调头圆弧的的半径比
alpha = 2 * pi / 3; % 盘出点的圆心角
R = 4.5; % 掉头区域半径
x = 2.805639223409097; % 调头曲线大圆半径
y = 1.402819611704549; % 调头曲线小圆半径
angle1 = 2.450545798771806; % 半径较小的圆弧的圆心角
angle2 = 3.497743349968404; % 半径较大的圆弧的圆心角
L_best = 11.78205037983; % 掉头曲线的长度

l = 1.7; % 螺距
a = l / (2 * pi); % 螺线的参数a
theta_0 = 16 * pi;
theta_in = R / a; % 螺线盘入点的极角
theta_out = theta_in - alpha; % 螺线盘出点的极角
x_in = [R * cos(theta_in), R * sin(theta_in)]; % 螺线盘入点的坐标
x_out = [R * cos(theta_out), R * sin(theta_out)]; % 螺线盘出点的坐标
o1 = x_in * (R - x) / R; % 大圆圆心
o2 = x_out * (R - y) / R; % 小圆圆心

% 初始化笛卡尔坐标图
figure;
hold on;
axis equal;
grid on;

% 生成盘入螺线的角度和半径
theta_vals = linspace(theta_in, theta_0, 1000); % 生成螺线的角度
r_vals = a * theta_vals; % 生成对应的半径
% 将螺线从极坐标转换为笛卡尔坐标
[x_spiral, y_spiral] = pol2cart(theta_vals, r_vals); % 极坐标到笛卡尔坐标的转换
plot(x_spiral, y_spiral, 'k--', 'Color', 'black'); % 绘制参考螺线

theta_vals = linspace(theta_in - alpha, theta_0, 1000); % 生成螺线的角度
% 生成盘出螺线的角度和半径
r_vals = a * (theta_vals + alpha); % 生成对应的半径
% 将螺线从极坐标转换为笛卡尔坐标
[x_spiral, y_spiral] = pol2cart(theta_vals, r_vals); % 极坐标到笛卡尔坐标的转换
plot(x_spiral, y_spiral, 'k--', 'Color', 'b'); % 绘制参考螺线

% 设置笛卡尔坐标下的掉头区域
theta_circle = linspace(0, 2*pi, 100); % 0到2*pi生成圆的角度
r_circle = R * ones(size(theta_circle)); % 圆的半径为R
% 将掉头圆从极坐标转换为笛卡尔坐标
[x_circle, y_circle] = pol2cart(theta_circle, r_circle); % 极坐标到笛卡尔坐标的转换
plot(x_circle, y_circle, 'LineWidth', 2); % 绘制掉头区域的圆

% 生成极坐标下的大圆与小圆
% 大圆的坐标
theta_vals_large = linspace(theta_in - angle1, theta_in, 1000);
theta_vals_small = linspace(theta_out - angle2, theta_out, 1000);
x_large_circle = o1(1) + x * cos(theta_vals_large); % 大圆的x坐标
y_large_circle = o1(2) + x * sin(theta_vals_large); % 大圆的y坐标
% 小圆的坐标
x_small_circle = o2(1) + y * cos(theta_vals_small); % 小圆的x坐标
y_small_circle = o2(2) + y * sin(theta_vals_small); % 小圆的y坐标
plot(x_large_circle, y_large_circle, 'r-', 'LineWidth', 1);
plot(x_small_circle, y_small_circle, 'r-', 'LineWidth', 1);

% 动画循环
for t = 1:201
    % 设置坐标轴范围
    xlim([-10, 10]);  % 调整 X 轴的显示范围
    ylim([-10, 10]);  % 调整 Y 轴的显示范围
    % 清除之前绘制的图像
    cla;
    
    % 重新绘制背景圆与螺线
    plot(x_spiral, y_spiral, 'k--', 'Color', 'black');
    plot(x_circle, y_circle, 'LineWidth', 2);
    plot(x_large_circle, y_large_circle, 'r-', 'LineWidth', 1);
    plot(x_small_circle, y_small_circle, 'r-', 'LineWidth', 1);

    % 获取每个时间步 t 下的点的坐标
    x_coords = result_xy(1:2:end, t); % 获取所有点的横坐标
    y_coords = result_xy(2:2:end, t); % 获取所有点的纵坐标
    
    % 跳过 NaN 和 (0, 0) 坐标点
    valid_indices = ~isnan(x_coords) & ~isnan(y_coords) & ~(x_coords == 0 & y_coords == 0);
    x_coords = x_coords(valid_indices);
    y_coords = y_coords(valid_indices);
    
    % 绘制第一个点为红色点
    plot(x_coords(1), y_coords(1), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
    
    % 绘制其余点为蓝色点
    plot(x_coords(2:end), y_coords(2:end), 'bo', 'MarkerSize', 5, 'MarkerFaceColor', 'b');
    
    % 用绿色线连接相邻的点
    for i = 2:length(x_coords)
        line([x_coords(i-1) x_coords(i)], [y_coords(i-1) y_coords(i)], 'Color', 'g', 'LineWidth', 1.5);
    end
    
    % 显示当前时间
    time_text = ['Time: ', num2str(t-101), ' s'];
    text(-9, 9, time_text, 'FontSize', 12, 'Color', 'k', 'FontWeight', 'bold'); % 显示在左上角
    
    % 暂停一段时间以创建动画效果
    pause(0.1); % 可调整暂停时间以控制动画速度
end



