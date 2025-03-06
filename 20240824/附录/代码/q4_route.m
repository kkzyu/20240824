% 初始化参数
l = 1.7; % 螺距为1.7m
a_temp = l / (2 * pi); % 螺线参数为a
R = 4.5; % 直径为9m
theta = R / a_temp;
tan_angle = 1 / theta;
angle = atan(tan_angle); % 盘入螺线在盘入处的切线与圆的切线存在一个较小的夹角

% 定义参考半径比和旋转角列表
k_lst = [1, 1.5, 2, 2.5, 3, 4]; % 参考半径比
alpha_lst = linspace(pi/2, pi, 6); % 参考旋转角

% 预分配结果矩阵
L_matrix = zeros(length(k_lst), length(alpha_lst));

% 双重循环遍历所有可能的k和alpha组合
for i = 1:length(alpha_lst)
    alpha = alpha_lst(i);
    for j = 1:length(k_lst)
        k = k_lst(j);
        a_temp = (1 + cos(alpha)) * k;
        b_temp = (1 - cos(alpha)) * R * (k + 1);
        c_temp = (cos(alpha) - 1) * R^2;
        delta = b_temp^2 - 4 * a_temp * c_temp;
        x = (-b_temp + sqrt(delta)) / (2 * a_temp);
        y = k * x;
        side1 = R - y;
        side2 = R - x;
        side3 = x + y;
        angle1 = pi - acos((side1^2 + side3^2 - side2^2) / (2 * side1 * side3));
        angle2 = pi + angle1 - alpha;
        L_matrix(j, i) = angle1 * y + angle2 * x; % 将结果存入矩阵
    end
end

L_matrix(:, 6) = pi * R;

% 转换 alpha 从弧度到度
alpha_lst_deg = alpha_lst * 180 / pi;

% 创建折线图
figure;
hold on;
for i = 1:length(k_lst)
    plot(alpha_lst_deg, L_matrix(i, :), '-o', 'DisplayName', sprintf('k = %.1f', k_lst(i)));
end

% 设置图例和标签
legend('Location', 'northwest');
xticks(alpha_lst_deg);
xticklabels(arrayfun(@(alpha) sprintf('%.0f', alpha), alpha_lst_deg, 'UniformOutput', false));
xlabel('盘出点圆心角 ω (°)');
ylabel('掉头曲线长度 (m)');

% 显示网格
grid on;
hold off;

% 最短掉头曲线的计算
% 定义参考半径比和旋转角
L_min = 0.45; % q3中的最小螺距
k = 2; % 用户输入的半径比
alpha = 2 * pi * L_min / l;

% 掉头曲线的各项参数
a_temp = (1 + cos(alpha)) * k;
b_temp = (1 - cos(alpha)) * R * (k + 1);
c_temp = (cos(alpha) - 1) * R^2;
delta = b_temp^2 - 4 * a_temp * c_temp;
x = (-b_temp + sqrt(delta)) / (2 * a_temp); % 较小圆弧的半径
y = k * x; % 较大圆弧的半径
side1 = R - y;
side2 = R - x;
side3 = x + y;
angle1 = pi - acos((side1^2 + side3^2 - side2^2) / (2 * side1 * side3));
angle2 = pi + angle1 - alpha

L_shortest = angle1 * y + angle2 * x;

% 参数初始化
k = 2; % 调头圆弧的的半径比
alpha = 2 * pi / 3; % 盘出点的圆心角
R = 4.5; % 掉头区域半径
x_radius = 2.805639223409097; % 调头曲线大圆半径
y_radius = 1.402819611704549; % 调头曲线小圆半径
angle1 = 2.450545798771806; % 半径较大的圆弧的圆心角
angle2 = 3.497743349968404; % 半径较小的圆弧的圆心角
L_best = 11.782050379839399; % 掉头曲线的长度

% 盘入螺线方程:r=a*theta
% 盘出螺线方程:r=a*(theta+alpha)

d1 = 2.86; % 第一个和第二个把手的距离
d2 = 1.65; % 第i个把手与第i+1个把手的距离
l = 1.7; % 螺距
a = l / (2 * pi); % 螺线的参数a
theta_0 = 16 * pi;
theta_in = R / a; % 螺线盘入点的极角
theta_out = theta_in - alpha; % 螺线盘出点的极角
x_in = [R * cos(theta_in), R * sin(theta_in)]; % 螺线盘入点的坐标
x_out = [R * cos(theta_out), R * sin(theta_out)]; % 螺线盘出点的坐标
o1 = x_in * (R - x_radius) / R; % 大圆圆心
o2 = x_out * (R - y_radius) / R; % 小圆圆心


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

theta_vals = linspace(theta - alpha, theta_0, 1000); % 生成螺线的角度
% 生成盘出螺线的角度和半径
r_vals = a * (theta_vals + alpha); % 生成对应的半径
% 将螺线从极坐标转换为笛卡尔坐标
[x_spiral, y_spiral] = pol2cart(theta_vals, r_vals); % 极坐标到笛卡尔坐标的转换
plot(x_spiral, y_spiral, 'k-', 'Color', 'b'); % 绘制参考螺线

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
x_large_circle = o1(1) + x_radius * cos(theta_vals_large); % 大圆的x坐标
y_large_circle = o1(2) + x_radius * sin(theta_vals_large); % 大圆的y坐标
% 小圆的坐标
x_small_circle = o2(1) + y_radius * cos(theta_vals_small); % 小圆的x坐标
y_small_circle = o2(2) + y_radius * sin(theta_vals_small); % 小圆的y坐标
plot(x_large_circle, y_large_circle, 'r-', 'LineWidth', 1);
plot(x_small_circle, y_small_circle, 'r-', 'LineWidth', 1);
xlim([-8, 8]);
ylim([-8, 8]);
% 在图上添加x_in点的标注
text(x_in(1)+2, x_in(2), '盘入点', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
% 在图上添加x_out点的标注
text(x_out(1)-0.5, x_out(2)-0.3, '盘出点', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');

