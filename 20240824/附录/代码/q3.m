R = 4.5;
L_upper = 0.70;
L_lower = 0.28;
L_lst = linspace(L_lower, L_upper, 47);
r_lst = [];
for L = L_lst
    r = calc_collision_r(L);
    r_lst = [r_lst, r];
end

L_upper = 0.46;
L_lower = 0.43;
L_lst_fine = L_lower:0.0001:L_upper;
r_lst_fine = [];
for L = L_lst_fine
    r = calc_collision_r(L);
    r_lst_fine = [r_lst_fine, r];
end


% 创建一个新的图形窗口
figure;
% 绘制第一个子图
plot(L_lst, r_lst, 'o-', 'MarkerSize', 2);  % 绘制折线
hold on
plot([0.25, 0.75], [4.5, 4.5], 'Color', 'r', 'LineWidth', 1, 'LineStyle', '--')
xlim([0.28, 0.70])
xlabel('螺距 (m)');        % 设置X轴标签
ylabel('最小碰撞半径');        % 设置Y轴标签
grid on;                   % 显示网格

figure
% 绘制第一个子图
plot(L_lst_fine, r_lst_fine, 'o-', 'MarkerSize', 2);
hold on
plot([0.25, 0.75], [4.5, 4.5], 'Color', 'r', 'LineWidth', 1, 'LineStyle', '--')
xlim([0.43, 0.46])
xlabel('螺距 (m)');        % 设置X轴标签
ylabel('最小碰撞半径');        % 设置Y轴标签
grid on;                  % 显示网格