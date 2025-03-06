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

% 计算螺线长度函数
l_theta = @(theta) a/2 * (theta .* sqrt(1 + theta.^2) + log(theta + sqrt(1 + theta.^2))); 
l_theta_0 = l_theta(theta_in);
t_lst = 1:100;
% 初始化龙头把手位置
result_xy_after_delta = zeros(448,100);


delta_time = 0.05;
% 初始化正时间内龙头把手位置
l_lst_after = 0.95:1:99.955;
for t = 1:100
    t_ = l_lst_after(t);
    if t_ < L_best
        position = route2position(t_);
        result_xy_after_delta(1, t) = position(1);
        result_xy_after_delta(2, t) = position(2);
    else
        l_temp = t_ - L_best + l_theta_0;
        theta_t = fzero(@(theta) l_theta(theta) - l_temp, theta_0);
        theta_temp = theta_t - alpha;
        r_temp = theta_t * a;
        x_temp = r_temp * cos(theta_temp);
        y_temp = r_temp * sin(theta_temp);
        result_xy_after_delta(1, t) = x_temp;
        result_xy_after_delta(2, t) = y_temp;
    end
end




tolerance = 0.1;
% 储存正时刻位置
for t = 1:100
    current_position = result_xy_after_delta(1:2, t)';
    [theta_temp, r_temp] = cart2pol(current_position(1), current_position(2));
    flag = 1;
    for i = 2:224
        d = (i == 2) * d1 + (i > 2) * d2; % 判断使用d1还是d2

        % 计算是否在盘出螺线上，若在盘出螺线上，计算极角
        theta_real = r_temp / a - alpha;
        delta_theta = theta_real - theta_temp - fix((theta_real - theta_temp) / (2 * pi)) * 2 * pi;
        
   
        % 若出现一次不在盘出螺线上，则后续点不会在盘出螺线上，不用考虑盘出螺线的检验
        if abs(delta_theta) > tolerance
            flag = 0;
        end

        if r_temp <= R
            % 计算当前位置对应的路径长度
            l_answer = position2route(current_position);
            points12 = find_circle_intersections(current_position, d, o1, x_radius);
            points34 = find_circle_intersections(current_position, d, o2, y_radius);
            points1 = points12(1, :);
            points2 = points12(2, :);
            points3 = points34(1, :);
            points4 = points34(2, :);
            l1 = position2route(points1);
            l2 = position2route(points2);
            l3 = position2route(points3);
            l4 = position2route(points4);
            l_list = [l1, l2, l3, l4];
            for j = 1:4
                if l_list(j) == -1 || l_list(j) >= l_answer
                    l_list(j) = 0;
                end
            end
            if max(l_list) > 0
                x_cartesian = route2position(max(l_list));
                % 存储当前把手的直角坐标到 result
                result_xy_after_delta(2*i-1, t) = x_cartesian(1); % 第i节龙身x
                result_xy_after_delta(2*i, t) = x_cartesian(2);   % 第i节龙身y
            else
                % 查找失败
                x_temp = current_position(1);
                y_temp = current_position(2);
                theta_temp = fzero(@(x) (a*x*cos(x)-x_temp)^2+(a*x*sin(x)-y_temp)^2-d^2,theta_in);
                r_temp = a * theta_temp;
                x_cartesian = [r_temp * cos(theta_temp), r_temp * sin(theta_temp)];
                result_xy_after_delta(2*i-1, t) = x_cartesian(1); % 第i节龙身x
                result_xy_after_delta(2*i, t) = x_cartesian(2);   % 第i节龙身y
            end
        elseif r_temp > R && flag
            % 计算 beta_i
            beta = fzero(@(beta) (r_temp^2 + (r_temp-a*beta)^2 - d^2) - 2 * r_temp * (r_temp+a*beta) * cos(beta), [0, pi/2]);
            theta_ = theta_real - beta;
            r_ = (theta_ + alpha) * a;
            if r_ >= R
                theta_temp = theta_;
                r_temp = r_;
                x_cartesian = [r_temp * cos(theta_temp), r_temp * sin(theta_temp)];
                result_xy_after_delta(2*i-1, t) = x_cartesian(1); % 第i节龙身x
                result_xy_after_delta(2*i, t) = x_cartesian(2);   % 第i节龙身y
            else
                % 计算当前位置对应的路径长度
                l_answer = L_best;
                points12 = find_circle_intersections(current_position, d, o1, x_radius);
                points34 = find_circle_intersections(current_position, d, o2, y_radius);
                points1 = points12(1, :);
                points2 = points12(2, :);
                points3 = points34(1, :);
                points4 = points34(2, :);
                l1 = position2route(points1);
                l2 = position2route(points2);
                l3 = position2route(points3);
                l4 = position2route(points4);
                l_list = [l1, l2, l3, l4];
                for j = 1:4
                    if l_list(j) == -1 
                        l_list(j) = 0;
                    end
                end
                if max(l_list) > 0
                    x_cartesian = route2position(max(l_list));
                    % 存储当前把手的直角坐标到 result
                    result_xy_after_delta(2*i-1, t) = x_cartesian(1); % 第i节龙身x
                    result_xy_after_delta(2*i, t) = x_cartesian(2);   % 第i节龙身y
                end
            end
        elseif r_temp > R
            % 计算 beta_i
            beta = fzero(@(beta) (r_temp^2 + (r_temp-a*beta)^2 - d^2) - 2 * r_temp * (r_temp+a*beta) * cos(beta), [0, pi/2]);
            theta_temp = theta_temp + beta;
            [x_temp, y_temp] = pol2cart(theta_temp, a * theta_temp);
            % 存储当前把手的直角坐标到 result
            result_xy_after_delta(2*i-1, t) = x_temp; % 第i节龙身x
            result_xy_after_delta(2*i, t) = y_temp;   % 第i节龙身y
        end
        current_position(1) = result_xy_after_delta(2*i-1, t);
        current_position(2) = result_xy_after_delta(2*i, t);
        r_temp = sqrt(dot(current_position, current_position));
    end
end

