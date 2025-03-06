function intersection_points = find_circle_intersections(center1, r1, center2, r2)
    % center1 = [x1, y1]: 第一个圆的圆心坐标
    % r1: 第一个圆的半径
    % center2 = [x2, y2]: 第二个圆的圆心坐标
    % r2: 第二个圆的半径

    % 解的输出初始化
    intersection_points = [-100, -100; -100, -100];

    % 获取圆心的坐标
    x1 = center1(1);
    y1 = center1(2);
    x2 = center2(1);
    y2 = center2(2);

    % 计算圆心之间的距离
    d = sqrt((x2 - x1)^2 + (y2 - y1)^2);

    % 检查是否有解的条件
    if d > (r1 + r2) || d < abs(r1 - r2)
        % 圆相离或一个圆在另一个圆内，没有交点
        return;
    elseif d == 0 && r1 == r2
        % 两个圆完全重合，无穷多个交点
        return;
    end

    % 计算交点的公式
    a = (r1^2 - r2^2 + d^2) / (2 * d);
    h = sqrt(r1^2 - a^2);

    % 圆心到交点的直线上的中间点
    x0 = x1 + a * (x2 - x1) / d;
    y0 = y1 + a * (y2 - y1) / d;

    % 交点1
    x_intersect1 = x0 + h * (y2 - y1) / d;
    y_intersect1 = y0 - h * (x2 - x1) / d;

    % 交点2
    x_intersect2 = x0 - h * (y2 - y1) / d;
    y_intersect2 = y0 + h * (x2 - x1) / d;

    % 输出解
    if d == r1 + r2 || d == abs(r1 - r2)
        % 圆相切，只有一个交点
        intersection_points = [x_intersect1, y_intersect1; -100, -100];
    else
        % 圆相交，有两个交点
        intersection_points = [x_intersect1, y_intersect1; x_intersect2, y_intersect2];
    end
end