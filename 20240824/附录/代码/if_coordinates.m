function flag = if_coordinates(x1, y1, x2, y2, x3, y3, x4, y4)
    % 输入四个点的坐标，判定两块板是否会相撞
    % x1, y1, x2, y2 表示内侧第一块板的两端坐标
    % x3, y3, x4, y4 表示外侧第二块板的两端坐标

    % 计算向量
    vector1 = [x1 - x2, y1 - y2]; % 第一块板的向量
    d = sqrt(dot(vector1, vector1)); % 向量长度
    D1 = [x1, y1] + vector1 * 0.275 / d; % D1 为第一块板的一端
    D2 = [x2, y2] - vector1 * 0.275 / d; % D2 为第一块板的另一端
    vector2 = [y2 - y1, x1 - x2]; % 法向量
    vector2 = vector2 * 0.15 / d; % 调整法向量的大小

    % 初始化两个碰撞标志
    flag1 = 0;
    flag2 = 0;

    % 计算D1和D2距离
    distance_D1 = dot(D1, D1);
    distance_D2 = dot(D2, D2);

    % 计算C1
    if dot(D1 + vector2, D1 + vector2) > distance_D1
        C1 = D1 + vector2;
    else
        C1 = D1 - vector2;
    end

    % 计算C2
    if dot(D2 + vector2, D2 + vector2) > distance_D2
        C2 = D2 + vector2;
    else
        C2 = D2 - vector2;
    end

    % 第二块板的坐标
    B1 = [x3, y3];
    B2 = [x4, y4];
    d_B1B2 = 2.2 - 0.275 * 2; % 第二块板的端点距离

    % 计算C1能否碰撞
    vector1_cross = (B1(1) - C1(1)) * (B2(2) - C1(2)) - (B1(2) - C1(2)) * (B2(1) - C1(1)); % 计算外积
    distance1 = abs(vector1_cross) / d_B1B2; % 距离计算
    if distance1 < 0.15
        B1C1 = C1 - B1;
        B1B2 = B2 - B1;
        lamda = dot(B1C1, B1B2) / dot(B1B2, B1B2);
        B1H = lamda * B1B2;
        B2H = -B1B2 + B1H;
        dot1 = dot(B1H, B2H);
        if dot1 < 0.275 * (0.275 + d_B1B2)
            flag1 = 1;
        end
    end

    % 计算C2能否碰撞
    vector2_cross = (B1(1) - C2(1)) * (B2(2) - C2(2)) - (B1(2) - C2(2)) * (B2(1) - C2(1)); % 计算外积
    distance2 = abs(vector2_cross) / d_B1B2; % 距离计算
    if distance2 < 0.15
        B1C2 = C2 - B1;
        B1B2 = B2 - B1;
        lamda = dot(B1C2, B1B2) / dot(B1B2, B1B2);
        B1H = lamda * B1B2;
        B2H = -B1B2 + B1H;
        dot1 = dot(B1H, B2H);
        if dot1 < 0.275 * (0.275 + d_B1B2)
            flag2 = 1;
        end
    end

    flag = flag1 || flag2;
end
