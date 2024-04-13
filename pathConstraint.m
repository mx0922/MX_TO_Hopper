function [c, ceq] = pathConstraint(t, x, u, p)

% 这里主要就是运动学约束:即脚的位置不能超过那个cube
% 其实也可以将接触力与脚的位置加入到路径约束中

Nt = length(t);

rBase = x(1:3, :);
qBase = x(4:6, :);

P_foot = u(1:3, :);
F_foot = u(4:6, :);

%% 1.脚的运动学约束
p_hat = p.p_hat;
L_cube = p.L_cube;

deltaP = zeros(Nt, 1);
for ii = 1:Nt
    % 利用欧拉角得到Base到World的旋转矩阵
    R = getRotMatFromBaseToWorld(qBase(:, ii));
    
    % 脚相对于base的位置在base坐标系下的表示
    p_Base = R \ (P_foot(:, ii) - rBase(:, ii));
    
    % 得到与nominal位置的距离差
    delta_p = p_Base - p_hat;
    
    % 多处理一步：得到其模的平方
    deltaP(ii) = delta_p' * delta_p;
end

% deltaP <= L_cube^2
c1 = deltaP - L_cube^2;

%% 接触力与脚的位置
miu = p.miu;

fx = F_foot(1, :);
fy = F_foot(2, :);
fz = F_foot(3, :);

% fz >= 0
c2 = reshape(-fz, Nt, 1);

% 摩擦锥
c3 = [-fx - miu * fz, fx - miu * fz, -fy - miu * fz, fy - miu * fz];
c3 = reshape(c3, numel(c3), 1);

% p_foot_z >= 0
c4 = reshape(-P_foot(3, :), Nt, 1);

%% 整合
c = [c1; c2; c3; c4];
ceq = [];

end