function dx = dynamics(t, x, u, p)

% 我这里定义 x = [q; dq] 也即是Base的位置和姿态及其导数
% u 这里认为是[P_foot; F_foot]也即是脚的位置与接触力，默认我已经利用相关函数得到t时刻对应的脚的位置与接触力了

Nt = length(t);

ddrBase = zeros(3, Nt);
ddqBase = zeros(3, Nt);

rBase = x(1:3, :);
qBase = x(4:6, :);

drBase = x(7:9, :);
dqBase = x(10:12, :);

% 脚的位置与接触力
P_foot = u(1:3, :);
F_foot = u(4:6, :);

m = p.m;
g = p.g;
G = [0; 0; g];
I = p.I;

for ii = 1:Nt
    % m * ddr = f - m*G
    ddrBase(:, ii) = F_foot(:, ii) / m - G;
    
    % I*dw + cross(w, I*w) = cross(f, r - p)
    q_ii = qBase(:, ii);
    dq_ii = dqBase(:, ii);
    
    % Euler converter
    % w = Cq * dq; dw = dCq * dq + Cq * ddq
    Cq = getEulerConverter(q_ii);
    dCq = getEulerConverterDiv(q_ii, dq_ii);
    
    % w = Cq * dq
    w_ii = Cq * dq_ii;
    
    rhs_mx = cross(F_foot(:, ii), rBase(:, ii) - P_foot(:, ii));
    dw_ii = I \ (rhs_mx - cross(w_ii, I * w_ii));
    
    % dw = dCq * dq + Cq * ddq
    ddq_ii = Cq \ (dw_ii - dCq * dq_ii);
    ddqBase(:, ii) = ddq_ii;
    
end

dx = [drBase; dqBase; ddrBase; ddqBase];

end