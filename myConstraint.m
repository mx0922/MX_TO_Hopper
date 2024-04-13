function [c, ceq] = myConstraint(z, pack, dynamics, bndCst, pathCst)

global param

% 解构出约束要用的相关量
[t, p, f, x] = unPackDecVar(z, pack);

% 等式约束1：总时间 = T = param.T
% 不等式约束1：每一个phase的duration >= 0.1s(param.minT)
T = sum(t);
ceq1 = T - param.T;
c1 = reshape(param.minT - t, [], 1);

% 等式约束2：地形相关，由于假设为平地，因此不需要与x,y相关，直接使其两次落脚点的高度等于0即可
ceq2 = [p(28 * 2 + 11); p(28 * 2 + 23)]; % 脚的运动：z方向的两次接触的位置为0(平地)

% 利用得到的p,f进行三次多项式插值得到细分时刻的p和f
pNode = getFootMotionNodes(t, p, param.Nt);
fNode = getFootForceNodes(t, f, param.Nt);

% 利用细分时刻的p,f计算动力学及约束
u = [pNode; fNode];
tSpan = linspace(0, T, param.Nt);

% 得到边界约束与路径约束
[c_bnd, ceq_bnd] = bndCst(x(:, 1), x(:, end));

[c_path, ceq_path] = pathCst(tSpan, x, u);

% 动力学得到dx
dx = dynamics(tSpan, x, u);
[c_dyn, ceq_dyn] = getDynamicsCsts(x, dx, T);

% 整合
c = [c1; c_bnd; c_path; c_dyn];
ceq = [ceq1; ceq2; ceq_bnd; ceq_path; ceq_dyn];

end