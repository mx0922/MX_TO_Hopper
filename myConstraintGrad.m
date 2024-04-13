function [c, ceq, cGrad, ceqGrad] = myConstraintGrad(z, pack, dynamics, bndCst, pathCst)

[c, ceq] = myConstraint(z, pack, dynamics, bndCst, pathCst);

% 解构成约束要用到的形式
[t, p, f, x] = unPackDecVar(z, pack);

nDecVar = numel(z); % nDecVar = 5 + 10 * 3 + 13 * 3 + 12 * 31

% 等式约束1：总时间 = T = 1.5s
% 不等式约束1：每一个phase的duration >= 0.1s
% T = sum(t);
% ceq1 = T - param.T;
ceq1Grad = zeros(1, nDecVar);
ceq1Grad(1:5) = 1;

% c1 = reshape(0.1 - t, [], 1);
c1Grad = zeros(numel(t), nDecVar);
c1Grad(1, 1) = -1;
c1Grad(2, 2) = -1;
c1Grad(3, 3) = -1;
c1Grad(4, 4) = -1;
c1Grad(5, 5) = -1;

% 利用得到的p,f进行多项式插值得到细分时刻的p和f
ceq2 = [p(28 * 2 + 11); p(28 * 2 + 23)]; % 脚的运动：z方向的两次接触的位置为0(平地)
ceq2Grad = zeros(numel(ceq2), nDecVar);
ceq2Grad(1, 5 + 20 + 5) = 1;
ceq2Grad(2, 5 + 20 + 10) = 1;

% 再得到边界约束与路径约束
[~, ceq_bnd] = bndCst(x(:, 1), x(:, end));

ceq_bndGrad = zeros(numel(ceq_bnd), nDecVar);
ceq_bndGrad(1, 74 + 1) = -1;
ceq_bndGrad(2, 74 + 2) = -1;
ceq_bndGrad(3, 74 + 3) = -1;
ceq_bndGrad(4, 74 + 4) = -1;
ceq_bndGrad(5, 74 + 5) = -1;
ceq_bndGrad(6, 74 + 6) = -1;

ceq_bndGrad( 7, 434 + 1) = -1;
ceq_bndGrad( 8, 434 + 2) = -1;
ceq_bndGrad( 9, 434 + 3) = -1;
ceq_bndGrad(10, 434 + 4) = -1;
ceq_bndGrad(11, 434 + 5) = -1;
ceq_bndGrad(12, 434 + 6) = -1;

% pNode = getFootMotionNodes(t, p, param.Nt);
% fNode = getFootForceNodes(t, f, param.Nt);

% pNode = getFootMotionNodes(t, p, 1501);
% fNode = getFootForceNodes(t, f, 1501);

% % 利用细分时刻的p,f计算动力学及约束
% u = [pNode; fNode];
% tSpan = linspace(0, T, param.Nt);

% % 比较复杂的就数pathCst以及动力学约束了
% [c_path, ceq_path] = pathCst(tSpan, x, u);

[c_pathGrad, ceq_pathGrad] = getPathCstGrad(z, pack);

% % 动力学得到dx
% dx = dynamics(tSpan, x, u);
% [c_dyn, ceq_dyn] = getDynamicsCsts(x, dx, T);

[~, ~, c_dynGrad, ceq_dynGrad, ~] = getDynamicsCstGrad(z, pack);

% [c_dyn, ceq_dyn, c_dynGrad, ceq_dynGrad, dx_New] = getDynamicsCstGrad_New_New(z, pack);

% c =         [c1;     c_bnd;     c_path;     c_dyn];
% ceq =         [ceq1;     ceq2;     ceq_bnd;     ceq_path;     ceq_dyn];

cGrad = [c1Grad; c_pathGrad; c_dynGrad]';
ceqGrad = [ceq1Grad; ceq2Grad; ceq_bndGrad; ceq_pathGrad; ceq_dynGrad]';

end