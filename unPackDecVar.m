function [t, p, f, x] = unPackDecVar(z, pack)

nState = pack.nState;
nTime = pack.nTime;
tIdx = pack.tIdx;
pIdx = pack.pIdx;
fIdx = pack.fIdx;
xIdx = pack.xIdx;

% 开始解开包，得到相对简单的t,x
t = z(tIdx)';
x = z(xIdx);
x = reshape(x, nState, nTime);

% 比较复杂一点儿的就数p和f了
% 先得到p,f的优化变量的列向量形式
pCol_DV = z(pIdx);
fCol_DV = z(fIdx);

% 再将列向量进行填充得到原本的形式
p = zeros(1, pack.np);
f = zeros(1, pack.nf);

idx_p_DV = pack.idx_p_DV;
idx_f_DV = pack.idx_f_DV;

p(idx_p_DV) = pCol_DV;
f(idx_f_DV) = fCol_DV;

% p中的一些数据要与优化变量保持一致
p([13, 15, 17]) = p(11);                        p([25, 27]) = p(23);
p(28 + [13, 15, 17]) = p(28 + 11);              p(28 + [25, 27]) = p(28 + 23);
p(28 * 2 + [13, 15, 17]) = p(28 * 2 + 11);      p(28 * 2 + [25, 27]) = p(28 * 2 + 23);

% fz的第一个应该等于m*g
f(32 * 2 + 1) = 20 * 9.81;

end