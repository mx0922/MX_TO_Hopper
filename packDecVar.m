function [z, pack] = packDecVar(t, p, f, x)

nt = numel(t);  np = numel(p);      nf = numel(f);      nx = numel(x);

% 先转换为列向量
tCol = reshape(t, nt, 1);
pCol = reshape(p, np, 1);
fCol = reshape(f, nf, 1);
xCol = reshape(x, nx, 1);

%% 对p,f的处理
% 脚的位置、速度作为DecVar在p中的序号
idx_px_DV = [7 8 9 10 11 19 20 21 22 23];
idx_py_DV = [7 8 9 10 11 19 20 21 22 23] + 28;
idx_pz_DV = [7 8 9 10 11 19 20 21 22 23] + 28 * 2;

idx_p_DV = [idx_px_DV, idx_py_DV, idx_pz_DV];

% 从pCol中取出相应的DecVar
pCol_DV = pCol(idx_p_DV);

% p中优化变量的个数
n_p_DV = numel(pCol_DV);

% 接触力的大小及导数
idx_fx_DV = [3 4 5 6 15 16 17 18 27 28 29 30 31];
idx_fy_DV = [3 4 5 6 15 16 17 18 27 28 29 30 31] + 32;
idx_fz_DV = [3 4 5 6 15 16 17 18 27 28 29 30 31] + 32 * 2;

idx_f_DV = [idx_fx_DV, idx_fy_DV, idx_fz_DV];

% 从fCol中取出相应的DecVar
fCol_DV = fCol(idx_f_DV);

% f中优化变量的个数
n_f_DV = numel(fCol_DV);

%% 对t,x的处理 
nState = size(x, 1);
nTime = size(x, 2);

indz = reshape(nt + n_p_DV + n_f_DV + (1:nx), nState, nTime);

tIdx = 1:nt;
pIdx = (nt + 1):(nt + n_p_DV);
fIdx = (nt + n_p_DV + 1):(nt + n_p_DV + n_f_DV);
xIdx = indz(1:nState, :);

% 将对应的数据塞到对应序号的位置上
z = zeros(nt + n_p_DV + n_f_DV + nx, 1);
z(tIdx(:), 1) = tCol;
z(pIdx(:), 1) = pCol_DV;
z(fIdx(:), 1) = fCol_DV;
z(xIdx(:), 1) = xCol;

% pack
pack.nState = nState;
pack.nTime = nTime;
pack.tIdx = tIdx;
pack.pIdx = pIdx;
pack.fIdx = fIdx;
pack.xIdx = xIdx;

pack.idx_p_DV = idx_p_DV;
pack.idx_f_DV = idx_f_DV;

pack.np = np;
pack.nf = nf;

end