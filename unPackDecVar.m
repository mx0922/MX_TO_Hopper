function [t, p, f, x] = unPackDecVar(z, pack)

nState = pack.nState;
nTime = pack.nTime;
tIdx = pack.tIdx;
pIdx = pack.pIdx;
fIdx = pack.fIdx;
xIdx = pack.xIdx;

% ��ʼ�⿪�����õ���Լ򵥵�t,x
t = z(tIdx)';
x = z(xIdx);
x = reshape(x, nState, nTime);

% �Ƚϸ���һ����ľ���p��f��
% �ȵõ�p,f���Ż���������������ʽ
pCol_DV = z(pIdx);
fCol_DV = z(fIdx);

% �ٽ��������������õ�ԭ������ʽ
p = zeros(1, pack.np);
f = zeros(1, pack.nf);

idx_p_DV = pack.idx_p_DV;
idx_f_DV = pack.idx_f_DV;

p(idx_p_DV) = pCol_DV;
f(idx_f_DV) = fCol_DV;

% p�е�һЩ����Ҫ���Ż���������һ��
p([13, 15, 17]) = p(11);                        p([25, 27]) = p(23);
p(28 + [13, 15, 17]) = p(28 + 11);              p(28 + [25, 27]) = p(28 + 23);
p(28 * 2 + [13, 15, 17]) = p(28 * 2 + 11);      p(28 * 2 + [25, 27]) = p(28 * 2 + 23);

% fz�ĵ�һ��Ӧ�õ���m*g
f(32 * 2 + 1) = 20 * 9.81;

end