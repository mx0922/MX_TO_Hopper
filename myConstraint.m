function [c, ceq] = myConstraint(z, pack, dynamics, bndCst, pathCst)

global param

% �⹹��Լ��Ҫ�õ������
[t, p, f, x] = unPackDecVar(z, pack);

% ��ʽԼ��1����ʱ�� = T = param.T
% ����ʽԼ��1��ÿһ��phase��duration >= 0.1s(param.minT)
T = sum(t);
ceq1 = T - param.T;
c1 = reshape(param.minT - t, [], 1);

% ��ʽԼ��2��������أ����ڼ���Ϊƽ�أ���˲���Ҫ��x,y��أ�ֱ��ʹ��������ŵ�ĸ߶ȵ���0����
ceq2 = [p(28 * 2 + 11); p(28 * 2 + 23)]; % �ŵ��˶���z��������νӴ���λ��Ϊ0(ƽ��)

% ���õõ���p,f�������ζ���ʽ��ֵ�õ�ϸ��ʱ�̵�p��f
pNode = getFootMotionNodes(t, p, param.Nt);
fNode = getFootForceNodes(t, f, param.Nt);

% ����ϸ��ʱ�̵�p,f���㶯��ѧ��Լ��
u = [pNode; fNode];
tSpan = linspace(0, T, param.Nt);

% �õ��߽�Լ����·��Լ��
[c_bnd, ceq_bnd] = bndCst(x(:, 1), x(:, end));

[c_path, ceq_path] = pathCst(tSpan, x, u);

% ����ѧ�õ�dx
dx = dynamics(tSpan, x, u);
[c_dyn, ceq_dyn] = getDynamicsCsts(x, dx, T);

% ����
c = [c1; c_bnd; c_path; c_dyn];
ceq = [ceq1; ceq2; ceq_bnd; ceq_path; ceq_dyn];

end