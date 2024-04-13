function [c, ceq] = pathConstraint(t, x, u, p)

% ������Ҫ�����˶�ѧԼ��:���ŵ�λ�ò��ܳ����Ǹ�cube
% ��ʵҲ���Խ��Ӵ�����ŵ�λ�ü��뵽·��Լ����

Nt = length(t);

rBase = x(1:3, :);
qBase = x(4:6, :);

P_foot = u(1:3, :);
F_foot = u(4:6, :);

%% 1.�ŵ��˶�ѧԼ��
p_hat = p.p_hat;
L_cube = p.L_cube;

deltaP = zeros(Nt, 1);
for ii = 1:Nt
    % ����ŷ���ǵõ�Base��World����ת����
    R = getRotMatFromBaseToWorld(qBase(:, ii));
    
    % �������base��λ����base����ϵ�µı�ʾ
    p_Base = R \ (P_foot(:, ii) - rBase(:, ii));
    
    % �õ���nominalλ�õľ����
    delta_p = p_Base - p_hat;
    
    % �ദ��һ�����õ���ģ��ƽ��
    deltaP(ii) = delta_p' * delta_p;
end

% deltaP <= L_cube^2
c1 = deltaP - L_cube^2;

%% �Ӵ�����ŵ�λ��
miu = p.miu;

fx = F_foot(1, :);
fy = F_foot(2, :);
fz = F_foot(3, :);

% fz >= 0
c2 = reshape(-fz, Nt, 1);

% Ħ��׶
c3 = [-fx - miu * fz, fx - miu * fz, -fy - miu * fz, fy - miu * fz];
c3 = reshape(c3, numel(c3), 1);

% p_foot_z >= 0
c4 = reshape(-P_foot(3, :), Nt, 1);

%% ����
c = [c1; c2; c3; c4];
ceq = [];

end