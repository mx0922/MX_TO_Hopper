% ��Ҫ�Ƶ��ŵ��˶���Ӵ����Ľ����ݶȡ���path constraints�Ľ����ݶ�

clear; close all; clc

%% p&f��x,y,z��ά�������ζ���ʽ��ʾ
% px, fx
syms a b c d 'real'
syms x0 dx0 x1 dx1 'real'
syms t T 'real'

x = a + b * t + c * t^2 + d * t^3;
dx = diff(x, t);

t = 0;
eq_x0 = subs(x) == x0;
eq_dx0 = subs(dx) == dx0;

t = T;
eq_x1 = subs(x) == x1;
eq_dx1 = subs(dx) == dx1;

S1 = solve([eq_x0, eq_dx0, eq_x1, eq_dx1], [a,b,c,d]);
      
a = S1.a;
b = S1.b;
c = S1.c;
d = S1.d;

syms t 'real'
px = subs(x);

% % ��ô���������У�����
% syms t y0 dy0 y1 dy1 'real'
% py = subs(x, [x0 dx0 x1 dx1], [y0 dy0 y1 dy1]);
% 
% syms z0 dz0 z1 dz1 'real'
% pz = subs(x, [x0 dx0 x1 dx1], [z0 dz0 z1 dz1]);

% py,fy
syms a b c d 'real'
syms y0 dy0 y1 dy1 'real'
syms t T 'real'

y = a + b * t + c * t^2 + d * t^3;
dy = diff(y, t);

t = 0;
eq_y0 = subs(y) == y0;
eq_dy0 = subs(dy) == dy0;

t = T;
eq_y1 = subs(y) == y1;
eq_dy1 = subs(dy) == dy1;

S2 = solve([eq_y0, eq_dy0, eq_y1, eq_dy1], [a,b,c,d]);
      
a = S2.a;
b = S2.b;
c = S2.c;
d = S2.d;

syms t 'real'
py = subs(y);

% pz,fz
syms a b c d 'real'
syms z0 dz0 z1 dz1 'real'
syms t T 'real'

z = a + b * t + c * t^2 + d * t^3;
dz = diff(z, t);

t = 0;
eq_z0 = subs(z) == z0;
eq_dz0 = subs(dz) == dz0;

t = T;
eq_z1 = subs(z) == z1;
eq_dz1 = subs(dz) == dz1;

S3 = solve([eq_z0, eq_dz0, eq_z1, eq_dz1], [a,b,c,d]);
      
a = S3.a;
b = S3.b;
c = S3.c;
d = S3.d;

syms t 'real'
pz = subs(z);

%% ���õõ���p&f������Լ��
% 1.�ŵ��˶����˶�ѧԼ�� pFoot ��һ����p_hatΪ���ģ��߳�Ϊ2*L_cube��cube��
syms rx ry rz qx qy qz 'real'

i = sym([1; 0; 0]);
j = sym([0; 1; 0]);
k = sym([0; 0; 1]);

% rotation matrix - SO(3)
Rot_x = @(q)( [i, cos(q)*j + sin(q)*k, cos(q)*k - sin(q)*j] );
Rot_y = @(q)( [cos(q)*i - sin(q)*k, j, cos(q)*k + sin(q)*i] );
Rot_z = @(q)( [cos(q)*i + sin(q)*j, cos(q)*j - sin(q)*i, k] );

% ��ת����R
R = Rot_z(qz) * Rot_y(qy) * Rot_x(qx);

pFoot = [px; py; pz];
pCOM = [rx; ry; rz];

% �������base��λ����base����ϵ�µı�ʾpBase
% pFoot = pCOM + R * pBase;
pBase = R \ (pFoot - pCOM);

% pBase��һ���ο�λ��
syms p_hat_x p_hat_y p_hat_z 'real'
p_hat = [p_hat_x, p_hat_y, p_hat_z]';

% ʵ��λ����ο�λ�õĲ�
delta_p = pBase - p_hat;

% mx = (delta_p' * delta_p)^(1/2) - 0.2;

% û�в���winkler������2��������ʽ���б�ʾ�����õ���ģ��ƽ��
syms L_cube 'real'
kin = (delta_p' * delta_p) - L_cube^2;

Jac_kin = jacobian(kin, [T, x0, dx0, x1, dx1, y0, dy0, y1, dy1, z0, dz0, z1, dz1, rx, ry, rz, qx, qy, qz, t]);

% matlabFunction(...
%     Jac_kin, ...
%     'file', 'autoGen_getJacobianKinCst.m', ...
%     'vars', {'T', 'x0', 'dx0', 'x1', 'dx1', 'y0', 'dy0', 'y1', 'dy1', 'z0', 'dz0', 'z1', 'dz1', ...
%     'rx', 'ry', 'rz', 'qx', 'qy', 'qz', 't', 'p_hat_x', 'p_hat_y', 'p_hat_z', 'L_cube'}, ...
%     'Optimize', false);

% ��ʱԼ625s������������С��Լ1/10
tic;
matlabFunction(...
    Jac_kin, ...
    'file', 'autoGen_getJacobianKinCst.m', ...
    'vars', {'T', 'x0', 'dx0', 'x1', 'dx1', 'y0', 'dy0', 'y1', 'dy1', 'z0', 'dz0', 'z1', 'dz1', ...
    'rx', 'ry', 'rz', 'qx', 'qy', 'qz', 't', 'p_hat_x', 'p_hat_y', 'p_hat_z', 'L_cube'});
toc;

% ��2�͵�3��Լ����pz >= 0 , fz >= 0
Jac_pz = jacobian(-pz, [T, z0, dz0, z1, dz1, t]);
matlabFunction(...
    Jac_pz, ...
    'file', 'autoGen_getJacobianFzCst.m', ...
    'vars', {T, z0, dz0, z1, dz1, t});

% ��4��Լ����Ħ��׶Լ��
syms miu 'real' % Ħ��ϵ��

% -miu * fz <= fx <= miu * fz (yͬ��)
CWC = [-px - miu * pz; px - miu * pz; -py - miu * pz; py - miu * pz];

Jac_CWC = jacobian(CWC, [T, x0, dx0, x1, dx1, y0, dy0, y1, dy1, z0, dz0, z1, dz1, t]);
matlabFunction(...
    Jac_CWC, ...
    'file', 'autoGen_getJacobianCWCCst.m', ...
    'vars', {T, x0, dx0, x1, dx1, y0, dy0, y1, dy1, z0, dz0, z1, dz1, miu, t});