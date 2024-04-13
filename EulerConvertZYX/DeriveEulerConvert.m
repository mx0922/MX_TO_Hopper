% �Ƶ���ŷ�����ٶ�ת������ת�ٶȵ�ת������C���䵼��dC
% ���ڶ���ѧ������dq, ddq��w,dw֮���ת����ϵ��
clear; close all; clc

%% necessary preparations
% unit vectors - R(3)
i = sym([1; 0; 0]);
j = sym([0; 1; 0]);
k = sym([0; 0; 1]);

% rotation matrix - SO(3)
Rot_x = @(q)( [i, cos(q)*j + sin(q)*k, cos(q)*k - sin(q)*j] );
Rot_y = @(q)( [cos(q)*i - sin(q)*k, j, cos(q)*k + sin(q)*i] );
Rot_z = @(q)( [cos(q)*i + sin(q)*j, cos(q)*j - sin(q)*i, k] );

syms x y z dx dy dz 'real'

C = [Rot_z(z)*Rot_y(y)*i, Rot_z(z)*j, k];

q = [x; y; z];
dq = [dx; dy; dz];

dC = sym(zeros(3));

for ii = 1:3
    dC(:, ii) = jacobian(C(:, ii), q) * dq;
end
