% �Ƶ����Ȼ�����ģ�͵��˶�ѧ��ϵ
clear; close all; clc

i = sym([1; 0; 0]);
j = sym([0; 1; 0]);
k = sym([0; 0; 1]);

% rotation matrix - SO(3)
Rot_x = @(q)( [i, cos(q)*j + sin(q)*k, cos(q)*k - sin(q)*j] );
Rot_y = @(q)( [cos(q)*i - sin(q)*k, j, cos(q)*k + sin(q)*i] );
Rot_z = @(q)( [cos(q)*i + sin(q)*j, cos(q)*j - sin(q)*i, k] );

% base��λ���Լ�ŷ����
syms rx ry rz qx qy qz 'real'

rBase = [rx; ry; rz];
qBase = [qx; qy; qz];

% ���Կ�����roll-pitch-yaw(ZYX)
R = Rot_z(qz) * Rot_y(qy) * Rot_x(qx);

% ��С�ȳ��ȣ�l1, l2
syms l1 l2 'real'
% �Ȳ��Ĺؽڽǣ�hip(z,x,y -- q1, q2, q3), knee(y -- q4)
syms q1 q2 q3 q4 'real'
% hip�ؽ������base���ĵľ���
syms j1z 'real'

J1 = -j1z * k;

L1 = -l1 * k;

L2 = -l2 * k;

R1 = Rot_z(q1) * Rot_x(q2) * Rot_y(q3);
R2 = Rot_y(q4);

P_knee_base = J1 + R1 * L1;
P_ankle_base = J1 + R1 * (L1 + R2 * L2);

P_ankle_world = rBase + R * P_ankle_base;
P_knee_world  = rBase + R * P_knee_base;
P_hip_world   = rBase + R * J1;

pPoints = [P_hip_world, P_knee_world, P_ankle_world];

matlabFunction(...
    pPoints, ...
    'file', 'autoGen_legPointsPos.m', ...
    'vars', {'rx', 'ry', 'rz', 'qx', 'qy', 'qz', 'q1', 'q2', 'q3', 'q4', 'j1z', 'l1', 'l2'});


% �������base��λ����base����ϵ�µı�ʾ���Լ�������Ȳ��ؽڵ�jacobian����
qLeg = [q1; q2; q3; q4];
Jac_foot_base = jacobian(P_ankle_base, qLeg);

matlabFunction(...
    Jac_foot_base, ...
    'file', 'autoGen_JacFootBase.m', ...
    'vars', {'q1', 'q2', 'q3', 'q4', 'j1z', 'l1', 'l2'});

% ���˶�ѧ��������base��λ����base����ϵ�ı�ʾ
matlabFunction(...
    P_ankle_base, ...
    'file', 'autoGen_PosFootBase.m', ...
    'vars', {'q1', 'q2', 'q3', 'q4', 'j1z', 'l1', 'l2'});