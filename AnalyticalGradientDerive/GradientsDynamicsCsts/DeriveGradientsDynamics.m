% 推导动力学相关约束的数值梯度，由于牵扯进了脚的运动与接触力，解析梯度求解比较繁杂，比较考验耐心

% 考虑到matlab的符号简化速度太慢，笔者也提供了利用maple求解动力学约束的解析梯度的程序，速度较快。

% Written by Meng Xiang in BIT on 2021/10/19. Copyrights Reserved.
clear; close all; clc

disp('Note: A little time-consuming, be patient...');

%% part 1: 脚的运动和接触力——前期准备
% px
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

% py
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

% pz
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

% 整合
pFoot = [px; py; pz];

disp('-->1-1.Foot motion obtained!');

%% force
% fx
syms a b c d 'real'
syms fx0 dfx0 fx1 dfx1 'real'
syms t T 'real'

x = a + b * t + c * t^2 + d * t^3;
dx = diff(x, t);

t = 0;
eq_x0 = subs(x) == fx0;
eq_dx0 = subs(dx) == dfx0;

t = T;
eq_x1 = subs(x) == fx1;
eq_dx1 = subs(dx) == dfx1;

S1 = solve([eq_x0, eq_dx0, eq_x1, eq_dx1], [a,b,c,d]);
      
a = S1.a;
b = S1.b;
c = S1.c;
d = S1.d;

syms t 'real'
fx = subs(x);

% fy
syms a b c d 'real'
syms fy0 dfy0 fy1 dfy1 'real'
syms t T 'real'

y = a + b * t + c * t^2 + d * t^3;
dy = diff(y, t);

t = 0;
eq_y0 = subs(y) == fy0;
eq_dy0 = subs(dy) == dfy0;

t = T;
eq_y1 = subs(y) == fy1;
eq_dy1 = subs(dy) == dfy1;

S2 = solve([eq_y0, eq_dy0, eq_y1, eq_dy1], [a,b,c,d]);
      
a = S2.a;
b = S2.b;
c = S2.c;
d = S2.d;

syms t 'real'
fy = subs(y);

% fz
syms a b c d 'real'
syms fz0 dfz0 fz1 dfz1 'real'
syms t T 'real'

z = a + b * t + c * t^2 + d * t^3;
dz = diff(z, t);

t = 0;
eq_z0 = subs(z) == fz0;
eq_dz0 = subs(dz) == dfz0;

t = T;
eq_z1 = subs(z) == fz1;
eq_dz1 = subs(dz) == dfz1;

S3 = solve([eq_z0, eq_dz0, eq_z1, eq_dz1], [a,b,c,d]);
      
a = S3.a;
b = S3.b;
c = S3.c;
d = S3.d;

syms t 'real'
fz = subs(z);

% 整合
fFoot = [fx; fy; fz];

disp('-->1-2. Foot force obtained!');

%% part 2: dynamics

syms rx ry rz qx qy qz drx dry drz dqx dqy dqz 'real'

i = sym([1; 0; 0]);
j = sym([0; 1; 0]);
k = sym([0; 0; 1]);

% rotation matrix - SO(3)
Rot_x = @(q)( [i, cos(q)*j + sin(q)*k, cos(q)*k - sin(q)*j] );
Rot_y = @(q)( [cos(q)*i - sin(q)*k, j, cos(q)*k + sin(q)*i] );
Rot_z = @(q)( [cos(q)*i + sin(q)*j, cos(q)*j - sin(q)*i, k] );

% 单刚体系统的动力学参数
syms m g Ix Iy Iz 'real'

G = g * k;
I = diag([Ix Iy Iz]);

rBase = [rx; ry; rz];
qBase = [qx; qy; qz];

drBase = [drx; dry; drz];
dqBase = [dqx; dqy; dqz];

ddrBase = fFoot / m - G;

% 欧拉角与角速度之间的转换关系
Cq = [Rot_z(qz)*Rot_y(qy)*i, Rot_z(qz)*j, k];

dCq = sym(zeros(3));

for ii = 1:3
    dCq(:, ii) = jacobian(Cq(:, ii), qBase) * dqBase;
end

rhs = cross(fFoot, rBase - pFoot);

w = Cq * dqBase;

dw = I \ (rhs - cross(w, I * w));

ddqBase = Cq \ (dw - dCq * dqBase);

disp('-->2. Dynamics equations obtained!');

%% part 3: 得到dynamics之后，利用四次多项式得到约束

% 2021.10.16注：
% 还将T,t用在四次多项式中是不对的，你细品
% 这里用T1，t1表示

syms rx1 ry1 rz1 qx1 qy1 qz1 drx1 dry1 drz1 dqx1 dqy1 dqz1 'real' 

% rx
syms a b c d e f 'real'
syms t1 T1 'real'

pos = a + b * t1 + c * t1^2 + d * t1^3 + e * t1^4;
vel = diff(pos, t1);
acc = diff(pos, t1, 2);

t1 = 0;
eq_p0 = subs(pos) == rx;
eq_v0 = subs(vel) == drx;
eq_a0 = subs(acc) == ddrBase(1);

t1 = T1;
eq_p1 = subs(pos) == rx1;
eq_v1 = subs(vel) == drx1;

S = solve([eq_p0, eq_v0, eq_a0, eq_p1, eq_v1], [a,b,c,d,e]);

a = S.a;
b = S.b;
c = S.c;
d = S.d;
e = S.e;

syms t1 'real'
ddrx = subs(acc);

% ry
syms a b c d e f 'real'
syms t1 T1 'real'

pos = a + b * t1 + c * t1^2 + d * t1^3 + e * t1^4;
vel = diff(pos, t1);
acc = diff(pos, t1, 2);

t1 = 0;
eq_p0 = subs(pos) == ry;
eq_v0 = subs(vel) == dry;
eq_a0 = subs(acc) == ddrBase(2);

t1 = T1;
eq_p1 = subs(pos) == ry1;
eq_v1 = subs(vel) == dry1;

S = solve([eq_p0, eq_v0, eq_a0, eq_p1, eq_v1], [a,b,c,d,e]);

a = S.a;
b = S.b;
c = S.c;
d = S.d;
e = S.e;

syms t1 'real'
ddry = subs(acc);

% rz
syms a b c d e f 'real'
syms t1 T1 'real'

pos = a + b * t1 + c * t1^2 + d * t1^3 + e * t1^4;
vel = diff(pos, t1);
acc = diff(pos, t1, 2);

t1 = 0;
eq_p0 = subs(pos) == rz;
eq_v0 = subs(vel) == drz;
eq_a0 = subs(acc) == ddrBase(3);

t1 = T1;
eq_p1 = subs(pos) == rz1;
eq_v1 = subs(vel) == drz1;

S = solve([eq_p0, eq_v0, eq_a0, eq_p1, eq_v1], [a,b,c,d,e]);

a = S.a;
b = S.b;
c = S.c;
d = S.d;
e = S.e;

syms t1 'real'
ddrz = subs(acc);

% qx
syms a b c d e f 'real'
syms t1 T1 'real'

assume(T1 > 0);

pos = a + b * t1 + c * t1^2 + d * t1^3 + e * t1^4;
vel = diff(pos, t1);
acc = diff(pos, t1, 2);

t1 = 0;
eq_p0 = subs(pos) == qx;
eq_v0 = subs(vel) == dqx;
eq_a0 = subs(acc) == ddqBase(1);

t1 = T1;
eq_p1 = subs(pos) == qx1;
eq_v1 = subs(vel) == dqx1;

S = solve([eq_p0, eq_v0, eq_a0, eq_p1, eq_v1], [a,b,c,d,e]);

a = S.a;
b = S.b;
c = S.c;
d = S.d;
e = S.e;

syms t1 'real'
ddqx = subs(acc);

% qy
syms a b c d e f 'real'
syms t1 T1 'real'

assume(T1 > 0);

pos = a + b * t1 + c * t1^2 + d * t1^3 + e * t1^4;
vel = diff(pos, t1);
acc = diff(pos, t1, 2);

t1 = 0;
eq_p0 = subs(pos) == qy;
eq_v0 = subs(vel) == dqy;
eq_a0 = subs(acc) == ddqBase(2);

t1 = T1;
eq_p1 = subs(pos) == qy1;
eq_v1 = subs(vel) == dqy1;

S = solve([eq_p0, eq_v0, eq_a0, eq_p1, eq_v1], [a,b,c,d,e]);

a = S.a;
b = S.b;
c = S.c;
d = S.d;
e = S.e;

syms t1 'real'
ddqy = subs(acc);

% qz
syms a b c d e f 'real'
syms t1 T1 'real'

assume(T1 > 0);

pos = a + b * t1 + c * t1^2 + d * t1^3 + e * t1^4;
vel = diff(pos, t1);
acc = diff(pos, t1, 2);

t1 = 0;
eq_p0 = subs(pos) == qz;
eq_v0 = subs(vel) == dqz;
eq_a0 = subs(acc) == ddqBase(3);

t1 = T1;
eq_p1 = subs(pos) == qz1;
eq_v1 = subs(vel) == dqz1;

S = solve([eq_p0, eq_v0, eq_a0, eq_p1, eq_v1], [a,b,c,d,e]);

a = S.a;
b = S.b;
c = S.c;
d = S.d;
e = S.e;

syms t1 'real'
ddqz = subs(acc);

disp('-->3. Dynamics acceleration obtained!');

%% 整合一下
ddx = [ddrx; ddry; ddrz; ddqx; ddqy; ddqz];

% 这里得到的函数体量太大，可以考虑用maple进行化简——利用maple得到的简化版函数体量仅为matlab生成的1/50，大大加快了运算速度

% matlabFunction(...
%     ddx, ...
%     'file', 'autoGen_dynamicsCst.m', ...
%     'vars', {'T', 'x0', 'dx0', 'x1', 'dx1', 'y0', 'dy0', 'y1', 'dy1', 'z0', 'dz0', 'z1', 'dz1', ...
%     'fx0', 'dfx0', 'fx1', 'dfx1', 'fy0', 'dfy0', 'fy1', 'dfy1', 'fz0', 'dfz0', 'fz1', 'dfz1', ...
%     'rx', 'ry', 'rz', 'qx', 'qy', 'qz', 'drx', 'dry', 'drz', 'dqx', 'dqy', 'dqz', ...
%     'rx1', 'ry1', 'rz1', 'qx1', 'qy1', 'qz1', 'drx1', 'dry1', 'drz1', 'dqx1', 'dqy1', 'dqz1', ...
%     'm', 'g', 'Ix', 'Iy', 'Iz', 't', 't1', 'T1'}, ...
%     'Optimize', false);

Jac_ddx = jacobian(ddx, [T, x0, dx0, x1, dx1, y0, dy0, y1, dy1, z0, dz0, z1, dz1, fx0, dfx0, fx1, dfx1, fy0, dfy0, fy1, dfy1, fz0, dfz0, fz1, dfz1, rx, ry, rz, qx, qy, qz, drx, dry, drz, dqx, dqy, dqz, rx1, ry1, rz1, qx1, qy1, qz1, drx1, dry1, drz1, dqx1, dqy1, dqz1, t, t1, T1]);

% matlabFunction(...
%     Jac_ddx, ...
%     'file', 'autoGen_dynamicsCstGrad.m', ...
%     'vars', {'T', 'x0', 'dx0', 'x1', 'dx1', 'y0', 'dy0', 'y1', 'dy1', 'z0', 'dz0', 'z1', 'dz1', ...
%     'fx0', 'dfx0', 'fx1', 'dfx1', 'fy0', 'dfy0', 'fy1', 'dfy1', 'fz0', 'dfz0', 'fz1', 'dfz1', ...
%     'rx', 'ry', 'rz', 'qx', 'qy', 'qz', 'drx', 'dry', 'drz', 'dqx', 'dqy', 'dqz', ...
%     'rx1', 'ry1', 'rz1', 'qx1', 'qy1', 'qz1', 'drx1', 'dry1', 'drz1', 'dqx1', 'dqy1', 'dqz1', ...
%     'm', 'g', 'Ix', 'Iy', 'Iz', 't', 't1', 'T1'}, ...
%     'Optimize', false);

