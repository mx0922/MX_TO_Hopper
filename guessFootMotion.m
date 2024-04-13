function motion_foot = guessFootMotion()

% 接下来应该思考如何将guess值与提供的已知量结合起来自动对guess值进行赋值，而不是像目前这样手动来给
% 不过手动给也有手动给的好处：可以根据期望的运动给出一些比较合理的猜测值，更有利于求解器快速求解

%% pos and vel of X direction
% T1 stance phase
x11 = 0;    dx11 = 0;
x12 = 0;    dx12 = 0;

% T2 swing phase
x21 = x12;      dx21 = dx12;
x22 = 0.08;     dx22 = 1.0;
x23 = 0.17;     dx23 = 1.0;
x24 = 0.25;     dx24 = 0;

% T3 stance phase
x31 = x24;      dx31 = dx24;
x32 = x24;      dx32 = dx24;

% T4 swing phase
x41 = x32;      dx41 = dx32;
x42 = 0.33;     dx42 = 1.2;
x43 = 0.42;     dx43 = 1.0;
x44 = 0.5;      dx44 = 0;

% T5 stance phase
x51 = x44;      dx51 = dx44;
x52 = x44;      dx52 = dx44;

foot_x = [x11, dx11, x12, dx12, x21, dx21, x22, dx22, x23, dx23, x24, dx24, x31, dx31, x32, dx32, x41, dx41, x42, dx42, x43, dx43, x44, dx44, x51, dx51, x52, dx52];
% motion_foot.x = reshape(foot_x, 2, []);

%% pos and vel of Y direction
% T1 stance phase
y11 = 0;    dy11 = 0;
y12 = 0;    dy12 = 0;

% T2 swing phase
y21 = y12;      dy21 = dy12;
y22 = 0.01;     dy22 = 0.25;
y23 = 0.03;     dy23 = 0.25;
y24 = 0.05;     dy24 = 0;

% T3 stance phase
y31 = y24;      dy31 = dy24;
y32 = y24;      dy32 = dy24;

% T4 swing phase
y41 = y32;      dy41 = dy32;
y42 = 0.06;     dy42 = 0.25;
y43 = 0.08;     dy43 = 0.25;
y44 = 0.1;      dy44 = 0;

% T5 stance phase
y51 = y44;      dy51 = dy44;
y52 = y44;      dy52 = dy44;

foot_y = [y11, dy11, y12, dy12, y21, dy21, y22, dy22, y23, dy23, y24, dy24, y31, dy31, y32, dy32, y41, dy41, y42, dy42, y43, dy43, y44, dy44, y51, dy51, y52, dy52];
% motion_foot.y = reshape(foot_y, 2, []);

%% pos and vel of Z direction
% T1 stance phase
z11 = 0;    dz11 = 0;
z12 = 0;    dz12 = 0;

% T2 swing phase
z21 = z12;      dz21 = dz12;
z22 = 0.04;     dz22 = 0.5;
z23 = 0.04;     dz23 = -0.5;
z24 = 0;        dz24 = 0;

% T3 stance phase
z31 = z24;      dz31 = dz24;
z32 = z24;      dz32 = dz24;

% T4 swing phase
z41 = z32;      dz41 = dz32;
z42 = 0.04;     dz42 = 0.5;
z43 = 0.04;     dz43 = -0.5;
z44 = 0;        dz44 = 0;

% T5 stance phase
z51 = z44;      dz51 = dz44;
z52 = z44;      dz52 = dz44;

foot_z = [z11, dz11, z12, dz12, z21, dz21, z22, dz22, z23, dz23, z24, dz24, z31, dz31, z32, dz32, z41, dz41, z42, dz42, z43, dz43, z44, dz44, z51, dz51, z52, dz52];
% motion_foot.z = reshape(foot_z, 2, []);

%% 整合
motion_foot = [foot_x, foot_y, foot_z];

end