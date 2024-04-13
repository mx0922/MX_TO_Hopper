function force_foot = guessFootForce(p)

m = p.m;
g = p.g;

%% value and derivative of force of X direction
% T1 stance phase
x11 = 0;        dx11 = 0;
x12 = 20;       dx12 = 0.5;
x13 = 30;       dx13 = -0.5;
x14 = 0;        dx14 = 0;

% T2 swing phase
x21 = x14;    dx21 = dx14;
x22 = x14;    dx22 = dx14;

% T3 stance phase
x31 = x22;      dx31 = dx22;
x32 = 20;       dx32 = 0.5;
x33 = 30;       dx33 = -0.5;
x34 = 0;        dx34 = 0;

% T4 swing phase
x41 = x34;      dx41 = dx34;
x42 = x34;      dx42 = dx34;

% T5 stance phase
x51 = x42;      dx51 = dx42;
x52 = 10;       dx52 = 0.5;
x53 = -10;      dx53 = 0.5;
x54 = 0;        dx54 = 0;

force_x = [x11, dx11, x12, dx12, x13, dx13, x14, dx14, x21, dx21, x22, dx22, x31, dx31, x32, dx32, x33, dx33, x34, dx34, x41, dx41, x42, dx42, x51, dx51, x52, dx52, x53, dx53, x54, dx54];
% force_foot.x = reshape(force_x, 2, []);

%% value and derivative of force of Y direction
% T1 stance phase
y11 = 0;        dy11 = 0;
y12 = 5;        dy12 = 0.15;
y13 = 10;       dy13 = -0.15;
y14 = 0;        dy14 = 0;

% T2 swing phase
y21 = y14;    dy21 = dy14;
y22 = y14;    dy22 = dy14;

% T3 stance phase
y31 = y22;      dy31 = dy22;
y32 = 5;        dy32 = 0.15;
y33 = 10;       dy33 = -0.15;
y34 = 0;        dy34 = 0;

% T4 swing phase
y41 = y34;      dy41 = dy34;
y42 = y34;      dy42 = dy34;

% T5 stance phase
y51 = y42;      dy51 = dy42;
y52 = 7;        dy52 = 0.2;
y53 = -7;       dy53 = 0.2;
y54 = 0;        dy54 = 0;

force_y = [y11, dy11, y12, dy12, y13, dy13, y14, dy14, y21, dy21, y22, dy22, y31, dy31, y32, dy32, y33, dy33, y34, dy34, y41, dy41, y42, dy42, y51, dy51, y52, dy52, y53, dy53, y54, dy54];
% force_foot.y = reshape(force_y, 2, []);

%% value and derivative of force of Z direction
% T1 stance phase
z11 = m * g;            dz11 = 0;
z12 = 0.7 * m * g;      dz12 = 0.5;
z13 = 1.2 * m * g;      dz13 = -0.5;
z14 = 0;                dz14 = 0;

% T2 swing phase
z21 = z14;    dz21 = dz14;
z22 = z14;    dz22 = dz14;

% T3 stance phase
z31 = z22;              dz31 = dz22;
z32 = 0.7 * m * g;      dz32 = 0.5;
z33 = 1.1 * m * g;      dz33 = -0.5;
z34 = 0;                dz34 = 0;

% T4 swing phase
z41 = z34;      dz41 = dz34;
z42 = z34;      dz42 = dz34;

% T5 stance phase
z51 = z42;                  dz51 = dz42;
z52 = 0.8 * m * g;          dz52 = 0.5;
z53 = 1.2 * m * g;          dz53 = -0.5;
z54 = m * g;                dz54 = 0;

force_z = [z11, dz11, z12, dz12, z13, dz13, z14, dz14, z21, dz21, z22, dz22, z31, dz31, z32, dz32, z33, dz33, z34, dz34, z41, dz41, z42, dz42, z51, dz51, z52, dz52, z53, dz53, z54, dz54];
% force_foot.z = reshape(force_z, 2, []);

%% ÕûºÏ
force_foot = [force_x, force_y, force_z];

end