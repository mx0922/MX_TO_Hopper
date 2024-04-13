function fNode = getFootForceNodes(tPhase, f, Nt)

T1 = tPhase(1);
T2 = tPhase(2);
T3 = tPhase(3);
T4 = tPhase(4);
T5 = tPhase(5);

% T1 stance phase
x11 = f(1);         dx11 = f(2);
x12 = f(3);         dx12 = f(4);
x13 = f(5);         dx13 = f(6);
x14 = f(7);         dx14 = f(8);
S11_x = getHermite3PolyCoeff(T1/3, x11, dx11, x12, dx12);
S12_x = getHermite3PolyCoeff(T1/3, x12, dx12, x13, dx13);
S13_x = getHermite3PolyCoeff(T1/3, x13, dx13, x14, dx14);

y11 = f(1 + 32);         dy11 = f(2 + 32);
y12 = f(3 + 32);         dy12 = f(4 + 32);
y13 = f(5 + 32);         dy13 = f(6 + 32);
y14 = f(7 + 32);         dy14 = f(8 + 32);
S11_y = getHermite3PolyCoeff(T1/3, y11, dy11, y12, dy12);
S12_y = getHermite3PolyCoeff(T1/3, y12, dy12, y13, dy13);
S13_y = getHermite3PolyCoeff(T1/3, y13, dy13, y14, dy14);

z11 = f(1 + 32 * 2);         dz11 = f(2 + 32 * 2);
z12 = f(3 + 32 * 2);         dz12 = f(4 + 32 * 2);
z13 = f(5 + 32 * 2);         dz13 = f(6 + 32 * 2);
z14 = f(7 + 32 * 2);         dz14 = f(8 + 32 * 2);
S11_z = getHermite3PolyCoeff(T1/3, z11, dz11, z12, dz12);
S12_z = getHermite3PolyCoeff(T1/3, z12, dz12, z13, dz13);
S13_z = getHermite3PolyCoeff(T1/3, z13, dz13, z14, dz14);

% T2 swing phase
x21 = x14;    dx21 = dx14;
x22 = x14;    dx22 = dx14;
S2_x = getHermite3PolyCoeff(T2, x21, dx21, x22, dx22);

y21 = y14;    dy21 = dy14;
y22 = y14;    dy22 = dy14;
S2_y = getHermite3PolyCoeff(T2, y21, dy21, y22, dy22);

z21 = z14;    dz21 = dz14;
z22 = z14;    dz22 = dz14;
S2_z = getHermite3PolyCoeff(T2, z21, dz21, z22, dz22);

% T3 stance phase
x31 = f(13);        dx31 = f(14);
x32 = f(15);        dx32 = f(16);
x33 = f(17);        dx33 = f(18);
x34 = f(19);        dx34 = f(20);
S31_x = getHermite3PolyCoeff(T3/3, x31, dx31, x32, dx32);
S32_x = getHermite3PolyCoeff(T3/3, x32, dx32, x33, dx33);
S33_x = getHermite3PolyCoeff(T3/3, x33, dx33, x34, dx34);

y31 = f(13 + 32);        dy31 = f(14 + 32);
y32 = f(15 + 32);        dy32 = f(16 + 32);
y33 = f(17 + 32);        dy33 = f(18 + 32);
y34 = f(19 + 32);        dy34 = f(20 + 32);
S31_y = getHermite3PolyCoeff(T3/3, y31, dy31, y32, dy32);
S32_y = getHermite3PolyCoeff(T3/3, y32, dy32, y33, dy33);
S33_y = getHermite3PolyCoeff(T3/3, y33, dy33, y34, dy34);

z31 = f(13 + 32 * 2);        dz31 = f(14 + 32 * 2);
z32 = f(15 + 32 * 2);        dz32 = f(16 + 32 * 2);
z33 = f(17 + 32 * 2);        dz33 = f(18 + 32 * 2);
z34 = f(19 + 32 * 2);        dz34 = f(20 + 32 * 2);
S31_z = getHermite3PolyCoeff(T3/3, z31, dz31, z32, dz32);
S32_z = getHermite3PolyCoeff(T3/3, z32, dz32, z33, dz33);
S33_z = getHermite3PolyCoeff(T3/3, z33, dz33, z34, dz34);

% T4 swing phase
x41 = x34;      dx41 = dx34;
x42 = x34;      dx42 = dx34;
S4_x = getHermite3PolyCoeff(T4, x41, dx41, x42, dx42);

y41 = y34;      dy41 = dy34;
y42 = y34;      dy42 = dy34;
S4_y = getHermite3PolyCoeff(T4, y41, dy41, y42, dy42);

z41 = z34;      dz41 = dz34;
z42 = z34;      dz42 = dz34;
S4_z = getHermite3PolyCoeff(T4, z41, dz41, z42, dz42);

% T5 stance phase
x51 = f(25);        dx51 = f(26);
x52 = f(27);        dx52 = f(28);
x53 = f(29);        dx53 = f(30);
x54 = f(31);        dx54 = f(32);
S51_x = getHermite3PolyCoeff(T5/3, x51, dx51, x52, dx52);
S52_x = getHermite3PolyCoeff(T5/3, x52, dx52, x53, dx53);
S53_x = getHermite3PolyCoeff(T5/3, x53, dx53, x54, dx54);

y51 = f(25 + 32);        dy51 = f(26 + 32);
y52 = f(27 + 32);        dy52 = f(28 + 32);
y53 = f(29 + 32);        dy53 = f(30 + 32);
y54 = f(31 + 32);        dy54 = f(32 + 32);
S51_y = getHermite3PolyCoeff(T5/3, y51, dy51, y52, dy52);
S52_y = getHermite3PolyCoeff(T5/3, y52, dy52, y53, dy53);
S53_y = getHermite3PolyCoeff(T5/3, y53, dy53, y54, dy54);

z51 = f(25 + 32 * 2);        dz51 = f(26 + 32 * 2);
z52 = f(27 + 32 * 2);        dz52 = f(28 + 32 * 2);
z53 = f(29 + 32 * 2);        dz53 = f(30 + 32 * 2);
z54 = f(31 + 32 * 2);        dz54 = f(32 + 32 * 2);
S51_z = getHermite3PolyCoeff(T5/3, z51, dz51, z52, dz52);
S52_z = getHermite3PolyCoeff(T5/3, z52, dz52, z53, dz53);
S53_z = getHermite3PolyCoeff(T5/3, z53, dz53, z54, dz54);

% 求整个时间区间内的脚的轨迹
T = T1 + T2 + T3 + T4 + T5;
dt = T / (Nt - 1);
% tSpan = linspace(0, T, Nt);
% 
% Nt = length(tSpan);

pos_x = zeros(1, Nt);
vel_x = zeros(1, Nt);

pos_y = zeros(1, Nt);
vel_y = zeros(1, Nt);

pos_z = zeros(1, Nt);
vel_z = zeros(1, Nt);

tol = 1e-6;

for i = 1:Nt
    t = (i - 1) * dt;
    
    if t <= T1 + tol
        tt = t;
        if tt <= T1/3 + tol
            [pos_x(i), vel_x(i)] = getHermite3Poly(S11_x, tt);
            [pos_y(i), vel_y(i)] = getHermite3Poly(S11_y, tt);
            [pos_z(i), vel_z(i)] = getHermite3Poly(S11_z, tt);
        elseif tt > T1/3 + tol && tt <= T1 * 2 / 3 + tol
            [pos_x(i), vel_x(i)] = getHermite3Poly(S12_x, tt - T1/3);
            [pos_y(i), vel_y(i)] = getHermite3Poly(S12_y, tt - T1/3);
            [pos_z(i), vel_z(i)] = getHermite3Poly(S12_z, tt - T1/3);
        else
            [pos_x(i), vel_x(i)] = getHermite3Poly(S13_x, tt - T1 * 2 /3);
            [pos_y(i), vel_y(i)] = getHermite3Poly(S13_y, tt - T1 * 2 /3);
            [pos_z(i), vel_z(i)] = getHermite3Poly(S13_z, tt - T1 * 2 /3);
        end        
        
    elseif t > T1 + tol && t <= T1 + T2 + tol
        tt = t - T1;
        [pos_x(i), vel_x(i)] = getHermite3Poly(S2_x, tt);
        [pos_y(i), vel_y(i)] = getHermite3Poly(S2_y, tt);
        [pos_z(i), vel_z(i)] = getHermite3Poly(S2_z, tt);
        
    elseif t > T1 + T2 + tol && t <= T1 + T2 + T3 + tol
        tt = t - T1 - T2;        
        if tt <= T3/3 + tol
            [pos_x(i), vel_x(i)] = getHermite3Poly(S31_x, tt);
            [pos_y(i), vel_y(i)] = getHermite3Poly(S31_y, tt);
            [pos_z(i), vel_z(i)] = getHermite3Poly(S31_z, tt);
        elseif tt > T3/3 + tol && tt <= T3 * 2 / 3 + tol
            [pos_x(i), vel_x(i)] = getHermite3Poly(S32_x, tt - T3/3);
            [pos_y(i), vel_y(i)] = getHermite3Poly(S32_y, tt - T3/3);
            [pos_z(i), vel_z(i)] = getHermite3Poly(S32_z, tt - T3/3);
        else
            [pos_x(i), vel_x(i)] = getHermite3Poly(S33_x, tt - T3 * 2 /3);
            [pos_y(i), vel_y(i)] = getHermite3Poly(S33_y, tt - T3 * 2 /3);
            [pos_z(i), vel_z(i)] = getHermite3Poly(S33_z, tt - T3 * 2 /3);
        end        
        
    elseif t > T1 + T2 + T3 + tol && t <= T1 + T2 + T3 + T4 + tol
        tt = t - T1 - T2 - T3;
        [pos_x(i), vel_x(i)] = getHermite3Poly(S4_x, tt);
        [pos_y(i), vel_y(i)] = getHermite3Poly(S4_y, tt);
        [pos_z(i), vel_z(i)] = getHermite3Poly(S4_z, tt);
        
    else
        tt = t - T1 - T2 - T3 - T4;
        if tt <= T5/3 + tol
            [pos_x(i), vel_x(i)] = getHermite3Poly(S51_x, tt);
            [pos_y(i), vel_y(i)] = getHermite3Poly(S51_y, tt);
            [pos_z(i), vel_z(i)] = getHermite3Poly(S51_z, tt);
        elseif tt > T5/3 + tol && tt <= T5 * 2 / 3 + tol
            [pos_x(i), vel_x(i)] = getHermite3Poly(S52_x, tt - T5/3);
            [pos_y(i), vel_y(i)] = getHermite3Poly(S52_y, tt - T5/3);
            [pos_z(i), vel_z(i)] = getHermite3Poly(S52_z, tt - T5/3);
        else
            [pos_x(i), vel_x(i)] = getHermite3Poly(S53_x, tt - T5 * 2 /3);
            [pos_y(i), vel_y(i)] = getHermite3Poly(S53_y, tt - T5 * 2 /3);
            [pos_z(i), vel_z(i)] = getHermite3Poly(S53_z, tt - T5 * 2 /3);
        end
    end
    
end

fNode = [pos_x; pos_y; pos_z];

end