function pNode = getFootMotionNodes(tPhase, p, Nt)

T1 = tPhase(1);
T2 = tPhase(2);
T3 = tPhase(3);
T4 = tPhase(4);
T5 = tPhase(5);

% T1 stance phase
x11 = p(1);    dx11 = p(2);
x12 = p(3);    dx12 = p(4);
S1_x = getHermite3PolyCoeff(T1, x11, dx11, x12, dx12);

y11 = p(1 + 28);    dy11 = p(2 + 28);
y12 = p(2 + 28);    dy12 = p(2 + 28);
S1_y = getHermite3PolyCoeff(T1, y11, dy11, y12, dy12);

z11 = p(1 + 28 * 2);    dz11 = p(2 + 28 * 2);
z12 = p(3 + 28 * 2);    dz12 = p(4 + 28 * 2);
S1_z = getHermite3PolyCoeff(T1, z11, dz11, z12, dz12);

% T2 swing phase
x21 = p(5);     dx21 = p(6);
x22 = p(7);     dx22 = p(8);
x23 = p(9);     dx23 = p(10);
x24 = p(11);    dx24 = p(12);
S21_x = getHermite3PolyCoeff(T2/3, x21, dx21, x22, dx22);
S22_x = getHermite3PolyCoeff(T2/3, x22, dx22, x23, dx23);
S23_x = getHermite3PolyCoeff(T2/3, x23, dx23, x24, dx24);

y21 = p(5 + 28);     dy21 = p(6 + 28);
y22 = p(7 + 28);     dy22 = p(8 + 28);
y23 = p(9 + 28);     dy23 = p(10 + 28);
y24 = p(11 + 28);    dy24 = p(12 + 28);
S21_y = getHermite3PolyCoeff(T2/3, y21, dy21, y22, dy22);
S22_y = getHermite3PolyCoeff(T2/3, y22, dy22, y23, dy23);
S23_y = getHermite3PolyCoeff(T2/3, y23, dy23, y24, dy24);

z21 = p(5 + 28 * 2);     dz21 = p(6 + 28 * 2);
z22 = p(7 + 28 * 2);     dz22 = p(8 + 28 * 2);
z23 = p(9 + 28 * 2);     dz23 = p(10 + 28 * 2);
z24 = p(11 + 28 * 2);    dz24 = p(12 + 28 * 2);
S21_z = getHermite3PolyCoeff(T2/3, z21, dz21, z22, dz22);
S22_z = getHermite3PolyCoeff(T2/3, z22, dz22, z23, dz23);
S23_z = getHermite3PolyCoeff(T2/3, z23, dz23, z24, dz24);

% T3 stance phase
x31 = x24;      dx31 = dx24;
x32 = x24;      dx32 = dx24;
S3_x = getHermite3PolyCoeff(T3, x31, dx31, x32, dx32);

y31 = y24;      dy31 = dy24;
y32 = y24;      dy32 = dy24;
S3_y = getHermite3PolyCoeff(T3, y31, dy31, y32, dy32);

z31 = z24;      dz31 = dz24;
z32 = z24;      dz32 = dz24;
S3_z = getHermite3PolyCoeff(T3, z31, dz31, z32, dz32);

% T4 swing phase
x41 = p(17);        dx41 = p(18);
x42 = p(19);        dx42 = p(20);
x43 = p(21);        dx43 = p(22);
x44 = p(23);        dx44 = p(24);
S41_x = getHermite3PolyCoeff(T4/3, x41, dx41, x42, dx42);
S42_x = getHermite3PolyCoeff(T4/3, x42, dx42, x43, dx43);
S43_x = getHermite3PolyCoeff(T4/3, x43, dx43, x44, dx44);

y41 = p(17 + 28);        dy41 = p(18 + 28);
y42 = p(19 + 28);        dy42 = p(20 + 28);
y43 = p(21 + 28);        dy43 = p(22 + 28);
y44 = p(23 + 28);        dy44 = p(24 + 28);
S41_y = getHermite3PolyCoeff(T4/3, y41, dy41, y42, dy42);
S42_y = getHermite3PolyCoeff(T4/3, y42, dy42, y43, dy43);
S43_y = getHermite3PolyCoeff(T4/3, y43, dy43, y44, dy44);

z41 = p(17 + 28 * 2);        dz41 = p(18 + 28 * 2);
z42 = p(19 + 28 * 2);        dz42 = p(20 + 28 * 2);
z43 = p(21 + 28 * 2);        dz43 = p(22 + 28 * 2);
z44 = p(23 + 28 * 2);        dz44 = p(24 + 28 * 2);
S41_z = getHermite3PolyCoeff(T4/3, z41, dz41, z42, dz42);
S42_z = getHermite3PolyCoeff(T4/3, z42, dz42, z43, dz43);
S43_z = getHermite3PolyCoeff(T4/3, z43, dz43, z44, dz44);

% T5 stance phase
x51 = x44;      dx51 = dx44;
x52 = x44;      dx52 = dx44;
S5_x = getHermite3PolyCoeff(T5, x51, dx51, x52, dx52);

y51 = y44;      dy51 = dy44;
y52 = y44;      dy52 = dy44;
S5_y = getHermite3PolyCoeff(T5, y51, dy51, y52, dy52);

z51 = z44;      dz51 = dz44;
z52 = z44;      dz52 = dz44;
S5_z = getHermite3PolyCoeff(T5, z51, dz51, z52, dz52);

% 求整个时间区间内的脚的轨迹
T = T1 + T2 + T3 + T4 + T5;
dt = T / (Nt - 1);

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
        [pos_x(i), vel_x(i)] = getHermite3Poly(S1_x, t);
        [pos_y(i), vel_y(i)] = getHermite3Poly(S1_y, t);
        [pos_z(i), vel_z(i)] = getHermite3Poly(S1_z, t);
    elseif t > T1 + tol && t <= T1 + T2 + tol
        tt = t - T1;
        if tt <= T2/3 + tol
            [pos_x(i), vel_x(i)] = getHermite3Poly(S21_x, tt);
            [pos_y(i), vel_y(i)] = getHermite3Poly(S21_y, tt);
            [pos_z(i), vel_z(i)] = getHermite3Poly(S21_z, tt);
        elseif tt > T2/3 + tol && tt <= T2 * 2 / 3 + tol
            [pos_x(i), vel_x(i)] = getHermite3Poly(S22_x, tt - T2/3);
            [pos_y(i), vel_y(i)] = getHermite3Poly(S22_y, tt - T2/3);
            [pos_z(i), vel_z(i)] = getHermite3Poly(S22_z, tt - T2/3);
        else
            [pos_x(i), vel_x(i)] = getHermite3Poly(S23_x, tt - T2 * 2 /3);
            [pos_y(i), vel_y(i)] = getHermite3Poly(S23_y, tt - T2 * 2 /3);
            [pos_z(i), vel_z(i)] = getHermite3Poly(S23_z, tt - T2 * 2 /3);
        end
    elseif t > T1 + T2 + tol && t <= T1 + T2 + T3 + tol
        tt = t - T1 - T2;
        [pos_x(i), vel_x(i)] = getHermite3Poly(S3_x, tt);
        [pos_y(i), vel_y(i)] = getHermite3Poly(S3_y, tt);
        [pos_z(i), vel_z(i)] = getHermite3Poly(S3_z, tt);
    elseif t > T1 + T2 + T3 + tol && t <= T1 + T2 + T3 + T4 + tol
        tt = t - T1 - T2 - T3;
        if tt <= T4/3 + tol
            [pos_x(i), vel_x(i)] = getHermite3Poly(S41_x, tt);
            [pos_y(i), vel_y(i)] = getHermite3Poly(S41_y, tt);
            [pos_z(i), vel_z(i)] = getHermite3Poly(S41_z, tt);
        elseif tt > T4/3 + tol && tt <= T4 * 2 / 3 + tol
            [pos_x(i), vel_x(i)] = getHermite3Poly(S42_x, tt - T4/3);
            [pos_y(i), vel_y(i)] = getHermite3Poly(S42_y, tt - T4/3);
            [pos_z(i), vel_z(i)] = getHermite3Poly(S42_z, tt - T4/3);
        else
            [pos_x(i), vel_x(i)] = getHermite3Poly(S43_x, tt - T4 * 2 /3);
            [pos_y(i), vel_y(i)] = getHermite3Poly(S43_y, tt - T4 * 2 /3);
            [pos_z(i), vel_z(i)] = getHermite3Poly(S43_z, tt - T4 * 2 /3);
        end
        
    else
        tt = t - T1 - T2 - T3 - T4;
        [pos_x(i), vel_x(i)] = getHermite3Poly(S5_x, tt);
        [pos_y(i), vel_y(i)] = getHermite3Poly(S5_y, tt);
        [pos_z(i), vel_z(i)] = getHermite3Poly(S5_z, tt);
    end
    
end

pNode = [pos_x; pos_y; pos_z];

end