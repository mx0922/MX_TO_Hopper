clear; close all; clc

syms a b c d 'real'
syms p0 v0 p1 v1 'real'
syms t T 'real'

pos = a + b * t + c * t^2 + d * t^3;
vel = diff(pos, t);
acc = diff(pos, t, 2);

t = 0;
eq_p0 = subs(pos) == p0;
eq_v0 = subs(vel) == v0;

t = T;
eq_p1 = subs(pos) == p1;
eq_v1 = subs(vel) == v1;

S = solve([eq_p0, eq_v0, eq_p1, eq_v1], [a,b,c,d]);
      
a = S.a;
b = S.b;
c = S.c;
d = S.d;

% matlabFunction(...
%     a, b, c, d, ...
%     'file', 'autoGen_hermitePolyCoeff.m', ...
%     'vars', {T, p0, v0, p1, v1});

syms t 'real'
jac_p = jacobian(subs(pos), [T, p0, v0, p1, v1]);