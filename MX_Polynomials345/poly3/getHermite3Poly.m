function [p, v] = getHermite3Poly(S, t)

a = S.a;    b = S.b;    c = S.c;    d = S.d;

p = a + b * t +   c * t^2 +     d * t^3;
v =         b + 2 * c * t + 3 * d * t^2;

end