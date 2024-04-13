function [pos, vel, acc] = get4orderPoly(S, t)

a = S.a;    b = S.b;    c = S.c;    d = S.d;    e = S.e;

pos = a + b .* t +   c .* t.^2 +     d .* t.^3 +      e .* t.^4;
vel =         b + 2 .* c .* t + 3 .* d .* t.^2 +  4 .* e .* t.^3;
acc =                 2 .* c +   6 .* d .* t + 12 .* e .* t.^2;

end