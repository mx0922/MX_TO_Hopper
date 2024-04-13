function C = getEulerConverter(q)

y = q(2);   z = q(3);

sy = sin(y);    cy = cos(y);
sz = sin(z);    cz = cos(z);

C = [cy*cz, -sz, 0;...
     cy*sz,  cz, 0;...
       -sy,   0, 1];

end