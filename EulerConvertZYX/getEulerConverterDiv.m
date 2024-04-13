function dC = getEulerConverterDiv(q, dq)

y = q(2);       z = q(3);
dy = dq(2);     dz = dq(3);

sy = sin(y);    cy = cos(y);
sz = sin(z);    cz = cos(z);

dC = [-dy*cz*sy - dz*cy*sz, -dz*cz, 0;...
       dz*cy*cz - dy*sy*sz, -dz*sz, 0;...
                    -dy*cy,      0, 0];

end