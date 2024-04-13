function C = getRotMatFromBaseToWorld(q)

% �õ���base����ϵ�µ�����ת������������ϵ�±�ʾ����ת����

x = q(1);   y = q(2);   z = q(3);

sx = sin(x);    cx = cos(x);
sy = sin(y);    cy = cos(y);
sz = sin(z);    cz = cos(z);

C = [ cy*cz, cz*sx*sy - cx*sz, sx*sz + cx*cz*sy;
      cy*sz, cx*cz + sx*sy*sz, cx*sy*sz - cz*sx;
        -sy,            cy*sx,           cx*cy];

end