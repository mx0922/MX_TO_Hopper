function C = getRotMatFromBaseToWorld(q)

% 得到将base坐标系下的向量转换到世界坐标系下表示的旋转矩阵

x = q(1);   y = q(2);   z = q(3);

sx = sin(x);    cx = cos(x);
sy = sin(y);    cy = cos(y);
sz = sin(z);    cz = cos(z);

C = [ cy*cz, cz*sx*sy - cx*sz, sx*sz + cx*cz*sy;
      cy*sz, cx*cz + sx*sy*sz, cx*sy*sz - cz*sx;
        -sy,            cy*sx,           cx*cy];

end