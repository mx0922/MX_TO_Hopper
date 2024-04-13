function [c, ceq] = getDynamicsCsts(x, dx, T)

Nt = size(x, 2);

r = x(1:3, :);
q = x(4:6, :);

dr = x(7:9, :);
dq = x(10:12, :);

ddr = dx(7:9, :);
ddq = dx(10:12, :);

dt = T / (Nt - 1);

% px,py,pz,qx,qy,qz依次进行
[c_px, ceq_px] = getBasePolyCsts(r, dr, ddr, 'x', dt);
[c_py, ceq_py] = getBasePolyCsts(r, dr, ddr, 'y', dt);
[c_pz, ceq_pz] = getBasePolyCsts(r, dr, ddr, 'z', dt);

[c_qx, ceq_qx] = getBasePolyCsts(q, dq, ddq, 'x', dt);
[c_qy, ceq_qy] = getBasePolyCsts(q, dq, ddq, 'y', dt);
[c_qz, ceq_qz] = getBasePolyCsts(q, dq, ddq, 'z', dt);

% 整合
c = [c_px; c_py; c_pz; c_qx; c_qy; c_qz];
ceq = [ceq_px; ceq_py; ceq_pz; ceq_qx; ceq_qy; ceq_qz];

end