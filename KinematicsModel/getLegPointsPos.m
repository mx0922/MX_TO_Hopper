function pPoints = getLegPointsPos(rBase,qBase,q,p)
% �õ��Ȳ���������������ϵ�µ���άλ��

j1z = p.j1z;
l1 = p.l1;
l2 = p.l2;

rx = rBase(1);      ry = rBase(2);      rz = rBase(3);
qx = qBase(1);      qy = qBase(2);      qz = qBase(3);

pPoints = autoGen_legPointsPos(rx,ry,rz,qx,qy,qz,q(1),q(2),q(3),q(4),j1z,l1,l2);

end