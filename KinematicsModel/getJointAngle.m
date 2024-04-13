function q = getJointAngle(q0, pFootBase_ref, p)
% 本质就是利用数值逼近的方法求逆运动学

j1z = p.j1z;
l1 = p.l1;
l2 = p.l2;

% 得到实际的脚的位置
pFootBase_act = autoGen_PosFootBase(q0(1),q0(2),q0(3),q0(4),j1z,l1,l2);

% 期望与实际的误差
deltaP = pFootBase_ref - pFootBase_act;

% 容差设置
tol = 1e-9;

q_temp = q0;

if norm(deltaP) < tol
    q = q_temp;
end

while norm(deltaP) > tol
    
    Jac_foot_base = autoGen_JacFootBase(q_temp(1),q_temp(2),q_temp(3),q_temp(4),j1z,l1,l2);
    % KEY LINE !!!!
    q = q_temp + pinv(Jac_foot_base) * deltaP;
    
    pFootBase_act = autoGen_PosFootBase(q(1),q(2),q(3),q(4),j1z,l1,l2);
    deltaP = pFootBase_ref - pFootBase_act;
    
    q_temp = q;
end

end