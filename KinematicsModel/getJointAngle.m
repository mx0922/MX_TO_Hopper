function q = getJointAngle(q0, pFootBase_ref, p)
% ���ʾ���������ֵ�ƽ��ķ��������˶�ѧ

j1z = p.j1z;
l1 = p.l1;
l2 = p.l2;

% �õ�ʵ�ʵĽŵ�λ��
pFootBase_act = autoGen_PosFootBase(q0(1),q0(2),q0(3),q0(4),j1z,l1,l2);

% ������ʵ�ʵ����
deltaP = pFootBase_ref - pFootBase_act;

% �ݲ�����
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