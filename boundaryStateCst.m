function [c, ceq] = boundaryStateCst(x0, xF, p)

% ������������ʼ�ͽ���ʱ��baseλ��
q0 = p.q0;
qF = p.qF;

c = [];

% ��ʼ�ͽ�����λ�ú���̬�ĵ�ʽԼ��
ceq1 = q0 - x0(1:6);
ceq2 = qF - xF(1:6);

ceq = [ceq1; ceq2];

end