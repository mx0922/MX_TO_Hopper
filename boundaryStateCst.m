function [c, ceq] = boundaryStateCst(x0, xF, p)

% 给定的期望初始和结束时的base位姿
q0 = p.q0;
qF = p.qF;

c = [];

% 初始和结束的位置和姿态的等式约束
ceq1 = q0 - x0(1:6);
ceq2 = qF - xF(1:6);

ceq = [ceq1; ceq2];

end