function P_ankle_base = autoGen_PosFootBase(q1,q2,q3,q4,j1z,l1,l2)
%AUTOGEN_POSFOOTBASE
%    P_ANKLE_BASE = AUTOGEN_POSFOOTBASE(Q1,Q2,Q3,Q4,J1Z,L1,L2)

%    This function was generated by the Symbolic Math Toolbox version 8.4.
%    18-Oct-2021 14:23:07

t2 = cos(q1);
t3 = cos(q2);
t4 = cos(q3);
t5 = cos(q4);
t6 = sin(q1);
t7 = sin(q2);
t8 = sin(q3);
t9 = sin(q4);
t10 = l2.*t5;
t11 = l1+t10;
P_ankle_base = [-t11.*(t2.*t8+t4.*t6.*t7)-l2.*t9.*(t2.*t4-t6.*t7.*t8);-t11.*(t6.*t8-t2.*t4.*t7)-l2.*t9.*(t4.*t6+t2.*t7.*t8);-j1z-t3.*t4.*t11+l2.*t3.*t8.*t9];