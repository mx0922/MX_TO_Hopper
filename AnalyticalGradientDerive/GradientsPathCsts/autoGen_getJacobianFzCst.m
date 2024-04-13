function Jac_pz = autoGen_getJacobianFzCst(T,z0,dz0,z1,dz1,t)
%AUTOGEN_GETJACOBIANFZCST
%    JAC_PZ = AUTOGEN_GETJACOBIANFZCST(T,Z0,DZ0,Z1,DZ1,T)

%    This function was generated by the Symbolic Math Toolbox version 8.4.
%    18-Oct-2021 15:45:22

t2 = T.*dz0;
t3 = T.*dz1;
t4 = t.^2;
t5 = t.^3;
t6 = z0.*2.0;
t7 = z0.*3.0;
t8 = z1.*2.0;
t9 = z1.*3.0;
t11 = 1.0./T;
t10 = t2.*2.0;
t12 = t11.^2;
t13 = t11.^3;
t14 = -t8;
t15 = -t9;
t16 = t5.*t12;
t17 = t4.*t12.*3.0;
t18 = t5.*t13.*2.0;
t20 = t2+t3+t6+t14;
t21 = t3+t7+t10+t15;
t19 = -t16;
Jac_pz = [-t5.*t13.*(dz0+dz1)-t4.*t13.*t21.*2.0+t12.*t16.*t20.*3.0+t4.*t12.*(dz0.*2.0+dz1),t17-t18-1.0,-t+t19+t4.*t11.*2.0,-t17+t18,t19+t4.*t11,-dz0+t.*t12.*t21.*2.0-t4.*t13.*t20.*3.0];