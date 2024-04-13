function Jac_CWC = getJacobianCWCCst(T, f, miu, t)

x0 = f(1);      dx0 = f(2);     x1 = f(3);      dx1 = f(4);
y0 = f(5);      dy0 = f(6);     y1 = f(7);      dy1 = f(8);
z0 = f(9);      dz0 = f(10);    z1 = f(11);     dz1 = f(12);

Jac_CWC = autoGen_getJacobianCWCCst(T,x0,dx0,x1,dx1,y0,dy0,y1,dy1,z0,dz0,z1,dz1,miu,t);

end