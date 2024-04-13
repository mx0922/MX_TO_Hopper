function Jac_kin = getJacobianKinCst(T, p, x, t)

global param

x0 = p(1);      dx0 = p(2);         x1 = p(3);      dx1 = p(4);
y0 = p(5);      dy0 = p(6);         y1 = p(7);      dy1 = p(8); 
z0 = p(9);      dz0 = p(10);        z1 = p(11);     dz1 = p(12);

rx = x(1);      ry = x(2);      rz = x(3); 
qx = x(4);      qy = x(5);      qz = x(6);

p_hat_x = param.p_hat(1);
p_hat_y = param.p_hat(2);
p_hat_z = param.p_hat(3);
L_cube = param.L_cube;

Jac_kin = autoGen_getJacobianKinCst(T,x0,dx0,x1,dx1,y0,dy0,y1,dy1,z0,dz0,z1,dz1,rx,ry,rz,qx,qy,qz,t,p_hat_x,p_hat_y,p_hat_z,L_cube);

end