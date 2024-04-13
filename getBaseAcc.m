function ddxBase = getBaseAcc(T, p, f, xLow, xUpp, param, t, t1, T1)

x0 = p(1);      dx0 = p(2);         x1 = p(3);      dx1 = p(4);
y0 = p(5);      dy0 = p(6);         y1 = p(7);      dy1 = p(8);
z0 = p(9);      dz0 = p(10);        z1 = p(11);     dz1 = p(12);

fx0 = f(1);      dfx0 = f(2);         fx1 = f(3);      dfx1 = f(4);
fy0 = f(5);      dfy0 = f(6);         fy1 = f(7);      dfy1 = f(8);
fz0 = f(9);      dfz0 = f(10);        fz1 = f(11);     dfz1 = f(12);

rx = xLow(1);      ry = xLow(2);      rz = xLow(3);      qx = xLow(4);      qy = xLow(5);      qz = xLow(6); 
drx = xLow(7);     dry = xLow(8);     drz = xLow(9);     dqx = xLow(10);    dqy = xLow(11);    dqz = xLow(12); 

rx1 = xUpp(1);      ry1 = xUpp(2);      rz1 = xUpp(3);      qx1 = xUpp(4);      qy1 = xUpp(5);      qz1 = xUpp(6); 
drx1 = xUpp(7);     dry1 = xUpp(8);     drz1 = xUpp(9);     dqx1 = xUpp(10);    dqy1 = xUpp(11);    dqz1 = xUpp(12);

m = param.m;
g = param.g;
Ix = param.I(1, 1);
Iy = param.I(2, 2);
Iz = param.I(3, 3);

ddxBase = autoGen_dynamicsCst_Maple(T,x0,dx0,x1,dx1,y0,dy0,y1,dy1,z0,dz0,z1,dz1,...
    fx0,dfx0,fx1,dfx1,fy0,dfy0,fy1,dfy1,fz0,dfz0,fz1,dfz1,...
    rx,ry,rz,qx,qy,qz,drx,dry,drz,dqx,dqy,dqz,...
    rx1,ry1,rz1,qx1,qy1,qz1,drx1,dry1,drz1,dqx1,dqy1,dqz1,...
    m,g,Ix,Iy,Iz,t,t1,T1);

end