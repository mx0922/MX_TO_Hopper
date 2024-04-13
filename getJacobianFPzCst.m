function Jac_pz = getJacobianFPzCst(T, fp, t)

z0 = fp(1);     dz0 = fp(2);    z1 = fp(3);     dz1 = fp(4);

Jac_pz = autoGen_getJacobianFzCst(T,z0,dz0,z1,dz1,t);

end