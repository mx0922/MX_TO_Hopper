function [cGrad, ceqGrad] = getPathCstGrad(z, pack)

global param

[tPhase, p, f, x] = unPackDecVar(z, pack);

T1 = tPhase(1);      T2 = tPhase(2);      T3 = tPhase(3);      T4 = tPhase(4);      T5 = tPhase(5);

T = sum(tPhase);
Nt = size(x, 2);

dt = T / (Nt - 1);

nDecVar = numel(z);

%% 脚的运动学约束
c1Grad = zeros(Nt, nDecVar);
c3Grad = zeros(Nt, nDecVar);

tol = 1e-6;

for ii = 1:Nt
    t = (ii - 1) * dt;
    
    if t <= T1 + tol
        tt = t;
        
        % kinematic csts
        idxP = [1:4, 28+(1:4), 56+(1:4)];
        Jac_kin = getJacobianKinCst(T1, p(idxP), x(1:6, ii), tt);
        % c1Grad(ii, [1, 75, 76, 77, 78, 79, 80]) = Jac_mx([1, 14, 15, 16, 17, 18, 19]);
        idxDecVarKin = [1, 12*ii+(63:68)];
        idxJacKin = [1, 13+(1:6)];
        c1Grad(ii, idxDecVarKin) = Jac_kin(idxJacKin);
        
        dtGrad = (ii - 1)/(Nt - 1) * ones(1, 5);
        c1Grad(ii, 1:5) = c1Grad(ii, 1:5) + Jac_kin(end) * dtGrad;
        
        % pz >= 0  <==> -pz <= 0
        idxPz = 56 + (1:4);
        Jac_pz = getJacobianFPzCst(T1,p(idxPz),tt);
        
        idxDecVarPz = 1;
        idxJacPz = 1;
        c3Grad(ii, idxDecVarPz) = Jac_pz(idxJacPz);
        
        c3Grad(ii, 1:5) = c3Grad(ii, 1:5) + Jac_pz(end) * dtGrad;
        
    elseif t > T1 + tol && t <= T1 + T2 + tol
        tt = t - T1;
        if tt <= T2/3 + tol
            
            % kinematic csts
            idxP = 4 + [1:4, 28+(1:4), 56+(1:4)];
            Jac_kin = getJacobianKinCst(T2/3, p(idxP), x(1:6, ii), tt);
            
            Jac_kin(1) = Jac_kin(1) / 3;

            idxDecVarKin = [2, 5+(1:2), 15+(1:2), 25+(1:2), 12*ii+(63:68)];
            idxJacKin = [1, 1+(3:4), 5+(3:4), 9+(3:4), 13+(1:6)];
            c1Grad(ii, idxDecVarKin) = Jac_kin(idxJacKin);
            
            dtGrad = (ii - 1)/(Nt - 1) * ones(1, 5) - [1, 0, 0, 0, 0];
            c1Grad(ii, 1:5) = c1Grad(ii, 1:5) + Jac_kin(end) * dtGrad;

            % pz >= 0  <==> -pz <= 0
            idxPz = 56 + 4 + (1:4);
            Jac_pz = getJacobianFPzCst(T2/3,p(idxPz),tt);
            
            Jac_pz(1) = Jac_pz(1) / 3;

            idxDecVarPz = [2, 25+(1:2)];
            idxJacPz = [1, 1+(3:4)];
            c3Grad(ii, idxDecVarPz) = Jac_pz(idxJacPz);
            
            c3Grad(ii, 1:5) = c3Grad(ii, 1:5) + Jac_pz(end) * dtGrad;
            
        elseif tt > T2/3 + tol && tt <= T2 * 2 / 3 + tol
                        
            t_temp = tt - T2/3;            
            % kinematic csts
            idxP = 4 + 2 + [1:4, 28+(1:4), 56+(1:4)];
            Jac_kin = getJacobianKinCst(T2/3, p(idxP), x(1:6, ii), t_temp);
            
            Jac_kin(1) = Jac_kin(1) / 3;

            idxDecVarKin = [2, 5+(1:4), 15+(1:4), 25+(1:4), 12*ii+(63:68)];
            idxJacKin = [1, 1+(1:4), 5+(1:4), 9+(1:4), 13+(1:6)];
            c1Grad(ii, idxDecVarKin) = Jac_kin(idxJacKin);
            
            dtGrad = (ii - 1)/(Nt - 1) * ones(1, 5) - [1, 1/3, 0, 0, 0];
            c1Grad(ii, 1:5) = c1Grad(ii, 1:5) + Jac_kin(end) * dtGrad;

            % pz >= 0  <==> -pz <= 0
            idxPz = 56 + 4 + 2 + (1:4);
            Jac_pz = getJacobianFPzCst(T2/3,p(idxPz),t_temp);
            
            Jac_pz(1) = Jac_pz(1) / 3;

            idxDecVarPz = [2, 25+(1:4)];
            idxJacPz = [1, 1+(1:4)];
            c3Grad(ii, idxDecVarPz) = Jac_pz(idxJacPz);
            
            c3Grad(ii, 1:5) = c3Grad(ii, 1:5) + Jac_pz(end) * dtGrad;
            
        else
            
            t_temp = tt - 2 * T2/3;            
            % kinematic csts
            idxP = 4 + 2 + 2 + [1:4, 28+(1:4), 56+(1:4)];
            Jac_kin = getJacobianKinCst(T2/3, p(idxP), x(1:6, ii), t_temp);
            
            Jac_kin(1) = Jac_kin(1) / 3;

            idxDecVarKin = [2, 5+2+(1:3), 15+2+(1:3), 25+2+(1:3), 12*ii+(63:68)];
            idxJacKin = [1, 1+(1:3), 5+(1:3), 9+(1:3), 13+(1:6)];
            c1Grad(ii, idxDecVarKin) = Jac_kin(idxJacKin);
            
            dtGrad = (ii - 1)/(Nt - 1) * ones(1, 5) - [1, 2/3, 0, 0, 0];
            c1Grad(ii, 1:5) = c1Grad(ii, 1:5) + Jac_kin(end) * dtGrad;

            % pz >= 0  <==> -pz <= 0
            idxPz = 56 + 4 + 2 + 2 + (1:4);
            Jac_pz = getJacobianFPzCst(T2/3,p(idxPz),t_temp);
            
            Jac_pz(1) = Jac_pz(1) / 3;

            idxDecVarPz = [2, 25+2+(1:3)];
            idxJacPz = [1, 1+(1:3)];
            c3Grad(ii, idxDecVarPz) = Jac_pz(idxJacPz);
            
            c3Grad(ii, 1:5) = c3Grad(ii, 1:5) + Jac_pz(end) * dtGrad;
        end
    elseif t > T1 + T2 + tol && t <= T1 + T2 + T3 + tol
        tt = t - T1 - T2;
        
        % kinematic csts
        idxP = 12 + [1:4, 28+(1:4), 56+(1:4)];
        Jac_kin = getJacobianKinCst(T3, p(idxP), x(1:6, ii), tt);
        idxDecVarKin = [3, 10, 20, 30, 12*ii+(63:68)];
        
        idxJacKin_nom = [1, 2, 4, 6, 8, 10, 12, 13+(1:6)];
        
        % 对Jac_kin中的元素进行操作
        Jac_kin(idxJacKin_nom(2)) = Jac_kin(idxJacKin_nom(2)) + Jac_kin(idxJacKin_nom(3));
        Jac_kin(idxJacKin_nom(4)) = Jac_kin(idxJacKin_nom(4)) + Jac_kin(idxJacKin_nom(5));
        Jac_kin(idxJacKin_nom(6)) = Jac_kin(idxJacKin_nom(6)) + Jac_kin(idxJacKin_nom(7));
        
        idxJacKin = [1, 2, 6, 10, 13+(1:6)];
        
        c1Grad(ii, idxDecVarKin) = Jac_kin(idxJacKin);     
                
        % pz >= 0  <==> -pz <= 0
        idxPz = 56 + 12 + (1:4);
        Jac_pz = getJacobianFPzCst(T3,p(idxPz),tt);
        
        idxDecVarPz = [3, 30];
        
        idxJacPz_nom = [1, 2, 4];
        % 对Jac_pz中的元素进行操作
        Jac_pz(idxJacPz_nom(2)) = Jac_pz(idxJacPz_nom(2)) + Jac_pz(idxJacPz_nom(3));
        
        idxJacPz = [1, 2];
        
        c3Grad(ii, idxDecVarPz) = Jac_pz(idxJacPz);
        
        dtGrad = (ii - 1)/(Nt - 1) * ones(1, 5) - [1, 1, 0, 0, 0];
        c1Grad(ii, 1:5) = c1Grad(ii, 1:5) + Jac_kin(end) * dtGrad;
        c3Grad(ii, 1:5) = c3Grad(ii, 1:5) + Jac_pz(end) * dtGrad;
        
    elseif t > T1 + T2 + T3 + tol && t <= T1 + T2 + T3 + T4 + tol
        tt = t - T1 - T2 - T3;
        if tt <= T4/3 + tol
                       
            % kinematic csts
            idxP = 4 + 12 + [1:4, 28+(1:4), 56+(1:4)];
            Jac_kin = getJacobianKinCst(T4/3, p(idxP), x(1:6, ii), tt);
            
            Jac_kin(1) = Jac_kin(1) / 3;

            idxDecVarKin = [4, 10, 5+5+(1:2), 20, 5+15+(1:2), 30, 5+25+(1:2), 12*ii+(63:68)];
            idxJacKin = [1, 1+1, 1+(3:4), 1+5, 5+(3:4), 1+9, 9+(3:4), 13+(1:6)];
            c1Grad(ii, idxDecVarKin) = Jac_kin(idxJacKin);

            % pz >= 0  <==> -pz <= 0
            idxPz = 56 + 4 + 12 + (1:4);
            Jac_pz = getJacobianFPzCst(T4/3,p(idxPz),tt);
            
            Jac_pz(1) = Jac_pz(1) / 3;

            idxDecVarPz = [4, 25+5, 25+5+(1:2)];
            idxJacPz = [1, 1+1, 1+(3:4)];
            c3Grad(ii, idxDecVarPz) = Jac_pz(idxJacPz);
            
            dtGrad = (ii - 1)/(Nt - 1) * ones(1, 5) - [1, 1, 1, 0, 0];
            c1Grad(ii, 1:5) = c1Grad(ii, 1:5) + Jac_kin(end) * dtGrad;
            c3Grad(ii, 1:5) = c3Grad(ii, 1:5) + Jac_pz(end) * dtGrad;
            
        elseif tt > T4/3 + tol && tt <= T4 * 2 / 3 + tol
                        
            t_temp = tt - T4/3;            
            % kinematic csts
            idxP = 4 + 12 + 2 + [1:4, 28+(1:4), 56+(1:4)];
            Jac_kin = getJacobianKinCst(T4/3, p(idxP), x(1:6, ii), t_temp);
            
            Jac_kin(1) = Jac_kin(1) / 3;

            idxDecVarKin = [4, 5+5+(1:4), 15+5+(1:4), 25+5+(1:4), 12*ii+(63:68)];
            idxJacKin = [1, 1+(1:4), 5+(1:4), 9+(1:4), 13+(1:6)];
            c1Grad(ii, idxDecVarKin) = Jac_kin(idxJacKin);

            % pz >= 0  <==> -pz <= 0
            idxPz = 56 + 4 + 12 + 2 + (1:4);
            Jac_pz = getJacobianFPzCst(T4/3,p(idxPz),t_temp);
            
            Jac_pz(1) = Jac_pz(1) / 3;

            idxDecVarPz = [2, 25+5+(1:4)];
            idxJacPz = [1, 1+(1:4)];
            c3Grad(ii, idxDecVarPz) = Jac_pz(idxJacPz);
            
            dtGrad = (ii - 1)/(Nt - 1) * ones(1, 5) - [1, 1, 1, 1/3, 0];
            c1Grad(ii, 1:5) = c1Grad(ii, 1:5) + Jac_kin(end) * dtGrad;
            c3Grad(ii, 1:5) = c3Grad(ii, 1:5) + Jac_pz(end) * dtGrad;
            
        else
                       
            t_temp = tt - 2 * T4/3;            
            % kinematic csts
            idxP = 4 + 12 + 2 + 2 + [1:4, 28+(1:4), 56+(1:4)];
            Jac_kin = getJacobianKinCst(T4/3, p(idxP), x(1:6, ii), t_temp);
            
            Jac_kin(1) = Jac_kin(1) / 3;

            idxDecVarKin = [4, 5+5+2+(1:3), 15+5+2+(1:3), 25+5+2+(1:3), 12*ii+(63:68)];
            idxJacKin = [1, 1+(1:3), 5+(1:3), 9+(1:3), 13+(1:6)];
            c1Grad(ii, idxDecVarKin) = Jac_kin(idxJacKin);

            % pz >= 0  <==> -pz <= 0
            idxPz = 56 + 12 + 4 + 2 + 2 + (1:4);
            Jac_pz = getJacobianFPzCst(T4/3,p(idxPz),t_temp);
            
            Jac_pz(1) = Jac_pz(1) / 3;

            idxDecVarPz = [4, 25+5+2+(1:3)];
            idxJacPz = [1, 1+(1:3)];
            c3Grad(ii, idxDecVarPz) = Jac_pz(idxJacPz);
            
            dtGrad = (ii - 1)/(Nt - 1) * ones(1, 5) - [1, 1, 1, 2/3, 0];
            c1Grad(ii, 1:5) = c1Grad(ii, 1:5) + Jac_kin(end) * dtGrad;
            c3Grad(ii, 1:5) = c3Grad(ii, 1:5) + Jac_pz(end) * dtGrad;
            
        end
        
    else
        tt = t - T1 - T2 - T3 - T4;
                
        % kinematic csts
        idxP = 12 + 12 + [1:4, 28+(1:4), 56+(1:4)];
        Jac_kin = getJacobianKinCst(T5, p(idxP), x(1:6, ii), tt);
        idxDecVarKin = [5, 10+5, 20+5, 30+5, 12*ii+(63:68)];
        
        idxJacKin_nom = [1, 2, 4, 6, 8, 10, 12, 13+(1:6)];
        
        % 对Jac_kin中的元素进行操作
        Jac_kin(idxJacKin_nom(2)) = Jac_kin(idxJacKin_nom(2)) + Jac_kin(idxJacKin_nom(3));
        Jac_kin(idxJacKin_nom(4)) = Jac_kin(idxJacKin_nom(4)) + Jac_kin(idxJacKin_nom(5));
        Jac_kin(idxJacKin_nom(6)) = Jac_kin(idxJacKin_nom(6)) + Jac_kin(idxJacKin_nom(7));
        
        idxJacKin = [1, 2, 6, 10, 13+(1:6)];
        
        c1Grad(ii, idxDecVarKin) = Jac_kin(idxJacKin);
        
        % pz >= 0  <==> -pz <= 0
        idxPz = 56 + 12 + 12 + (1:4);
        Jac_pz = getJacobianFPzCst(T5,p(idxPz),tt);
        
        idxDecVarPz = [5, 30+5];
        
        idxJacPz_nom = [1, 2, 4];
        % 对Jac_pz中的元素进行操作
        Jac_pz(idxJacPz_nom(2)) = Jac_pz(idxJacPz_nom(2)) + Jac_pz(idxJacPz_nom(3));
        
        idxJacPz = [1, 2];
        
        c3Grad(ii, idxDecVarPz) = Jac_pz(idxJacPz);
        
        dtGrad = (ii - 1)/(Nt - 1) * ones(1, 5) - [1, 1, 1, 1, 0];
        c1Grad(ii, 1:5) = c1Grad(ii, 1:5) + Jac_kin(end) * dtGrad;
        c3Grad(ii, 1:5) = c3Grad(ii, 1:5) + Jac_pz(end) * dtGrad;
        
    end
    
end

%% 接触力约束
% -fz <= 0
c2Grad = zeros(Nt, nDecVar);

miu = param.miu;

cwcGrad1 = zeros(Nt, nDecVar);
cwcGrad2 = zeros(Nt, nDecVar);
cwcGrad3 = zeros(Nt, nDecVar);
cwcGrad4 = zeros(Nt, nDecVar);

for ii = 1:Nt
    t = (ii - 1) * dt;
    
    if t <= T1 + tol        
        tt = t;
        if tt <= T1/3 + tol
            
            % -fz <= 0
            idxFz = 64 + (1:4);            
            Jac_fz = getJacobianFPzCst(T1/3,f(idxFz),tt);
            
            Jac_fz(1) = Jac_fz(1)/3;
            
            idxDecVarFz = [1, 35+26+(1:2)];
            idxJacFz = [1, 1+(3:4)];
            c2Grad(ii, idxDecVarFz) = Jac_fz(idxJacFz);  
            
            % CWC
            idxCWC = [1:4, 32+(1:4), 64+(1:4)];            
            Jac_CWC = getJacobianCWCCst(T1/3, f(idxCWC), miu, tt);
            
            Jac_CWC(:, 1) = Jac_CWC(:, 1)/3;
            
            idxDecVarCWCx = [1,    35+(1:2), 35+26+(1:2)];
            idxDecVarCWCy = [1, 35+13+(1:2), 35+26+(1:2)];
            
            idxJacFx = [1,   1+(3:4), 1+8+(3:4)];
            idxJacFy = [1, 1+4+(3:4), 1+8+(3:4)];
            
            cwcGrad1(ii, idxDecVarCWCx) = Jac_CWC(1, idxJacFx);
            cwcGrad2(ii, idxDecVarCWCx) = Jac_CWC(2, idxJacFx);
            cwcGrad3(ii, idxDecVarCWCy) = Jac_CWC(3, idxJacFy);
            cwcGrad4(ii, idxDecVarCWCy) = Jac_CWC(4, idxJacFy);
            
            dtGrad = (ii - 1)/(Nt - 1) * ones(1, 5) - [0, 0, 0, 0, 0];
            c2Grad(ii, 1:5) = c2Grad(ii, 1:5) + Jac_fz(end) * dtGrad;
            cwcGrad1(ii, 1:5) = cwcGrad1(ii, 1:5) + Jac_CWC(1, end) * dtGrad;
            cwcGrad2(ii, 1:5) = cwcGrad2(ii, 1:5) + Jac_CWC(2, end) * dtGrad;
            cwcGrad3(ii, 1:5) = cwcGrad3(ii, 1:5) + Jac_CWC(3, end) * dtGrad;
            cwcGrad4(ii, 1:5) = cwcGrad4(ii, 1:5) + Jac_CWC(4, end) * dtGrad;
            
        elseif tt > T1/3 + tol && tt <= T1 * 2 / 3 + tol
            
            t_temp = tt - T1/3;
            % -fz <= 0
            idxFz = 64 + 2 + (1:4);            
            Jac_fz = getJacobianFPzCst(T1/3,f(idxFz),t_temp);
            
            Jac_fz(1) = Jac_fz(1)/3;
            
            idxDecVarFz = [1, 35+26+(1:4)];
            idxJacFz = [1, 1+(1:4)];
            c2Grad(ii, idxDecVarFz) = Jac_fz(idxJacFz);  
            
            % CWC
            idxCWC = 2+[1:4, 32+(1:4), 64+(1:4)];            
            Jac_CWC = getJacobianCWCCst(T1/3, f(idxCWC), miu, t_temp);
            
            Jac_CWC(:, 1) = Jac_CWC(:, 1)/3;
            
            idxDecVarCWCx = [1,    35+(1:4), 35+26+(1:4)];
            idxDecVarCWCy = [1, 35+13+(1:4), 35+26+(1:4)];
            
            idxJacFx = [1,   1+(1:4), 1+8+(1:4)];
            idxJacFy = [1, 1+4+(1:4), 1+8+(1:4)];
            
            cwcGrad1(ii, idxDecVarCWCx) = Jac_CWC(1, idxJacFx);
            cwcGrad2(ii, idxDecVarCWCx) = Jac_CWC(2, idxJacFx);
            cwcGrad3(ii, idxDecVarCWCy) = Jac_CWC(3, idxJacFy);
            cwcGrad4(ii, idxDecVarCWCy) = Jac_CWC(4, idxJacFy);
            
            dtGrad = (ii - 1)/(Nt - 1) * ones(1, 5) - [1/3, 0, 0, 0, 0];
            c2Grad(ii, 1:5) = c2Grad(ii, 1:5) + Jac_fz(end) * dtGrad;
            cwcGrad1(ii, 1:5) = cwcGrad1(ii, 1:5) + Jac_CWC(1, end) * dtGrad;
            cwcGrad2(ii, 1:5) = cwcGrad2(ii, 1:5) + Jac_CWC(2, end) * dtGrad;
            cwcGrad3(ii, 1:5) = cwcGrad3(ii, 1:5) + Jac_CWC(3, end) * dtGrad;
            cwcGrad4(ii, 1:5) = cwcGrad4(ii, 1:5) + Jac_CWC(4, end) * dtGrad;
            
        else
                        
            t_temp = tt - 2 * T1/3;
            % -fz <= 0
            idxFz = 64 + 2 + 2 + (1:4);            
            Jac_fz = getJacobianFPzCst(T1/3,f(idxFz),t_temp);
            
            Jac_fz(1) = Jac_fz(1)/3;
            
            idxDecVarFz = [1, 35+26+(3:4)];
            idxJacFz = [1, 1+(1:2)];
            c2Grad(ii, idxDecVarFz) = Jac_fz(idxJacFz);  
            
            % CWC
            idxCWC = 2+2+[1:4, 32+(1:4), 64+(1:4)];            
            Jac_CWC = getJacobianCWCCst(T1/3, f(idxCWC), miu, t_temp);
            
            Jac_CWC(:, 1) = Jac_CWC(:, 1)/3;
            
            idxDecVarCWCx = [1,    35+(3:4), 35+26+(3:4)];
            idxDecVarCWCy = [1, 35+13+(3:4), 35+26+(3:4)];
            
            idxJacFx = [1,   1+(1:2), 1+8+(1:2)];
            idxJacFy = [1, 1+4+(1:2), 1+8+(1:2)];
            
            cwcGrad1(ii, idxDecVarCWCx) = Jac_CWC(1, idxJacFx);
            cwcGrad2(ii, idxDecVarCWCx) = Jac_CWC(2, idxJacFx);
            cwcGrad3(ii, idxDecVarCWCy) = Jac_CWC(3, idxJacFy);
            cwcGrad4(ii, idxDecVarCWCy) = Jac_CWC(4, idxJacFy);
            
            dtGrad = (ii - 1)/(Nt - 1) * ones(1, 5) - [2/3, 0, 0, 0, 0];
            c2Grad(ii, 1:5) = c2Grad(ii, 1:5) + Jac_fz(end) * dtGrad;
            cwcGrad1(ii, 1:5) = cwcGrad1(ii, 1:5) + Jac_CWC(1, end) * dtGrad;
            cwcGrad2(ii, 1:5) = cwcGrad2(ii, 1:5) + Jac_CWC(2, end) * dtGrad;
            cwcGrad3(ii, 1:5) = cwcGrad3(ii, 1:5) + Jac_CWC(3, end) * dtGrad;
            cwcGrad4(ii, 1:5) = cwcGrad4(ii, 1:5) + Jac_CWC(4, end) * dtGrad;
            
        end
        
    elseif t > T1 + tol && t <= T1 + T2 + tol
        tt = t - T1;
                
        % -fz <= 0
        idxFz = 64 + 8 + (1:4);            
        Jac_fz = getJacobianFPzCst(T2,f(idxFz),tt);

        idxDecVarFz = 2;
        idxJacFz = 1;
        c2Grad(ii, idxDecVarFz) = Jac_fz(idxJacFz);  

        % CWC
        idxCWC = 8 + [1:4, 32+(1:4), 64+(1:4)];            
        Jac_CWC = getJacobianCWCCst(T2, f(idxCWC), miu, tt);

        idxDecVarCWCx = 2;
        idxDecVarCWCy = 2;

        idxJacFx = 1;
        idxJacFy = 1;

        cwcGrad1(ii, idxDecVarCWCx) = Jac_CWC(1, idxJacFx);
        cwcGrad2(ii, idxDecVarCWCx) = Jac_CWC(2, idxJacFx);
        cwcGrad3(ii, idxDecVarCWCy) = Jac_CWC(3, idxJacFy);
        cwcGrad4(ii, idxDecVarCWCy) = Jac_CWC(4, idxJacFy);
        
        dtGrad = (ii - 1)/(Nt - 1) * ones(1, 5) - [1, 0, 0, 0, 0];
        c2Grad(ii, 1:5) = c2Grad(ii, 1:5) + Jac_fz(end) * dtGrad;
        cwcGrad1(ii, 1:5) = cwcGrad1(ii, 1:5) + Jac_CWC(1, end) * dtGrad;
        cwcGrad2(ii, 1:5) = cwcGrad2(ii, 1:5) + Jac_CWC(2, end) * dtGrad;
        cwcGrad3(ii, 1:5) = cwcGrad3(ii, 1:5) + Jac_CWC(3, end) * dtGrad;
        cwcGrad4(ii, 1:5) = cwcGrad4(ii, 1:5) + Jac_CWC(4, end) * dtGrad;
        
    elseif t > T1 + T2 + tol && t <= T1 + T2 + T3 + tol
        tt = t - T1 - T2;
        if tt <= T3/3 + tol
                        
            % -fz <= 0
            idxFz = 64 + 12 + (1:4);            
            Jac_fz = getJacobianFPzCst(T3/3,f(idxFz),tt);
            
            Jac_fz(1) = Jac_fz(1)/3;
            
            idxDecVarFz = [3, 35+26+4+(1:2)];
            idxJacFz = [1, 1+(3:4)];
            c2Grad(ii, idxDecVarFz) = Jac_fz(idxJacFz);  
            
            % CWC
            idxCWC = 12 + [1:4, 32+(1:4), 64+(1:4)];            
            Jac_CWC = getJacobianCWCCst(T3/3, f(idxCWC), miu, tt);
            
            Jac_CWC(:, 1) = Jac_CWC(:, 1)/3;
            
            idxDecVarCWCx = [3,    35+4+(1:2), 35+4+26+(1:2)];
            idxDecVarCWCy = [3, 35+4+13+(1:2), 35+4+26+(1:2)];
            
            idxJacFx = [1,   1+(3:4), 1+8+(3:4)];
            idxJacFy = [1, 1+4+(3:4), 1+8+(3:4)];
            
            cwcGrad1(ii, idxDecVarCWCx) = Jac_CWC(1, idxJacFx);
            cwcGrad2(ii, idxDecVarCWCx) = Jac_CWC(2, idxJacFx);
            cwcGrad3(ii, idxDecVarCWCy) = Jac_CWC(3, idxJacFy);
            cwcGrad4(ii, idxDecVarCWCy) = Jac_CWC(4, idxJacFy);
            
            dtGrad = (ii - 1)/(Nt - 1) * ones(1, 5) - [1, 1, 0, 0, 0];
            c2Grad(ii, 1:5) = c2Grad(ii, 1:5) + Jac_fz(end) * dtGrad;
            cwcGrad1(ii, 1:5) = cwcGrad1(ii, 1:5) + Jac_CWC(1, end) * dtGrad;
            cwcGrad2(ii, 1:5) = cwcGrad2(ii, 1:5) + Jac_CWC(2, end) * dtGrad;
            cwcGrad3(ii, 1:5) = cwcGrad3(ii, 1:5) + Jac_CWC(3, end) * dtGrad;
            cwcGrad4(ii, 1:5) = cwcGrad4(ii, 1:5) + Jac_CWC(4, end) * dtGrad;
            
        elseif tt > T3/3 + tol && tt <= T3 * 2 / 3 + tol
            
            t_temp = tt - T3/3;
            % -fz <= 0
            idxFz = 64 + 12 + 2 + (1:4);            
            Jac_fz = getJacobianFPzCst(T3/3,f(idxFz),t_temp);
            
            Jac_fz(1) = Jac_fz(1)/3;
            
            idxDecVarFz = [3, 35+26+4+(1:4)];
            idxJacFz = [1, 1+(1:4)];
            c2Grad(ii, idxDecVarFz) = Jac_fz(idxJacFz);  
            
            % CWC
            idxCWC = 12 + 2 + [1:4, 32+(1:4), 64+(1:4)];            
            Jac_CWC = getJacobianCWCCst(T3/3, f(idxCWC), miu, t_temp);
            
            Jac_CWC(:, 1) = Jac_CWC(:, 1)/3;
            
            idxDecVarCWCx = [3,    35+4+(1:4), 35+4+26+(1:4)];
            idxDecVarCWCy = [3, 35+4+13+(1:4), 35+4+26+(1:4)];
            
            idxJacFx = [1,   1+(1:4), 1+8+(1:4)];
            idxJacFy = [1, 1+4+(1:4), 1+8+(1:4)];
            
            cwcGrad1(ii, idxDecVarCWCx) = Jac_CWC(1, idxJacFx);
            cwcGrad2(ii, idxDecVarCWCx) = Jac_CWC(2, idxJacFx);
            cwcGrad3(ii, idxDecVarCWCy) = Jac_CWC(3, idxJacFy);
            cwcGrad4(ii, idxDecVarCWCy) = Jac_CWC(4, idxJacFy);
            
            dtGrad = (ii - 1)/(Nt - 1) * ones(1, 5) - [1, 1, 1/3, 0, 0];
            c2Grad(ii, 1:5) = c2Grad(ii, 1:5) + Jac_fz(end) * dtGrad;
            cwcGrad1(ii, 1:5) = cwcGrad1(ii, 1:5) + Jac_CWC(1, end) * dtGrad;
            cwcGrad2(ii, 1:5) = cwcGrad2(ii, 1:5) + Jac_CWC(2, end) * dtGrad;
            cwcGrad3(ii, 1:5) = cwcGrad3(ii, 1:5) + Jac_CWC(3, end) * dtGrad;
            cwcGrad4(ii, 1:5) = cwcGrad4(ii, 1:5) + Jac_CWC(4, end) * dtGrad;
            
        else
            
            t_temp = tt - 2 * T3/3;
            % -fz <= 0
            idxFz = 64 + 12 + 2 + 2 + (1:4);            
            Jac_fz = getJacobianFPzCst(T3/3,f(idxFz),t_temp);
            
            Jac_fz(1) = Jac_fz(1)/3;
            
            idxDecVarFz = [3, 35+26+4+(3:4)];
            idxJacFz = [1, 1+(1:2)];
            c2Grad(ii, idxDecVarFz) = Jac_fz(idxJacFz);  
            
            % CWC
            idxCWC = 12 + 2 + 2 + [1:4, 32+(1:4), 64+(1:4)];            
            Jac_CWC = getJacobianCWCCst(T3/3, f(idxCWC), miu, t_temp);
            
            Jac_CWC(:, 1) = Jac_CWC(:, 1)/3;
            
            idxDecVarCWCx = [3,    35+4+(3:4), 35+4+26+(3:4)];
            idxDecVarCWCy = [3, 35+4+13+(3:4), 35+4+26+(3:4)];
            
            idxJacFx = [1,   1+(1:2), 1+8+(1:2)];
            idxJacFy = [1, 1+4+(1:2), 1+8+(1:2)];
            
            cwcGrad1(ii, idxDecVarCWCx) = Jac_CWC(1, idxJacFx);
            cwcGrad2(ii, idxDecVarCWCx) = Jac_CWC(2, idxJacFx);
            cwcGrad3(ii, idxDecVarCWCy) = Jac_CWC(3, idxJacFy);
            cwcGrad4(ii, idxDecVarCWCy) = Jac_CWC(4, idxJacFy);
            
            dtGrad = (ii - 1)/(Nt - 1) * ones(1, 5) - [1, 1, 2/3, 0, 0];
            c2Grad(ii, 1:5) = c2Grad(ii, 1:5) + Jac_fz(end) * dtGrad;
            cwcGrad1(ii, 1:5) = cwcGrad1(ii, 1:5) + Jac_CWC(1, end) * dtGrad;
            cwcGrad2(ii, 1:5) = cwcGrad2(ii, 1:5) + Jac_CWC(2, end) * dtGrad;
            cwcGrad3(ii, 1:5) = cwcGrad3(ii, 1:5) + Jac_CWC(3, end) * dtGrad;
            cwcGrad4(ii, 1:5) = cwcGrad4(ii, 1:5) + Jac_CWC(4, end) * dtGrad;
            
        end
    elseif t > T1 + T2 + T3 + tol && t <= T1 + T2 + T3 + T4 + tol
        tt = t - T1 - T2 - T3;
                
        % -fz <= 0
        idxFz = 64 + 8 + 12 + (1:4);            
        Jac_fz = getJacobianFPzCst(T4,f(idxFz),tt);

        idxDecVarFz = 4;
        idxJacFz = 1;
        c2Grad(ii, idxDecVarFz) = Jac_fz(idxJacFz);  

        % CWC
        idxCWC = 8 + 12 + [1:4, 32+(1:4), 64+(1:4)];            
        Jac_CWC = getJacobianCWCCst(T4, f(idxCWC), miu, tt);

        idxDecVarCWCx = 4;
        idxDecVarCWCy = 4;

        idxJacFx = 1;
        idxJacFy = 1;

        cwcGrad1(ii, idxDecVarCWCx) = Jac_CWC(1, idxJacFx);
        cwcGrad2(ii, idxDecVarCWCx) = Jac_CWC(2, idxJacFx);
        cwcGrad3(ii, idxDecVarCWCy) = Jac_CWC(3, idxJacFy);
        cwcGrad4(ii, idxDecVarCWCy) = Jac_CWC(4, idxJacFy);
        
        dtGrad = (ii - 1)/(Nt - 1) * ones(1, 5) - [1, 1, 1, 0, 0];
        c2Grad(ii, 1:5) = c2Grad(ii, 1:5) + Jac_fz(end) * dtGrad;
        cwcGrad1(ii, 1:5) = cwcGrad1(ii, 1:5) + Jac_CWC(1, end) * dtGrad;
        cwcGrad2(ii, 1:5) = cwcGrad2(ii, 1:5) + Jac_CWC(2, end) * dtGrad;
        cwcGrad3(ii, 1:5) = cwcGrad3(ii, 1:5) + Jac_CWC(3, end) * dtGrad;
        cwcGrad4(ii, 1:5) = cwcGrad4(ii, 1:5) + Jac_CWC(4, end) * dtGrad;
        
    else
        tt = t - T1 - T2 - T3 - T4;
        if tt <= T5/3 + tol
                        
            % -fz <= 0
            idxFz = 64 + 12 + 12 + (1:4);            
            Jac_fz = getJacobianFPzCst(T5/3,f(idxFz),tt);
            
            Jac_fz(1) = Jac_fz(1)/3;
            
            idxDecVarFz = [5, 35+26+4+4+(1:2)];
            idxJacFz = [1, 1+(3:4)];
            c2Grad(ii, idxDecVarFz) = Jac_fz(idxJacFz);  
            
            % CWC
            idxCWC = 12 + 12 + [1:4, 32+(1:4), 64+(1:4)];            
            Jac_CWC = getJacobianCWCCst(T5/3, f(idxCWC), miu, tt);
            
            Jac_CWC(:, 1) = Jac_CWC(:, 1)/3;
            
            idxDecVarCWCx = [5,    35+4+4+(1:2), 35+4+4+26+(1:2)];
            idxDecVarCWCy = [5, 35+4+4+13+(1:2), 35+4+4+26+(1:2)];
            
            idxJacFx = [1,   1+(3:4), 1+8+(3:4)];
            idxJacFy = [1, 1+4+(3:4), 1+8+(3:4)];
            
            cwcGrad1(ii, idxDecVarCWCx) = Jac_CWC(1, idxJacFx);
            cwcGrad2(ii, idxDecVarCWCx) = Jac_CWC(2, idxJacFx);
            cwcGrad3(ii, idxDecVarCWCy) = Jac_CWC(3, idxJacFy);
            cwcGrad4(ii, idxDecVarCWCy) = Jac_CWC(4, idxJacFy);
            
            dtGrad = (ii - 1)/(Nt - 1) * ones(1, 5) - [1, 1, 1, 1, 0];
            c2Grad(ii, 1:5) = c2Grad(ii, 1:5) + Jac_fz(end) * dtGrad;
            cwcGrad1(ii, 1:5) = cwcGrad1(ii, 1:5) + Jac_CWC(1, end) * dtGrad;
            cwcGrad2(ii, 1:5) = cwcGrad2(ii, 1:5) + Jac_CWC(2, end) * dtGrad;
            cwcGrad3(ii, 1:5) = cwcGrad3(ii, 1:5) + Jac_CWC(3, end) * dtGrad;
            cwcGrad4(ii, 1:5) = cwcGrad4(ii, 1:5) + Jac_CWC(4, end) * dtGrad;

        elseif tt > T5/3 + tol && tt <= T5 * 2 / 3 + tol
            
            t_temp = tt - T5/3;
            % -fz <= 0
            idxFz = 64 + 12 + 12 + 2 + (1:4);            
            Jac_fz = getJacobianFPzCst(T5/3,f(idxFz),t_temp);
            
            Jac_fz(1) = Jac_fz(1)/3;
            
            idxDecVarFz = [5, 35+26+4+4+(1:4)];
            idxJacFz = [1, 1+(1:4)];
            c2Grad(ii, idxDecVarFz) = Jac_fz(idxJacFz);  
            
            % CWC
            idxCWC = 12 + 12 + 2 + [1:4, 32+(1:4), 64+(1:4)];            
            Jac_CWC = getJacobianCWCCst(T5/3, f(idxCWC), miu, t_temp);
            
            Jac_CWC(:, 1) = Jac_CWC(:, 1)/3;
            
            idxDecVarCWCx = [5,    35+4+4+(1:4), 35+4+4+26+(1:4)];
            idxDecVarCWCy = [5, 35+4+4+13+(1:4), 35+4+4+26+(1:4)];
            
            idxJacFx = [1,   1+(1:4), 1+8+(1:4)];
            idxJacFy = [1, 1+4+(1:4), 1+8+(1:4)];
            
            cwcGrad1(ii, idxDecVarCWCx) = Jac_CWC(1, idxJacFx);
            cwcGrad2(ii, idxDecVarCWCx) = Jac_CWC(2, idxJacFx);
            cwcGrad3(ii, idxDecVarCWCy) = Jac_CWC(3, idxJacFy);
            cwcGrad4(ii, idxDecVarCWCy) = Jac_CWC(4, idxJacFy);
            
            dtGrad = (ii - 1)/(Nt - 1) * ones(1, 5) - [1, 1, 1, 1, 1/3];
            c2Grad(ii, 1:5) = c2Grad(ii, 1:5) + Jac_fz(end) * dtGrad;
            cwcGrad1(ii, 1:5) = cwcGrad1(ii, 1:5) + Jac_CWC(1, end) * dtGrad;
            cwcGrad2(ii, 1:5) = cwcGrad2(ii, 1:5) + Jac_CWC(2, end) * dtGrad;
            cwcGrad3(ii, 1:5) = cwcGrad3(ii, 1:5) + Jac_CWC(3, end) * dtGrad;
            cwcGrad4(ii, 1:5) = cwcGrad4(ii, 1:5) + Jac_CWC(4, end) * dtGrad;
            
        else
            
            t_temp = tt - 2 * T5/3;
            % -fz <= 0
            idxFz = 64 + 12 + 12 + 2 + 2 + (1:4);            
            Jac_fz = getJacobianFPzCst(T5/3,f(idxFz),t_temp);
            
            Jac_fz(1) = Jac_fz(1)/3;
            
            idxDecVarFz = [5, 35+26+4+4+2+(1:3)];
            idxJacFz = [1, 1+(1:3)];
            c2Grad(ii, idxDecVarFz) = Jac_fz(idxJacFz);  
            
            % CWC
            idxCWC = 12 + 12 + 2 + 2 + [1:4, 32+(1:4), 64+(1:4)];            
            Jac_CWC = getJacobianCWCCst(T5/3, f(idxCWC), miu, t_temp);
            
            Jac_CWC(:, 1) = Jac_CWC(:, 1)/3;
            
            idxDecVarCWCx = [5,    35+4+4+2+(1:3), 35+4+4+2+26+(1:3)];
            idxDecVarCWCy = [5, 35+4+4+2+13+(1:3), 35+4+4+2+26+(1:3)];
            
            idxJacFx = [1,   1+(1:3), 1+8+(1:3)];
            idxJacFy = [1, 1+4+(1:3), 1+8+(1:3)];
            
            cwcGrad1(ii, idxDecVarCWCx) = Jac_CWC(1, idxJacFx);
            cwcGrad2(ii, idxDecVarCWCx) = Jac_CWC(2, idxJacFx);
            cwcGrad3(ii, idxDecVarCWCy) = Jac_CWC(3, idxJacFy);
            cwcGrad4(ii, idxDecVarCWCy) = Jac_CWC(4, idxJacFy);
            
            dtGrad = (ii - 1)/(Nt - 1) * ones(1, 5) - [1, 1, 1, 1, 2/3];
            c2Grad(ii, 1:5) = c2Grad(ii, 1:5) + Jac_fz(end) * dtGrad;
            cwcGrad1(ii, 1:5) = cwcGrad1(ii, 1:5) + Jac_CWC(1, end) * dtGrad;
            cwcGrad2(ii, 1:5) = cwcGrad2(ii, 1:5) + Jac_CWC(2, end) * dtGrad;
            cwcGrad3(ii, 1:5) = cwcGrad3(ii, 1:5) + Jac_CWC(3, end) * dtGrad;
            cwcGrad4(ii, 1:5) = cwcGrad4(ii, 1:5) + Jac_CWC(4, end) * dtGrad;
        end
    end
    
end

% 整合——要与pathCst里的顺序匹配
cGrad = [c1Grad; c2Grad; cwcGrad1; cwcGrad2; cwcGrad3; cwcGrad4; c3Grad];
ceqGrad = [];

end