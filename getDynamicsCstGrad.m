function [c, ceq, cGrad, ceqGrad, ddx] = getDynamicsCstGrad(z, pack)

global param

[tPhase, p, f, x] = unPackDecVar(z, pack);

T1 = tPhase(1);      T2 = tPhase(2);      T3 = tPhase(3);      T4 = tPhase(4);      T5 = tPhase(5);

T = sum(tPhase);
Nt = size(x, 2);

dt = T / (Nt - 1);

nDecVar = numel(z);

% param.m = 20.0;
% param.g = 9.81;
% param.I = diag([0.25, 0.18, 0.16]);

tol = 1e-6;

idxLow = 1:(Nt - 1);
idxUpp = 2:Nt;

dynGradRx = zeros(2 * (Nt - 1), nDecVar);
dynGradRy = zeros(2 * (Nt - 1), nDecVar);
dynGradRz = zeros(2 * (Nt - 1), nDecVar);

dynGradQx = zeros(2 * (Nt - 1), nDecVar);
dynGradQy = zeros(2 * (Nt - 1), nDecVar);
dynGradQz = zeros(2 * (Nt - 1), nDecVar);

ddrxBase = zeros(2 * (Nt - 1), 1);
ddryBase = zeros(2 * (Nt - 1), 1);
ddrzBase = zeros(2 * (Nt - 1), 1);

ddqxBase = zeros(2 * (Nt - 1), 1);
ddqyBase = zeros(2 * (Nt - 1), 1);
ddqzBase = zeros(2 * (Nt - 1), 1);

for ii = 1:(Nt-1)
    t = (ii - 1) * dt;
    
    xLow = x(:, idxLow(ii));
    xUpp = x(:, idxUpp(ii));
    
    if t <= T1 + tol
        tt = t;
        idxP = [1:4, 28+(1:4), 56+(1:4)]; % p只用T1的四个值
        T_temp = T1/3;
        if tt <= T1/3 + tol  
            t_temp = tt;
            
            idxF = [1:4, 32+(1:4), 64+(1:4)]; % f在T1阶段的前4个值
            
            % idxJac = [fx3, fx4, fy3, fy4, fz3, fz4, xLow, xUpp]
            % 取出只与Decision Variable里相关量的梯度
            idxJac = 1 + [12+(3:4), 12+4+(3:4), 12+8+(3:4), 24+(1:24)];
            
            % idxDecVar = [fx1, fx2, fy1, fy2, fz1, fz2, xLow, xUpp]
            % decision variables里对应jacobian的index数据的位置
            idxDecVar = 5 + [30+(1:2), 30+13+(1:2), 30+26+(1:2), (58+12*ii):(81+12*ii)];          
            
            % 初始时刻
            Jac_ddx_0 = getJacobianAcc(T1/3, p(idxP), f(idxF), xLow, xUpp, param, t_temp, 0, dt);
            % 结束时刻
            Jac_ddx_dt = getJacobianAcc(T1/3, p(idxP), f(idxF), xLow, xUpp, param, t_temp, dt, dt);
            
            % t, T相对于[T1, T2, T3, T4, T5]的偏导数
            dTGrad = 1/3 * [1, 0, 0, 0, 0];
            dtGrad = (ii - 1)/(Nt - 1) * ones(1, 5) - [0, 0, 0, 0, 0];
                         
        elseif tt > T1/3 + tol && tt <= T1 * 2 / 3 + tol
            idxF = 2 + [1:4, 32+(1:4), 64+(1:4)]; % % f在T1阶段的中间4个值——序号+2
            
            idxJac = 1 + [12+(1:4), 12+4+(1:4), 12+8+(1:4), 24+(1:24)];
            
            idxDecVar = 5 + [30+(1:4), 30+13+(1:4), 30+26+(1:4), (58+12*ii):(81+12*ii)];
                        
            t_temp = tt - T1/3;
            % 初始时刻
            Jac_ddx_0 = getJacobianAcc(T1/3, p(idxP), f(idxF), xLow, xUpp, param, t_temp, 0, dt);
            % 结束时刻
            Jac_ddx_dt = getJacobianAcc(T1/3, p(idxP), f(idxF), xLow, xUpp, param, t_temp, dt, dt);
            
            % t, T相对于[T1, T2, T3, T4, T5]的偏导数
            dTGrad = 1/3 * [1, 0, 0, 0, 0];
            dtGrad = (ii - 1)/(Nt - 1) * ones(1, 5) - [1/3, 0, 0, 0, 0];
            
        else
            idxF = 4 + [1:4, 32+(1:4), 64+(1:4)]; % f在T1阶段的后4个值——序号+2
            
            idxJac = 1 + [12+(1:2), 12+4+(1:2), 12+8+(1:2), 24+(1:24)];
         
            idxDecVar = 5 + [30+(3:4), 30+13+(3:4), 30+26+(3:4), (58+12*ii):(81+12*ii)];            
            
            t_temp = tt - 2 * T1/3;
            % 初始时刻
            Jac_ddx_0 = getJacobianAcc(T1/3, p(idxP), f(idxF), xLow, xUpp, param, t_temp, 0, dt);
            % 结束时刻
            Jac_ddx_dt = getJacobianAcc(T1/3, p(idxP), f(idxF), xLow, xUpp, param, t_temp, dt, dt);
            
            % t, T相对于[T1, T2, T3, T4, T5]的偏导数
            dTGrad = 1/3 * [1, 0, 0, 0, 0];
            dtGrad = (ii - 1)/(Nt - 1) * ones(1, 5) - [2/3, 0, 0, 0, 0];
                        
        end
        
    elseif t > T1 + tol && t <= T1 + T2 + tol
        tt = t - T1;
        idxF = 8 + [1:4, 32+(1:4), 64+(1:4)]; % f只用T2的4个值
        T_temp = T2/3;
        if tt <= T2/3 + tol
            idxP = 4 + [1:4, 28+(1:4), 56+(1:4)];
            
            idxJac = 1 + [3:4, 4+(3:4), 8+(3:4), 24+(1:24)];
            
            idxDecVar = 5 + [1:2, 10+(1:2), 20+(1:2), (58+12*ii):(81+12*ii)];
        
            t_temp = tt;
            % 初始时刻
            Jac_ddx_0 = getJacobianAcc(T2/3, p(idxP), f(idxF), xLow, xUpp, param, t_temp, 0, dt);
            % 结束时刻
            Jac_ddx_dt = getJacobianAcc(T2/3, p(idxP), f(idxF), xLow, xUpp, param, t_temp, dt, dt);
            
            % t, T相对于[T1, T2, T3, T4, T5]的偏导数
            dTGrad = 1/3 * [0, 1, 0, 0, 0];
            dtGrad = (ii - 1)/(Nt - 1) * ones(1, 5) - [1, 0, 0, 0, 0];
            
        elseif tt > T2/3 + tol && tt <= T2 * 2 / 3 + tol
            % idxP = 4 + [1:4, 28+(1:4), 56+(1:4)]; % 上面的顺序+2
            idxP = 4 + 2 + [1:4, 28+(1:4), 56+(1:4)];
            
            idxJac = 1 + [1:4, 4+(1:4), 8+(1:4), 24+(1:24)];
            
            idxDecVar = 5 + [1:4, 10+(1:4), 20+(1:4), (58+12*ii):(81+12*ii)];           
         
            t_temp = tt - T2/3;
            % 初始时刻
            Jac_ddx_0 = getJacobianAcc(T2/3, p(idxP), f(idxF), xLow, xUpp, param, t_temp, 0, dt);
            % 结束时刻
            Jac_ddx_dt = getJacobianAcc(T2/3, p(idxP), f(idxF), xLow, xUpp, param, t_temp, dt, dt);
            
            % t, T相对于[T1, T2, T3, T4, T5]的偏导数
            dTGrad = 1/3 * [0, 1, 0, 0, 0];
            dtGrad = (ii - 1)/(Nt - 1) * ones(1, 5) - [1, 1/3, 0, 0, 0];
            
        else
            % idxP = 4 + 2 + [1:4, 28+(1:4), 56+(1:4)];
            idxP = 4 + 4 + [1:4, 28+(1:4), 56+(1:4)];
            
            idxJac = 1 + [1:3, 4+(1:3), 8+(1:3), 24+(1:24)]; 
            
            idxDecVar = 5 + [2+(1:3), 10+2+(1:3), 20+2+(1:3), (58+12*ii):(81+12*ii)];           
   
            t_temp = tt - 2 * T2/3;
            % 初始时刻
            Jac_ddx_0 = getJacobianAcc(T2/3, p(idxP), f(idxF), xLow, xUpp, param, t_temp, 0, dt);
            % 结束时刻
            Jac_ddx_dt = getJacobianAcc(T2/3, p(idxP), f(idxF), xLow, xUpp, param, t_temp, dt, dt);
            
            % t, T相对于[T1, T2, T3, T4, T5]的偏导数
            dTGrad = 1/3 * [0, 1, 0, 0, 0];
            dtGrad = (ii - 1)/(Nt - 1) * ones(1, 5) - [1, 2/3, 0, 0, 0];
            
        end
        
    elseif t > T1 + T2 + tol && t <= T1 + T2 + T3 + tol
        tt = t - T1 - T2;
        idxP = 12 + [1:4, 28+(1:4), 56+(1:4)];
        T_temp = T3/3;
        if tt <= T3/3 + tol
            idxF = 12 + [1:4, 32+(1:4), 64+(1:4)];
            
            % 这里需要注意一下：脚的运动与决策变量有关
            idxJac_nom = 1 + [1, 3, 5, 7, 9, 11, 12+(3:4), 12+4+(3:4), 12+8+(3:4), 24+(1:24)];
            
            idxDecVar = 5 + [5, 15, 25, 30+4+(1:2), 30+4+13+(1:2), 30+4+26+(1:2), (58+12*ii):(81+12*ii)];            
            
            t_temp = tt;
            % 初始时刻
            Jac_ddx_0 = getJacobianAcc(T3/3, p(idxP), f(idxF), xLow, xUpp, param, t_temp, 0, dt);
            % 结束时刻
            Jac_ddx_dt = getJacobianAcc(T3/3, p(idxP), f(idxF), xLow, xUpp, param, t_temp, dt, dt);
            
            % 需要对Jacobian矩阵进行处理
            Jac_ddx_0(:, idxJac_nom(1)) = Jac_ddx_0(:, idxJac_nom(1)) + Jac_ddx_0(:, idxJac_nom(2));
            Jac_ddx_0(:, idxJac_nom(3)) = Jac_ddx_0(:, idxJac_nom(3)) + Jac_ddx_0(:, idxJac_nom(4));
            Jac_ddx_0(:, idxJac_nom(5)) = Jac_ddx_0(:, idxJac_nom(5)) + Jac_ddx_0(:, idxJac_nom(6));
            
            Jac_ddx_dt(:, idxJac_nom(1)) = Jac_ddx_dt(:, idxJac_nom(1)) + Jac_ddx_dt(:, idxJac_nom(2));
            Jac_ddx_dt(:, idxJac_nom(3)) = Jac_ddx_dt(:, idxJac_nom(3)) + Jac_ddx_dt(:, idxJac_nom(4));
            Jac_ddx_dt(:, idxJac_nom(5)) = Jac_ddx_dt(:, idxJac_nom(5)) + Jac_ddx_dt(:, idxJac_nom(6));
            
            idxJac = 1 + [1, 5, 9, 12+(3:4), 12+4+(3:4), 12+8+(3:4), 24+(1:24)];
            
            % t, T相对于[T1, T2, T3, T4, T5]的偏导数
            dTGrad = 1/3 * [0, 0, 1, 0, 0];
            dtGrad = (ii - 1)/(Nt - 1) * ones(1, 5) - [1, 1, 0, 0, 0];
            
        elseif tt > T3/3 + tol && tt <= 2 * T3/3 + tol
            % idxF = 12 + [1:4, 32+(1:4), 64+(1:4)];
            idxF = 12 + 2 + [1:4, 32+(1:4), 64+(1:4)];
            
            idxJac_nom = 1 + [1, 3, 5, 7, 9, 11, 12+(1:4), 12+4+(1:4), 12+8+(1:4), 24+(1:24)];
            
            idxDecVar = 5 + [5, 15, 25, 30+4+(1:4), 30+4+13+(1:4), 30+4+26+(1:4), (58+12*ii):(81+12*ii)];
            
            t_temp = tt - T3/3;
            % 初始时刻
            Jac_ddx_0 = getJacobianAcc(T3/3, p(idxP), f(idxF), xLow, xUpp, param, t_temp, 0, dt);
            % 结束时刻
            Jac_ddx_dt = getJacobianAcc(T3/3, p(idxP), f(idxF), xLow, xUpp, param, t_temp, dt, dt);
            
            % 需要对Jacobian矩阵进行处理
            Jac_ddx_0(:, idxJac_nom(1)) = Jac_ddx_0(:, idxJac_nom(1)) + Jac_ddx_0(:, idxJac_nom(2));
            Jac_ddx_0(:, idxJac_nom(3)) = Jac_ddx_0(:, idxJac_nom(3)) + Jac_ddx_0(:, idxJac_nom(4));
            Jac_ddx_0(:, idxJac_nom(5)) = Jac_ddx_0(:, idxJac_nom(5)) + Jac_ddx_0(:, idxJac_nom(6));
            
            Jac_ddx_dt(:, idxJac_nom(1)) = Jac_ddx_dt(:, idxJac_nom(1)) + Jac_ddx_dt(:, idxJac_nom(2));
            Jac_ddx_dt(:, idxJac_nom(3)) = Jac_ddx_dt(:, idxJac_nom(3)) + Jac_ddx_dt(:, idxJac_nom(4));
            Jac_ddx_dt(:, idxJac_nom(5)) = Jac_ddx_dt(:, idxJac_nom(5)) + Jac_ddx_dt(:, idxJac_nom(6));
            
            idxJac = 1 + [1, 5, 9, 12+(1:4), 12+4+(1:4), 12+8+(1:4), 24+(1:24)];
            
            % t, T相对于[T1, T2, T3, T4, T5]的偏导数
            dTGrad = 1/3 * [0, 0, 1, 0, 0];
            dtGrad = (ii - 1)/(Nt - 1) * ones(1, 5) - [1, 1, 1/3, 0, 0];
            
        else
            % idxF = 12 + 2 + [1:4, 32+(1:4), 64+(1:4)];
            idxF = 12 + 2 + 2 + [1:4, 32+(1:4), 64+(1:4)];
            
            idxJac_nom = 1 + [1, 3, 5, 7, 9, 11, 12+(1:2), 12+4+(1:2), 12+8+(1:2), 24+(1:24)];
            
            idxDecVar = 5 + [5, 15, 25, 30+4+(3:4), 30+4+13+(3:4), 30+4+26+(3:4), (58+12*ii):(81+12*ii)];
            
            t_temp = tt - 2 * T3/3;
            % 初始时刻
            Jac_ddx_0 = getJacobianAcc(T3/3, p(idxP), f(idxF), xLow, xUpp, param, t_temp, 0, dt);
            % 结束时刻
            Jac_ddx_dt = getJacobianAcc(T3/3, p(idxP), f(idxF), xLow, xUpp, param, t_temp, dt, dt);
            
            % 需要对Jacobian矩阵进行处理
            Jac_ddx_0(:, idxJac_nom(1)) = Jac_ddx_0(:, idxJac_nom(1)) + Jac_ddx_0(:, idxJac_nom(2));
            Jac_ddx_0(:, idxJac_nom(3)) = Jac_ddx_0(:, idxJac_nom(3)) + Jac_ddx_0(:, idxJac_nom(4));
            Jac_ddx_0(:, idxJac_nom(5)) = Jac_ddx_0(:, idxJac_nom(5)) + Jac_ddx_0(:, idxJac_nom(6));
            
            Jac_ddx_dt(:, idxJac_nom(1)) = Jac_ddx_dt(:, idxJac_nom(1)) + Jac_ddx_dt(:, idxJac_nom(2));
            Jac_ddx_dt(:, idxJac_nom(3)) = Jac_ddx_dt(:, idxJac_nom(3)) + Jac_ddx_dt(:, idxJac_nom(4));
            Jac_ddx_dt(:, idxJac_nom(5)) = Jac_ddx_dt(:, idxJac_nom(5)) + Jac_ddx_dt(:, idxJac_nom(6));
            
            idxJac = 1 + [1, 5, 9, 12+(1:2), 12+4+(1:2), 12+8+(1:2), 24+(1:24)];
            
            % t, T相对于[T1, T2, T3, T4, T5]的偏导数
            dTGrad = 1/3 * [0, 0, 1, 0, 0];
            dtGrad = (ii - 1)/(Nt - 1) * ones(1, 5) - [1, 1, 2/3, 0, 0];
            
        end
    elseif t > T1 + T2 + T3 + tol && t <= T1 + T2 + T3 + T4 + tol 
        tt = t - T1 - T2 - T3;
        idxF = 8 + 12 + [1:4, 32+(1:4), 64+(1:4)];
        T_temp = T4/3;
        if tt <= T4/3 + tol
            idxP = 4 + 12 + [1:4, 28+(1:4), 56+(1:4)];
            
            idxJac = 1 + [1, 3:4, 4+1, 4+(3:4), 8+1, 8+(3:4), 24+(1:24)];

            idxDecVar = 5 + [5, 5+(1:2), 10+5, 10+5+(1:2), 20+5, 20+5+(1:2), (58+12*ii):(81+12*ii)];         
            
            t_temp = tt;
            % 初始时刻
            Jac_ddx_0 = getJacobianAcc(T4/3, p(idxP), f(idxF), xLow, xUpp, param, t_temp, 0, dt);
            % 结束时刻
            Jac_ddx_dt = getJacobianAcc(T4/3, p(idxP), f(idxF), xLow, xUpp, param, t_temp, dt, dt);
            
            % t, T相对于[T1, T2, T3, T4, T5]的偏导数
            dTGrad = 1/3 * [0, 0, 0, 1, 0];
            dtGrad = (ii - 1)/(Nt - 1) * ones(1, 5) - [1, 1, 1, 0, 0];
            
        elseif tt > T4/3 + tol && tt <= 2 * T4/3 + tol
            idxP = 4 + 12 + 2 + [1:4, 28+(1:4), 56+(1:4)];
            
            idxJac = 1 + [1:4, 4+(1:4), 8+(1:4), 24+(1:24)];
            
            idxDecVar = 5 + [5+(1:4), 10+5+(1:4), 20+5+(1:4), (58+12*ii):(81+12*ii)];
            
            t_temp = tt - T4/3;
            % 初始时刻
            Jac_ddx_0 = getJacobianAcc(T4/3, p(idxP), f(idxF), xLow, xUpp, param, t_temp, 0, dt);
            % 结束时刻
            Jac_ddx_dt = getJacobianAcc(T4/3, p(idxP), f(idxF), xLow, xUpp, param, t_temp, dt, dt);
            
            % t, T相对于[T1, T2, T3, T4, T5]的偏导数
            dTGrad = 1/3 * [0, 0, 0, 1, 0];
            dtGrad = (ii - 1)/(Nt - 1) * ones(1, 5) - [1, 1, 1, 1/3, 0];
            
        else
            idxP = 4 + 12 + 4 + [1:4, 28+(1:4), 56+(1:4)];
            
            idxJac = 1 + [1:3, 4+(1:3), 8+(1:3), 24+(1:24)];
            
            idxDecVar = 5 + [5+2+(1:3), 10+5+2+(1:3), 20+5+2+(1:3), (58+12*ii):(81+12*ii)];
            
            t_temp = tt - 2 * T4/3;
            % 初始时刻
            Jac_ddx_0 = getJacobianAcc(T4/3, p(idxP), f(idxF), xLow, xUpp, param, t_temp, 0, dt);
            % 结束时刻
            Jac_ddx_dt = getJacobianAcc(T4/3, p(idxP), f(idxF), xLow, xUpp, param, t_temp, dt, dt);
            
            % t, T相对于[T1, T2, T3, T4, T5]的偏导数
            dTGrad = 1/3 * [0, 0, 0, 1, 0];
            dtGrad = (ii - 1)/(Nt - 1) * ones(1, 5) - [1, 1, 1, 2/3, 0];
            
        end
    else
        tt = t - T1 - T2 - T3 - T4;
        idxP = 12 + 12 + [1:4, 28+(1:4), 56+(1:4)];
        T_temp = T5/3;
        if tt <= T5/3 + tol
            idxF = 12 + 12 + [1:4, 32+(1:4), 64+(1:4)];
            
            idxJac_nom = 1 + [1, 3, 5, 7, 9, 11, 12+(3:4), 12+4+(3:4), 12+8+(3:4), 24+(1:24)];
            
            idxDecVar = 5 + [5+5, 15+5, 25+5, 30+4+4+(1:2), 30+4+4+13+(1:2), 30+4+4+26+(1:2), (58+12*ii):(81+12*ii)];
            
            t_temp = tt;
            % 初始时刻
            Jac_ddx_0 = getJacobianAcc(T5/3, p(idxP), f(idxF), xLow, xUpp, param, t_temp, 0, dt);
            % 结束时刻
            Jac_ddx_dt = getJacobianAcc(T5/3, p(idxP), f(idxF), xLow, xUpp, param, t_temp, dt, dt);
            
            % 需要对Jacobian矩阵进行处理
            Jac_ddx_0(:, idxJac_nom(1)) = Jac_ddx_0(:, idxJac_nom(1)) + Jac_ddx_0(:, idxJac_nom(2));
            Jac_ddx_0(:, idxJac_nom(3)) = Jac_ddx_0(:, idxJac_nom(3)) + Jac_ddx_0(:, idxJac_nom(4));
            Jac_ddx_0(:, idxJac_nom(5)) = Jac_ddx_0(:, idxJac_nom(5)) + Jac_ddx_0(:, idxJac_nom(6));
            
            Jac_ddx_dt(:, idxJac_nom(1)) = Jac_ddx_dt(:, idxJac_nom(1)) + Jac_ddx_dt(:, idxJac_nom(2));
            Jac_ddx_dt(:, idxJac_nom(3)) = Jac_ddx_dt(:, idxJac_nom(3)) + Jac_ddx_dt(:, idxJac_nom(4));
            Jac_ddx_dt(:, idxJac_nom(5)) = Jac_ddx_dt(:, idxJac_nom(5)) + Jac_ddx_dt(:, idxJac_nom(6));
            
            idxJac = 1 + [1, 5, 9, 12+(3:4), 12+4+(3:4), 12+8+(3:4), 24+(1:24)];
            
            % t, T相对于[T1, T2, T3, T4, T5]的偏导数
            dTGrad = 1/3 * [0, 0, 0, 0, 1];
            dtGrad = (ii - 1)/(Nt - 1) * ones(1, 5) - [1, 1, 1, 1, 0];
            
        elseif tt > T5/3 + tol && tt <= 2 * T5/3 + tol
            idxF = 12 + 12 + 2 + [1:4, 32+(1:4), 64+(1:4)];
            
            idxJac_nom = 1 + [1, 3, 5, 7, 9, 11, 12+(1:4), 12+4+(1:4), 12+8+(1:4), 24+(1:24)];
            
            idxDecVar = 5 + [5+5, 15+5, 25+5, 30+4+4+(1:4), 30+4+4+13+(1:4), 30+4+4+26+(1:4), (58+12*ii):(81+12*ii)];
            
            t_temp = tt - T5/3;
            % 初始时刻
            Jac_ddx_0 = getJacobianAcc(T5/3, p(idxP), f(idxF), xLow, xUpp, param, t_temp, 0, dt);
            % 结束时刻
            Jac_ddx_dt = getJacobianAcc(T5/3, p(idxP), f(idxF), xLow, xUpp, param, t_temp, dt, dt);
            
            % 需要对Jacobian矩阵进行处理
            Jac_ddx_0(:, idxJac_nom(1)) = Jac_ddx_0(:, idxJac_nom(1)) + Jac_ddx_0(:, idxJac_nom(2));
            Jac_ddx_0(:, idxJac_nom(3)) = Jac_ddx_0(:, idxJac_nom(3)) + Jac_ddx_0(:, idxJac_nom(4));
            Jac_ddx_0(:, idxJac_nom(5)) = Jac_ddx_0(:, idxJac_nom(5)) + Jac_ddx_0(:, idxJac_nom(6));
            
            Jac_ddx_dt(:, idxJac_nom(1)) = Jac_ddx_dt(:, idxJac_nom(1)) + Jac_ddx_dt(:, idxJac_nom(2));
            Jac_ddx_dt(:, idxJac_nom(3)) = Jac_ddx_dt(:, idxJac_nom(3)) + Jac_ddx_dt(:, idxJac_nom(4));
            Jac_ddx_dt(:, idxJac_nom(5)) = Jac_ddx_dt(:, idxJac_nom(5)) + Jac_ddx_dt(:, idxJac_nom(6));
            
            idxJac = 1 + [1, 5, 9, 12+(1:4), 12+4+(1:4), 12+8+(1:4), 24+(1:24)];
            
            % t, T相对于[T1, T2, T3, T4, T5]的偏导数
            dTGrad = 1/3 * [0, 0, 0, 0, 1];
            dtGrad = (ii - 1)/(Nt - 1) * ones(1, 5) - [1, 1, 1, 1, 1/3];
            
        else
            idxF = 12 + 12 + 2 + 2 + [1:4, 32+(1:4), 64+(1:4)];
            
            idxJac_nom = 1 + [1, 3, 5, 7, 9, 11, 12+(1:3), 12+4+(1:3), 12+8+(1:3), 24+(1:24)];
            
            idxDecVar = 5 + [5+5, 15+5, 25+5, 30+4+4+2+(1:3), 30+4+4+2+13+(1:3), 30+4+4+2+26+(1:3), (58+12*ii):(81+12*ii)];
            
            t_temp = tt - 2 * T5/3;
            % 初始时刻
            Jac_ddx_0 = getJacobianAcc(T5/3, p(idxP), f(idxF), xLow, xUpp, param, t_temp, 0, dt);
            % 结束时刻
            Jac_ddx_dt = getJacobianAcc(T5/3, p(idxP), f(idxF), xLow, xUpp, param, t_temp, dt, dt);
            
            % 需要对Jacobian矩阵进行处理
            Jac_ddx_0(:, idxJac_nom(1)) = Jac_ddx_0(:, idxJac_nom(1)) + Jac_ddx_0(:, idxJac_nom(2));
            Jac_ddx_0(:, idxJac_nom(3)) = Jac_ddx_0(:, idxJac_nom(3)) + Jac_ddx_0(:, idxJac_nom(4));
            Jac_ddx_0(:, idxJac_nom(5)) = Jac_ddx_0(:, idxJac_nom(5)) + Jac_ddx_0(:, idxJac_nom(6));
            
            Jac_ddx_dt(:, idxJac_nom(1)) = Jac_ddx_dt(:, idxJac_nom(1)) + Jac_ddx_dt(:, idxJac_nom(2));
            Jac_ddx_dt(:, idxJac_nom(3)) = Jac_ddx_dt(:, idxJac_nom(3)) + Jac_ddx_dt(:, idxJac_nom(4));
            Jac_ddx_dt(:, idxJac_nom(5)) = Jac_ddx_dt(:, idxJac_nom(5)) + Jac_ddx_dt(:, idxJac_nom(6));
            
            idxJac = 1 + [1, 5, 9, 12+(1:3), 12+4+(1:3), 12+8+(1:3), 24+(1:24)];
            
            % t, T相对于[T1, T2, T3, T4, T5]的偏导数
            dTGrad = 1/3 * [0, 0, 0, 0, 1];
            dtGrad = (ii - 1)/(Nt - 1) * ones(1, 5) - [1, 1, 1, 1, 2/3];
        end
    end
    
    if 1
    Base_ddx_0 = getBaseAcc(T_temp, p(idxP), f(idxF), xLow, xUpp, param, t_temp, 0, dt);
    Base_ddx_dt = getBaseAcc(T_temp, p(idxP), f(idxF), xLow, xUpp, param, t_temp, dt, dt);
    
    ddrxBase(2 * ii - 1) = Base_ddx_0(1);       ddrxBase(2 * ii) = Base_ddx_dt(1);
    ddryBase(2 * ii - 1) = Base_ddx_0(2);       ddryBase(2 * ii) = Base_ddx_dt(2);
    ddrzBase(2 * ii - 1) = Base_ddx_0(3);       ddrzBase(2 * ii) = Base_ddx_dt(3);
    ddqxBase(2 * ii - 1) = Base_ddx_0(4);       ddqxBase(2 * ii) = Base_ddx_dt(4);
    ddqyBase(2 * ii - 1) = Base_ddx_0(5);       ddqyBase(2 * ii) = Base_ddx_dt(5);
    ddqzBase(2 * ii - 1) = Base_ddx_0(6);       ddqzBase(2 * ii) = Base_ddx_dt(6);
    
    dynGradRx(2 * ii - 1, idxDecVar) = Jac_ddx_0(1, idxJac);
    dynGradRy(2 * ii - 1, idxDecVar) = Jac_ddx_0(2, idxJac);
    dynGradRz(2 * ii - 1, idxDecVar) = Jac_ddx_0(3, idxJac);            
    dynGradQx(2 * ii - 1, idxDecVar) = Jac_ddx_0(4, idxJac);
    dynGradQy(2 * ii - 1, idxDecVar) = Jac_ddx_0(5, idxJac);
    dynGradQz(2 * ii - 1, idxDecVar) = Jac_ddx_0(6, idxJac);

    dynGradRx(2 * ii, idxDecVar) = Jac_ddx_dt(1, idxJac);
    dynGradRy(2 * ii, idxDecVar) = Jac_ddx_dt(2, idxJac);
    dynGradRz(2 * ii, idxDecVar) = Jac_ddx_dt(3, idxJac);            
    dynGradQx(2 * ii, idxDecVar) = Jac_ddx_dt(4, idxJac);
    dynGradQy(2 * ii, idxDecVar) = Jac_ddx_dt(5, idxJac);
    dynGradQz(2 * ii, idxDecVar) = Jac_ddx_dt(6, idxJac);
    end
    
    dt1Grad = 1/(Nt - 1) * [1, 1, 1, 1, 1];
    dT1Grad = 1/(Nt - 1) * [1, 1, 1, 1, 1];

    dynGradRx(2 * ii - 1, 1:5) = Jac_ddx_0(1, 1) * dTGrad + Jac_ddx_0(1, end-2) * dtGrad + Jac_ddx_0(1, end) * dT1Grad;
    dynGradRy(2 * ii - 1, 1:5) = Jac_ddx_0(2, 1) * dTGrad + Jac_ddx_0(2, end-2) * dtGrad + Jac_ddx_0(2, end) * dT1Grad;
    dynGradRz(2 * ii - 1, 1:5) = Jac_ddx_0(3, 1) * dTGrad + Jac_ddx_0(3, end-2) * dtGrad + Jac_ddx_0(3, end) * dT1Grad;
    dynGradQx(2 * ii - 1, 1:5) = Jac_ddx_0(4, 1) * dTGrad + Jac_ddx_0(4, end-2) * dtGrad + Jac_ddx_0(4, end) * dT1Grad;
    dynGradQy(2 * ii - 1, 1:5) = Jac_ddx_0(5, 1) * dTGrad + Jac_ddx_0(5, end-2) * dtGrad + Jac_ddx_0(5, end) * dT1Grad;
    dynGradQz(2 * ii - 1, 1:5) = Jac_ddx_0(6, 1) * dTGrad + Jac_ddx_0(6, end-2) * dtGrad + Jac_ddx_0(6, end) * dT1Grad;

    dynGradRx(2 * ii, 1:5) = Jac_ddx_dt(1, 1) * dTGrad + Jac_ddx_dt(1, end-2) * dtGrad + Jac_ddx_dt(1, end-1) * dt1Grad + Jac_ddx_dt(1, end) * dT1Grad;
    dynGradRy(2 * ii, 1:5) = Jac_ddx_dt(2, 1) * dTGrad + Jac_ddx_dt(2, end-2) * dtGrad + Jac_ddx_dt(2, end-1) * dt1Grad + Jac_ddx_dt(2, end) * dT1Grad;
    dynGradRz(2 * ii, 1:5) = Jac_ddx_dt(3, 1) * dTGrad + Jac_ddx_dt(3, end-2) * dtGrad + Jac_ddx_dt(3, end-1) * dt1Grad + Jac_ddx_dt(3, end) * dT1Grad;
    dynGradQx(2 * ii, 1:5) = Jac_ddx_dt(4, 1) * dTGrad + Jac_ddx_dt(4, end-2) * dtGrad + Jac_ddx_dt(4, end-1) * dt1Grad + Jac_ddx_dt(4, end) * dT1Grad;
    dynGradQy(2 * ii, 1:5) = Jac_ddx_dt(5, 1) * dTGrad + Jac_ddx_dt(5, end-2) * dtGrad + Jac_ddx_dt(5, end-1) * dt1Grad + Jac_ddx_dt(5, end) * dT1Grad;
    dynGradQz(2 * ii, 1:5) = Jac_ddx_dt(6, 1) * dTGrad + Jac_ddx_dt(6, end-2) * dtGrad + Jac_ddx_dt(6, end-1) * dt1Grad + Jac_ddx_dt(6, end) * dT1Grad;

end

%% 以上便得到了动力学中所有ddx相对于decision variables的梯度，接下来是对这些梯度进行处理以对应于约束的顺序
% 先来整rx方向的梯度，之后同理
ceq1GradRx = zeros(Nt - 2, nDecVar);
ceq1GradRy = zeros(Nt - 2, nDecVar);
ceq1GradRz = zeros(Nt - 2, nDecVar);
ceq1GradQx = zeros(Nt - 2, nDecVar);
ceq1GradQy = zeros(Nt - 2, nDecVar);
ceq1GradQz = zeros(Nt - 2, nDecVar);

ceq1Rx = zeros(Nt - 2, 1);
ceq1Ry = zeros(Nt - 2, 1);
ceq1Rz = zeros(Nt - 2, 1);
ceq1Qx = zeros(Nt - 2, 1);
ceq1Qy = zeros(Nt - 2, 1);
ceq1Qz = zeros(Nt - 2, 1);

for ii = 1:Nt-2
    ceq1GradRx(ii, :) = dynGradRx(2*ii, :) - dynGradRx(2*ii + 1, :);
    ceq1GradRy(ii, :) = dynGradRy(2*ii, :) - dynGradRy(2*ii + 1, :);
    ceq1GradRz(ii, :) = dynGradRz(2*ii, :) - dynGradRz(2*ii + 1, :);
    
    ceq1GradQx(ii, :) = dynGradQx(2*ii, :) - dynGradQx(2*ii + 1, :);
    ceq1GradQy(ii, :) = dynGradQy(2*ii, :) - dynGradQy(2*ii + 1, :);
    ceq1GradQz(ii, :) = dynGradQz(2*ii, :) - dynGradQz(2*ii + 1, :);
    
    ceq1Rx(ii, :) = ddrxBase(2*ii, :) - ddrxBase(2*ii + 1, :);
    ceq1Ry(ii, :) = ddryBase(2*ii, :) - ddryBase(2*ii + 1, :);
    ceq1Rz(ii, :) = ddrzBase(2*ii, :) - ddrzBase(2*ii + 1, :);
    ceq1Qx(ii, :) = ddqxBase(2*ii, :) - ddqxBase(2*ii + 1, :);
    ceq1Qy(ii, :) = ddqyBase(2*ii, :) - ddqyBase(2*ii + 1, :);
    ceq1Qz(ii, :) = ddqzBase(2*ii, :) - ddqzBase(2*ii + 1, :);
end

ceq2GradRx = zeros(2, nDecVar);
ceq2GradRx(1, :) = dynGradRx(1, :);
ceq2GradRx(2, :) = dynGradRx(end, :);

ceq2GradRy = zeros(2, nDecVar);
ceq2GradRy(1, :) = dynGradRy(1, :);
ceq2GradRy(2, :) = dynGradRy(end, :);

ceq2GradRz = zeros(2, nDecVar);
ceq2GradRz(1, :) = dynGradRz(1, :);
ceq2GradRz(2, :) = dynGradRz(end, :);

ceq2GradQx = zeros(2, nDecVar);
ceq2GradQx(1, :) = dynGradQx(1, :);
ceq2GradQx(2, :) = dynGradQx(end, :);

ceq2GradQy = zeros(2, nDecVar);
ceq2GradQy(1, :) = dynGradQy(1, :);
ceq2GradQy(2, :) = dynGradQy(end, :);

ceq2GradQz = zeros(2, nDecVar);
ceq2GradQz(1, :) = dynGradQz(1, :);
ceq2GradQz(2, :) = dynGradQz(end, :);

ceq2Rx = zeros(2, 1);
ceq2Rx(1, :) = ddrxBase(1, :);
ceq2Rx(2, :) = ddrxBase(end, :);

ceq2Ry = zeros(2, 1);
ceq2Ry(1, :) = ddryBase(1, :);
ceq2Ry(2, :) = ddryBase(end, :);

ceq2Rz = zeros(2, 1);
ceq2Rz(1, :) = ddrzBase(1, :);
ceq2Rz(2, :) = ddrzBase(end, :);

ceq2Qx = zeros(2, 1);
ceq2Qx(1, :) = ddqxBase(1, :);
ceq2Qx(2, :) = ddqxBase(end, :);

ceq2Qy = zeros(2, 1);
ceq2Qy(1, :) = ddqyBase(1, :);
ceq2Qy(2, :) = ddqyBase(end, :);

ceq2Qz = zeros(2, 1);
ceq2Qz(1, :) = ddqzBase(1, :);
ceq2Qz(2, :) = ddqzBase(end, :);

cGrad = [];
ceqGrad = [ceq1GradRx; ceq2GradRx; ceq1GradRy; ceq2GradRy; ceq1GradRz; ceq2GradRz; ...
    ceq1GradQx; ceq2GradQx; ceq1GradQy; ceq2GradQy; ceq1GradQz; ceq2GradQz];

c = [];
ceq = [ceq1Rx; ceq2Rx; ceq1Ry; ceq2Ry; ceq1Rz; ceq2Rz; ...
    ceq1Qx; ceq2Qx; ceq1Qy; ceq2Qy; ceq1Qz; ceq2Qz];

idxAcc = 1:2:2*(Nt-1);
ddx = [ddrxBase(idxAcc), ddryBase(idxAcc), ddrzBase(idxAcc), ddqxBase(idxAcc), ddqyBase(idxAcc), ddqzBase(idxAcc)]';

end