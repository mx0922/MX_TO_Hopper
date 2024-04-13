function soln = MX_trajOptimizer(problem)

% 首先需要对problem里的guess进行处理
G = problem.guess;
F = problem.func;
Opt = problem.options;

% 提取优化变量Decision Variables
[zGuess, pack] = packDecVar(G.time, G.motion_foot, G.force_foot, G.base);

% % 验证一下解构的正确性
% [t, p, f, x] = unPackDecVar(zGuess, pack);
% if ~isempty(find(t - G.time)) || ~isempty(find(p - G.motion_foot)) || ~isempty(find(f - G.force_foot)) || ~isempty(find(x - G.base))
%     error('Check your pack decision variables functions!');
% end

% 提供了约束的解析梯度
flagGradCst = strcmp(Opt.nlpOpt.GradConstr,'on');

% 目标函数没有，但又不能设为empty,因此在这里仅使其等于0
P.objective = @(z)( ...
    myObjective(z) );

% 非线性约束
if flagGradCst
    P.nonlcon = @(z)(...
        myConstraintGrad(z, pack, F.dynamics, F.bndCst, F.pathCst) );
    [~,~,cstIneqInit,cstEqInit] = P.nonlcon(zGuess);
    sparsityPattern.equalityConstraint = (cstEqInit~=0)';  % Only used for visualization! 
    sparsityPattern.inequalityConstraint = (cstIneqInit~=0)';  % Only used for visualization! 
else
    P.nonlcon = @(z)(...
        myConstraint(z, pack, F.dynamics, F.bndCst, F.pathCst) );    
end

P.x0 = zGuess;
P.lb = [];
P.ub = [];
P.Aineq = []; P.bineq = [];
P.Aeq = []; P.beq = [];
P.options = Opt.nlpOpt;
P.solver = 'fmincon';

% 开始求解
tic; 
[zSoln, ~, exitFlag, output] = fmincon(P);
[tSoln, pSoln, fSoln, xSoln] = unPackDecVar(zSoln, pack);
nlpTime = toc;

% 输出结果
soln.t = tSoln;
soln.p = pSoln;
soln.f = fSoln;
soln.x = xSoln;
soln.nlpTime = nlpTime;

soln.info.exitFlag = exitFlag;
soln.info.output = output;

if flagGradCst
    [~,~,cstIneqInit,cstEqInit] = P.nonlcon(zSoln);
    sparsityPattern.equalityConstraint = (cstEqInit~=0)';
    sparsityPattern.inequalityConstraint = (cstIneqInit~=0)';
end
soln.info.sparsityPattern = sparsityPattern;

end