function soln = MX_trajOptimizer(problem)

% ������Ҫ��problem���guess���д���
G = problem.guess;
F = problem.func;
Opt = problem.options;

% ��ȡ�Ż�����Decision Variables
[zGuess, pack] = packDecVar(G.time, G.motion_foot, G.force_foot, G.base);

% % ��֤һ�½⹹����ȷ��
% [t, p, f, x] = unPackDecVar(zGuess, pack);
% if ~isempty(find(t - G.time)) || ~isempty(find(p - G.motion_foot)) || ~isempty(find(f - G.force_foot)) || ~isempty(find(x - G.base))
%     error('Check your pack decision variables functions!');
% end

% �ṩ��Լ���Ľ����ݶ�
flagGradCst = strcmp(Opt.nlpOpt.GradConstr,'on');

% Ŀ�꺯��û�У����ֲ�����Ϊempty,����������ʹ�����0
P.objective = @(z)( ...
    myObjective(z) );

% ������Լ��
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

% ��ʼ���
tic; 
[zSoln, ~, exitFlag, output] = fmincon(P);
[tSoln, pSoln, fSoln, xSoln] = unPackDecVar(zSoln, pack);
nlpTime = toc;

% ������
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