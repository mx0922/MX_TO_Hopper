%% Program Profile
% 本程序的主要功能：利用直接法优化单刚体动力学的运动x、迈步时间T、脚的运动p以及接触力f的变化
% 动力学模型：单刚体动力学模型(6维运动：x, y, z, roll, pitch, yaw)
% 迈步时间：支撑相、飞行相时间
% 脚的运动：根据不同相选取不同个数的三次曲线
% 接触力变化：根据不同相选取不同个数的三次曲线(与脚的运动类似)

% 本程序主要参考A. Winkler于2018年在RAL发表的论文：《Gait and Trajectory Optimization for
% Legged Systems through Phase-based End-Effector Parameterization》

% 本程序提供所有约束的数值梯度，求解速度更快。

% Written by Meng Xiang in BIT on 2021/10/18
% Copyrights Reserved.

%% Main program
clear; close all; clc

global param

% invariant parameters
param = getModelParameters();

% 初始和结束的位置和姿态: q = [x; y; z; qx; qy; qz]——笛卡尔坐标系
% 欧拉角顺序：ZYX 
param.q0 = [  0;   0; 0.5; 0; 0; 0]; 
param.qF = [0.5; 0.1; 0.5; 0; 0; 0];

%% 初始猜测值
UsePreviousResults = 'false'; % 'true' or 'false'
% 先给一个guess，guess里的量并不全是优化变量，需要在优化之前对guess里的数据进行处理，手动满足类似于M.
% Posa(IJRR)提出的接触力与接触点运动的线性互补约束(LCP),如有疑问，可以联系我进行探讨。
% 这是我对整个运动的初始猜测值
% phase duration: T1 T2 T3 T4 T5
problem.guess.time = 1.0 * [0.45, 0.2, 0.3, 0.15, 0.4];
% foot motion and force
problem.guess.motion_foot = guessFootMotion();
problem.guess.force_foot = guessFootForce(param);
% q and dq
problem.guess.base = guessBaseMotion(param, problem.guess.time);

if UsePreviousResults
% 这是基于上述的guess得到的优化结果，后续可以对优化结果进行可视化
% cd ./OptimResults
% 
% load('tOptim.mat');
% problem.guess.time = tOptim;
% % foot motion and force
% load('pOptim.mat');
% load('fOptim.mat');
% problem.guess.motion_foot = pOptim;
% problem.guess.force_foot = fOptim;
% % q and dq
% load('xOptim.mat');
% problem.guess.base = xOptim;
% 
% cd ..
end

%% 子函数，例如动力学的函数,边界约束，路径约束等，以及非线性求解器的选项设置

addpath(genpath('./AnalyticalGradientDerive'));
addpath ./EulerConvertZYX
addpath ./KinematicsModel
addpath(genpath('./MX_Polynomials345')); % 将该文件夹及其子文件夹全都加入到路径

% sub-functions
problem.func.dynamics = @(t,x,u)( dynamics(t, x, u, param) );

problem.func.bndCst = @(x0, xF)( boundaryStateCst(x0, xF, param) );

problem.func.pathCst = @(t,x,u)( pathConstraint(t, x, u, param) );

% fmincon的选项设置：
problem.options.nlpOpt = optimset(...
    'Display','iter',...   % {'iter','final','off'}
    'TolFun',1e-6,...
    'MaxFunEvals',5e4,...
    'MaxIter', 600);   %options for fmincon

problem.options.nlpOpt.GradConstr = 'on';
problem.options.nlpOpt.GradObj = 'on';
problem.options.nlpOpt.DerivativeCheck = 'off';

% 如果想要检查你提供的解析梯度与fmincon使用的数值梯度之间是否存在较大误差，建议将‘FinDiffType’设置成‘central’
% 如果‘central’情况下仍然会失败，那很有可能说明你提供的解析梯度有一定的问题，建议进行检查
% problem.options.nlpOpt = optimset(...
%     'Display','iter',...   % {'iter','final','off'}
%     'TolFun',1e-6,...
%     'MaxFunEvals',5e4,...
%     'FinDiffType', 'central');   %options for fmincon
% problem.options.nlpOpt.DerivativeCheck = 'on';

%% 进入轨迹优化器做transcription,转换成NLP求解器可以求解的非线性问题
soln = MX_trajOptimizer(problem);

if isfield(soln.info,'sparsityPattern')
   figure(101); clf;
   spy(soln.info.sparsityPattern.equalityConstraint);
   axis equal
   title('Sparsity pattern in equality constraints');
   
   figure(102); clf;
   spy(soln.info.sparsityPattern.inequalityConstraint);
   axis equal
   title('Sparsity pattern in inequality constraints');
end

showSolutionAnimation(soln);
% AnimationVideo(soln);