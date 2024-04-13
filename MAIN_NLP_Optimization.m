%% Program Profile
% ���������Ҫ���ܣ�����ֱ�ӷ��Ż������嶯��ѧ���˶�x������ʱ��T���ŵ��˶�p�Լ��Ӵ���f�ı仯
% ����ѧģ�ͣ������嶯��ѧģ��(6ά�˶���x, y, z, roll, pitch, yaw)
% ����ʱ�䣺֧���ࡢ������ʱ��
% �ŵ��˶������ݲ�ͬ��ѡȡ��ͬ��������������
% �Ӵ����仯�����ݲ�ͬ��ѡȡ��ͬ��������������(��ŵ��˶�����)

% ��������Ҫ�ο�A. Winkler��2018����RAL��������ģ���Gait and Trajectory Optimization for
% Legged Systems through Phase-based End-Effector Parameterization��

% �������ṩ����Լ������ֵ�ݶȣ�����ٶȸ��졣

% Written by Meng Xiang in BIT on 2021/10/18
% Copyrights Reserved.

%% Main program
clear; close all; clc

global param

% invariant parameters
param = getModelParameters();

% ��ʼ�ͽ�����λ�ú���̬: q = [x; y; z; qx; qy; qz]�����ѿ�������ϵ
% ŷ����˳��ZYX 
param.q0 = [  0;   0; 0.5; 0; 0; 0]; 
param.qF = [0.5; 0.1; 0.5; 0; 0; 0];

%% ��ʼ�²�ֵ
UsePreviousResults = 'false'; % 'true' or 'false'
% �ȸ�һ��guess��guess���������ȫ���Ż���������Ҫ���Ż�֮ǰ��guess������ݽ��д����ֶ�����������M.
% Posa(IJRR)����ĽӴ�����Ӵ����˶������Ի���Լ��(LCP),�������ʣ�������ϵ�ҽ���̽�֡�
% �����Ҷ������˶��ĳ�ʼ�²�ֵ
% phase duration: T1 T2 T3 T4 T5
problem.guess.time = 1.0 * [0.45, 0.2, 0.3, 0.15, 0.4];
% foot motion and force
problem.guess.motion_foot = guessFootMotion();
problem.guess.force_foot = guessFootForce(param);
% q and dq
problem.guess.base = guessBaseMotion(param, problem.guess.time);

if UsePreviousResults
% ���ǻ���������guess�õ����Ż�������������Զ��Ż�������п��ӻ�
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

%% �Ӻ��������綯��ѧ�ĺ���,�߽�Լ����·��Լ���ȣ��Լ��������������ѡ������

addpath(genpath('./AnalyticalGradientDerive'));
addpath ./EulerConvertZYX
addpath ./KinematicsModel
addpath(genpath('./MX_Polynomials345')); % �����ļ��м������ļ���ȫ�����뵽·��

% sub-functions
problem.func.dynamics = @(t,x,u)( dynamics(t, x, u, param) );

problem.func.bndCst = @(x0, xF)( boundaryStateCst(x0, xF, param) );

problem.func.pathCst = @(t,x,u)( pathConstraint(t, x, u, param) );

% fmincon��ѡ�����ã�
problem.options.nlpOpt = optimset(...
    'Display','iter',...   % {'iter','final','off'}
    'TolFun',1e-6,...
    'MaxFunEvals',5e4,...
    'MaxIter', 600);   %options for fmincon

problem.options.nlpOpt.GradConstr = 'on';
problem.options.nlpOpt.GradObj = 'on';
problem.options.nlpOpt.DerivativeCheck = 'off';

% �����Ҫ������ṩ�Ľ����ݶ���fminconʹ�õ���ֵ�ݶ�֮���Ƿ���ڽϴ������齫��FinDiffType�����óɡ�central��
% �����central���������Ȼ��ʧ�ܣ��Ǻ��п���˵�����ṩ�Ľ����ݶ���һ�������⣬������м��
% problem.options.nlpOpt = optimset(...
%     'Display','iter',...   % {'iter','final','off'}
%     'TolFun',1e-6,...
%     'MaxFunEvals',5e4,...
%     'FinDiffType', 'central');   %options for fmincon
% problem.options.nlpOpt.DerivativeCheck = 'on';

%% ����켣�Ż�����transcription,ת����NLP������������ķ���������
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