function p = getModelParameters()

% mass and gravity
p.m = 20.0;
p.g = 9.81;
p.I = diag([0.25, 0.18, 0.16]);

% time 
p.T = 1.5;
p.minT = 0.1;

p.Nt = 31;

% kinematic model
p.p_hat = [0; 0; -0.5];
p.L_cube = 0.2;

% friction coefficient
p.miu = 0.6;

end