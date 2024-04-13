function base = guessBaseMotion(p, tArray)

% duration of phases
T1 = tArray(1); % stance phase
T2 = tArray(2);   T3 = tArray(3); % swing phase and stance phase
T4 = tArray(4);   T5 = tArray(5); % swing phase and stance phase

Nt = p.Nt;

T = T1 + T2 + T3 + T4 + T5;
tSpan = linspace(0, T, Nt);
dt = T / (Nt - 1);

t0 = 0;                     t1 = T1;                        t2 = T1 + T2;
t3 = T1 + T2 + T3;          t4 = T1 + T2 + T3 + T4;         t5 = T1 + T2 + T3 + T4 + T5;

r0 = p.q0(1:3);             r1 = [0.1; 0.02; 0.55];         r2 = [0.2; 0.04; 0.6];
r3 = [0.3; 0.07; 0.58];     r4 = [0.45; 0.09; 0.56];        r5 = p.qF(1:3);

q0 = p.q0(4:6);             q1 = deg2rad([5; -10; 5]);      q2 = deg2rad([-8; 8; -1]);
q3 = deg2rad([7; -9; 6]);   q4 = deg2rad([-5; 7; -2]);      q5 = p.qF(4:6);

guess.r = [r0, r1, r2, r3, r4, r5];
guess.q = [q0, q1, q2, q3, q4, q5];
guess.t = [t0, t1, t2, t3, t4, t5];
rBase = interp1(guess.t', guess.r', tSpan', 'pchip')';
qBase = interp1(guess.t', guess.q', tSpan', 'pchip')';

drBase = [zeros(3, 1), diff(rBase, 1, 2)/dt];
dqBase = [zeros(3, 1), diff(qBase, 1, 2)/dt];

drBase(:, end) = 0;
dqBase(:, end) = 0;

base = [rBase; qBase; drBase; dqBase];

end