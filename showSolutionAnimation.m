function showSolutionAnimation(soln)

global param

t = soln.t;
p = soln.p;
f = soln.f;
x = soln.x;

T = sum(t);

%% 先得到初始的各个时刻对应的状态量
pNode = getFootMotionNodes(t, p, param.Nt);
fNode = getFootForceNodes(t, f, param.Nt);

u = [pNode; fNode];
tSpan = linspace(0, T, param.Nt);

dx = dynamics(tSpan, x, u, param);

%% 对应于各个时间段多项式的更加细分时刻nMesh的状态量
nMesh = 10;
Nt1 = (param.Nt - 1) * nMesh + 1;

pNode1 = getFootMotionNodes(t, p, Nt1);
fNode1 = getFootForceNodes(t, f, Nt1);

xNode1 = getBaseMotionNodes(t, x, dx, Nt1);

tSpan1 = linspace(0, T, Nt1);

%% 将得到的结果套到单腿模型中利用逆运动学进行展示
lx = 0.3;   ly = 0.2;   lz = 0.2; % base 的长、宽、高

pLeg.j1z = 0.5 * lz; % 髋关节与Base重心的竖直高度
pLeg.l1 = 0.3; % 大腿长
pLeg.l2 = 0.3; % 小腿长

% 初始关节角
q0 = [0; 0; -0.8411; 0.8411];

% base的顶点及各个面，画图时用
base = makeRigidBody(lx, ly, lz);

% base的初始位姿
base.p = param.q0(1:3);
base.R = getRotMatFromBaseToWorld(param.q0(4:6));

figHandle = figure(1001); clf;
% figHandle.WindowState = 'fullscreen';
figHandle.Color = [1 1 1];

% 三维图各个轴的范围
AX = [-0.5, 1.0];   AY = [-0.5, 0.5];    AZ = [-0.1, 1.0];

for n = 1:Nt1
    
    rBase = xNode1(1:3, n);
    qBase = xNode1(4:6, n);
    
    R = getRotMatFromBaseToWorld(qBase);
    
    base.R = R;
    base.p = rBase;
    
    pFoot = pNode1(:, n);
    
    pFootBase = R' * (pFoot - rBase);
    
    % 数值解求逆运动学IK
    q = getJointAngle(q0, pFootBase, pLeg);
    
    % 得到髋、膝、踝的点在世界坐标系下的位置
    pPoints = getLegPointsPos(rBase,qBase,q,pLeg);
    
    pHip = pPoints(:, 1);
    pKnee = pPoints(:, 2);
    pAnkle = pPoints(:, 3);
    
    % base的顶点需要做一下姿态转换
    vert = base.R * base.vertex;
    for k = 1:3
        vert(k, :) = vert(k, :) + base.p(k);
    end
    
    hold off
    
    newplot;
    
    hh = patch('faces',base.face','vertices',vert','FaceColor',0.8 * [0.5 0.5 0.5]);
    
    hold on
    
    plot3([pHip(1) pKnee(1)], [pHip(2) pKnee(2)], [pHip(3) pKnee(3)], 'LineWidth',5, 'Color', [200,60,60]/255);
    plot3([pKnee(1) pAnkle(1)], [pKnee(2) pAnkle(2)], [pKnee(3) pAnkle(3)], 'LineWidth',5, 'Color', [60,60,200]/255);
    
    plot3(pHip(1), pHip(2), pHip(3),  'k.','MarkerSize',20,'LineWidth',2);
    plot3(pKnee(1), pKnee(2), pKnee(3),  'k.','MarkerSize',20,'LineWidth',2);
    plot3(pAnkle(1), pAnkle(2), pAnkle(3),  'k.','MarkerSize',20,'LineWidth',2);
    
    % ground
    xBnd = AX(1):0.01:AX(2);
    yBnd = AY(1):0.01:AY(2);

    [Gnd.x, Gnd.y] = meshgrid(xBnd, yBnd);
    Gnd.z = zeros(size(Gnd.x));
    GndColor = 1.3 * [115, 74, 18] / 255;

    surf(Gnd.x, Gnd.y, Gnd.z, 'EdgeColor', 'none', 'FaceColor', GndColor);
    
    % force vector
    arrowScale = 1/450;
    f1 = fNode1(:, n);
    quiver3(pFoot(1), pFoot(2), pFoot(3), f1(1), f1(2), f1(3), arrowScale, 'Color' , 'r', 'LineWidth', 2);
    
    view(25, 15); 
    
    title(sprintf('time : %6.3fs', tSpan1(n)));
    
    axis equal;  
    xlim(AX); ylim(AY); zlim(AZ);    
    xlabel('x');    ylabel('y');    zlabel('z');    
    axis off        
    drawnow
    
    q0 = q;
    
end

end

function base = makeRigidBody(w, d, h)

wdh = [w, d, h];

base.vertex = 0.5*[
-wdh(1) -wdh(2) -wdh(3);
-wdh(1) wdh(2) -wdh(3);
wdh(1) wdh(2) -wdh(3);
wdh(1) -wdh(2) -wdh(3);
-wdh(1) -wdh(2) wdh(3);
-wdh(1) wdh(2) wdh(3);
wdh(1) wdh(2) wdh(3);
wdh(1) -wdh(2) wdh(3);
]'; % vertex points
base.face = [
1 2 3 4; 2 6 7 3; 4 3 7 8; 1 5 8 4; 1 2 6 5; 5 6 7 8;
]';

end