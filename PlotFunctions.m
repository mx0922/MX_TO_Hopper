if 0
figure(1); clf;
pz = pNode1(3, :);
fz = fNode1(3, :);

xlabel('time [ s ]')

% 左纵轴，脚的位置
yyaxis left
plot(tSpan1, pz)
hold on
% 0参考线
line([0, T], [0, 0], 'linestyle', '--');
y1Bnd = [-0.1, 0.6];
ylim(y1Bnd);
% ylabel的名称及将名称旋转
y1 = ylabel('\it P^{z}');
set(y1, 'Rotation', 0, 'FontSize', 20)

% 右纵轴，接触力
yyaxis right
plot(tSpan1, fz)
hold on
% 0参考线
line([0, T], [0, 0], 'linestyle', '--');
y2Bnd = [-50, 350];
ylim(y2Bnd);
% ylabel的名称及将名称旋转
y2 = ylabel('\it F^{z}');
set(y2, 'Rotation', 0, 'FontSize', 20)

% grid on
T1 = t(1);
T2 = sum(t(1:2));
T3 = sum(t(1:3));
T4 = sum(t(1:4));
T5 = T;
% x轴按照自己想要显示的点的位置进行划分，每个点的值只取小数点后两位数
xticks([T1, T2, T3, T4, T5])
xticklabels([roundn(T1, -2), roundn(T2,-2), roundn(T3,-2), roundn(T4,-2), roundn(T5,-2)])

% xA = 0:0.01:T1;
% yA = 350 * ones(size(xA));
% area(xA, yA, 'FaceColor', [0.5, 0.5, 0.5]);

xA1 = [0, T1, T1, 0];
yA1 = [y2Bnd(1), y2Bnd(1), y2Bnd(2), y2Bnd(2)];
% fill(xA, yA, 'FaceColor', [0.5, 0.5, 0.5]); % 会报错，很奇怪
face1 = fill(xA1, yA1, 'b');
set(face1, 'FaceColor', [0.5, 0.5, 0.5], 'EdgeColor', 'none', 'FaceAlpha', 0.2);

xA2 = [T2, T3, T3, T2];
yA2 = [y2Bnd(1), y2Bnd(1), y2Bnd(2), y2Bnd(2)];
face2 = fill(xA2, yA2, 'b');
set(face2, 'FaceColor', [0.5, 0.5, 0.5], 'EdgeColor', 'none', 'FaceAlpha', 0.2);

xA3 = [T4, T5, T5, T4];
yA3 = [y2Bnd(1), y2Bnd(1), y2Bnd(2), y2Bnd(2)];
face3 = fill(xA3, yA3, 'b');
set(face3, 'FaceColor', [0.5, 0.5, 0.5], 'EdgeColor', 'none', 'FaceAlpha', 0.2);

% title('p^{z} & f^{z}')
end

figHandle = figure(2); clf;

figHandle.WindowState = 'fullscreen';
% figHandle.WindowStyle = 'docked';
figHandle.Color = [1 1 1];

zBase = xNode1(3, :);

% ax1 = subplot(2, 1, 1); clf;
tiledlayout(2, 1)

ax1 = nexttile;

plot(tSpan1, zBase, 'b', 'LineWidth', 3)
hold on
line([0, T], [0, 0], 'linestyle', '--');

for ii = 3:2:(length(tSpan)-2)
    line([tSpan(ii), tSpan(ii)], yy1Bnd, 'linestyle', '-', 'color', 1.1 * [0.8, 0.8, 0.8]);    
end

plot(tSpan, x(3, :), 'r.','MarkerSize',20,'LineWidth',2);

yy1Bnd = [-0.1, 0.9];

ylim(yy1Bnd);
xlim([0, T]);
% ax1 = gca;
set(gca, 'XTick', tSpan(1:2:end), 'XTickLabel', roundn(tSpan(1:2:end), -2))

grid on

xlabel('time[ s ]')
yy = ylabel('\it R^{z}');
set(yy, 'Rotation', 0, 'FontSize', 20);

%% ax2 = subplot(2, 1, 2); clf; 
ax2 = nexttile;

pz = pNode1(3, :);
fz = fNode1(3, :);

xlabel('time [ s ]');

% ax2 = gca;

% 左纵轴，脚的位置
yyaxis left
plot(tSpan1, pz)
hold on
% 0参考线
line([0, T], [0, 0], 'linestyle', '--');
y1Bnd = [-0.1, 0.6];
ylim(y1Bnd);
% ylabel的名称及将名称旋转
y1 = ylabel('\it P^{z}');
set(y1, 'Rotation', 0, 'FontSize', 20)

% 右纵轴，接触力
yyaxis right
plot(tSpan1, fz)
hold on
% 0参考线
line([0, T], [0, 0], 'linestyle', '--');
y2Bnd = [-50, 350];
ylim(y2Bnd);
% ylabel的名称及将名称旋转
y2 = ylabel('\it F^{z}');
set(y2, 'Rotation', 0, 'FontSize', 20)

% grid on
T1 = t(1);
T2 = sum(t(1:2));
T3 = sum(t(1:3));
T4 = sum(t(1:4));
T5 = T;
% x轴按照自己想要显示的点的位置进行划分，每个点的值只取小数点后两位数
xticks([T1, T2, T3, T4, T5])
xticklabels([roundn(T1, -2), roundn(T2,-2), roundn(T3,-2), roundn(T4,-2), roundn(T5,-2)])

xA1 = [0, T1, T1, 0];
yA1 = [y2Bnd(1), y2Bnd(1), y2Bnd(2), y2Bnd(2)];
% fill(xA, yA, 'FaceColor', [0.5, 0.5, 0.5]); % 会报错，很奇怪
face1 = fill(xA1, yA1, 'b');
set(face1, 'FaceColor', [0.5, 0.5, 0.5], 'EdgeColor', 'none', 'FaceAlpha', 0.2);

xA2 = [T2, T3, T3, T2];
yA2 = [y2Bnd(1), y2Bnd(1), y2Bnd(2), y2Bnd(2)];
face2 = fill(xA2, yA2, 'b');
set(face2, 'FaceColor', [0.5, 0.5, 0.5], 'EdgeColor', 'none', 'FaceAlpha', 0.2);

xA3 = [T4, T5, T5, T4];
yA3 = [y2Bnd(1), y2Bnd(1), y2Bnd(2), y2Bnd(2)];
face3 = fill(xA3, yA3, 'b');
set(face3, 'FaceColor', [0.5, 0.5, 0.5], 'EdgeColor', 'none', 'FaceAlpha', 0.2);

linkaxes([ax1, ax2], 'x');