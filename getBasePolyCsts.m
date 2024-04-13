function [c, ceq] = getBasePolyCsts(rBase, drBase, ddrBase, axis, dt)

Nt = size(rBase, 2);

idxLow = 1:(Nt - 1);
idxUpp = 2:Nt;

switch axis
    case 'x'
        idxRow = 1;
    case 'y'
        idxRow = 2;
    case 'z'
        idxRow = 3;
    otherwise
        error('Invalid axis!');
end

rx = rBase(idxRow, :);
drx = drBase(idxRow, :);
ddrx = ddrBase(idxRow, :);

rxLow = rx(:, idxLow);
rxUpp = rx(:, idxUpp);
drxLow = drx(:, idxLow);
drxUpp = drx(:, idxUpp);
ddrxLow = ddrx(:, idxLow);
% ddrxUpp = ddrx(:, idxUpp);

S = cell(Nt - 1, 1);
for ii = 1:(Nt-1)
    
    % ��i�ε�λ�á��ٶȡ����ٶ�
    x0 = rxLow(ii);
    dx0 = drxLow(ii);
    ddx0 = ddrxLow(ii);
    
    % ��i+1�ε�λ�á��ٶȡ����ٶ�
    x1 = rxUpp(ii);
    dx1 = drxUpp(ii);
%     ddx1 = ddrxUpp(ii);
    
    S{ii} = get4orderPolyCoeff(dt, x0, dx0, ddx0, x1, dx1);
    
end

% ��ʽԼ��1�� ��i��dt�Ľ���ʱ�̵ļ��ٶȵ��ڵ�ii+1��dt�Ŀ�ʼʱ�̵ļ��ٶ�
ceq1 = zeros(Nt - 2, 1);
for ii = 1:Nt-2
  
    [~, ~, ddS_dt] = get4orderPoly(S{ii}, dt);
    [~, ~, ddS_0] = get4orderPoly(S{ii+1}, 0);
    
    ceq1(ii) = ddS_dt - ddS_0;
end

% ��ʽԼ��2�� �������̿�ʼ�ͽ����ļ��ٶȵ���0
[~, ~,  ddS_0] = get4orderPoly(S{1}, 0);
[~, ~,  ddS_dt] = get4orderPoly(S{Nt - 1}, dt);
ceq2 = [ddS_0; ddS_dt];

% ����
c = [];
ceq = [ceq1; ceq2];

end