function xNode = getBaseMotionNodes(tPhase, x, dx, Nt)

T = sum(tPhase);

nState = size(x, 1);
nTime = size(x, 2);

tSpan = linspace(0, T, nTime);
dt = T / (nTime - 1);

tspan = linspace(0, T, Nt);

[~, bin] = histc(tspan, [-inf, tSpan, inf]);
bin = bin - 1;

xNode = zeros(nState/2 * 3, Nt);

for i = 1:(nTime - 1)
    
    xLow = x(:, i);
    xUpp = x(:, i+1);
    dxLow = dx(:, i);
    
    S_x = get4orderPolyCoeff(dt, xLow(1), xLow(7), dxLow(7), xUpp(1), xUpp(7));
    S_y = get4orderPolyCoeff(dt, xLow(2), xLow(8), dxLow(8), xUpp(2), xUpp(8));
    S_z = get4orderPolyCoeff(dt, xLow(3), xLow(9), dxLow(9), xUpp(3), xUpp(9));
    
    S_qx = get4orderPolyCoeff(dt, xLow(4), xLow(10), dxLow(10), xUpp(4), xUpp(10));
    S_qy = get4orderPolyCoeff(dt, xLow(5), xLow(11), dxLow(11), xUpp(5), xUpp(11));
    S_qz = get4orderPolyCoeff(dt, xLow(6), xLow(12), dxLow(12), xUpp(6), xUpp(12));
    
    idx = i==bin;
    if sum(idx) > 0
                
        delta = tspan(idx) - tSpan(i);
        
        [xNode(1, idx), xNode(7, idx), xNode(13, idx)] = get4orderPoly(S_x, delta);
        [xNode(2, idx), xNode(8, idx), xNode(14, idx)] = get4orderPoly(S_y, delta);
        [xNode(3, idx), xNode(9, idx), xNode(15, idx)] = get4orderPoly(S_z, delta);
        
        [xNode(4, idx), xNode(10, idx), xNode(16, idx)] = get4orderPoly(S_qx, delta);
        [xNode(5, idx), xNode(11, idx), xNode(17, idx)] = get4orderPoly(S_qy, delta);
        [xNode(6, idx), xNode(12, idx), xNode(18, idx)] = get4orderPoly(S_qz, delta);
    end
    
end

% Replace any out-of-bounds queries with NaN
outOfBounds = bin==0 | bin==(nTime+1);
xNode(:,outOfBounds) = nan;

% Check for any points that are exactly on the upper grid point:
if sum(tspan==tSpan(end))>0
    xNode(:,tspan==tSpan(end)) = [x(:,end); dx(7:end, end)];
end

end