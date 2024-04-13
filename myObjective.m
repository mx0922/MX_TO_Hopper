function [cost, costGrad] = myObjective(z)

cost = 0.0;

nDecVar = numel(z);

costGrad = zeros(1, nDecVar);

end