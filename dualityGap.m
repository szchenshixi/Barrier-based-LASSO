function gap = dualityGap(beta)
global g_X;   % nxp
global g_y;   % nx1
global g_lambda;  %scalar
tempDeno = abs(2*g_X'*g_y - 2*g_X'*g_X*beta);
s = min(g_lambda./tempDeno);
v = 2*s*(g_y-g_X*beta);
gap = objective(beta) - dualObjective(v);
end

function G = dualObjective(v)
global g_y;   % nx1
G = -(1/4)*v'*v - v'*g_y;
end

function f0 = objective(beta)
global g_X;   % nxp
global g_y;   % nx1
global g_lambda;  %scalar
firstTerm = (g_y-g_X*beta)' * (g_y-g_X*beta);
secondTerm = g_lambda * ones(size(beta))' * abs(beta);
f0 = firstTerm + secondTerm;
end