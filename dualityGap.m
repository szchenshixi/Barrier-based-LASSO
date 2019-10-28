function gap = dualityGap(beta)
global g_X;   % nxp
global g_y;   % nx1
global g_lambda;  %scalar
s = min(g_lambda./(abs(2*g_X'*(g_X*beta-g_y))));
v = 2*s*(g_X*beta-g_y);
gap = objective(beta) - dualObjective(v);
end