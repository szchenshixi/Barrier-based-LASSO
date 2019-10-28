function f0 = objective(beta)
global g_X;   % nxp
global g_y;   % nx1
global g_lambda;  %scalar
firstTerm = (g_y-g_X*beta)' * (g_y-g_X*beta);
secondTerm = g_lambda * sum(abs(beta));
f0 = firstTerm + secondTerm;
end