function f0 = augmentedObjective(beta, t, m)
global g_X;   % nxp
global g_y;   % nx1
global g_lambda;  %scalar
firstTerm = m*(g_y-g_X*beta)' * (g_y-g_X*beta) + m*g_lambda*sum(t);
secondTerm = -sum(log(t-beta)) - sum(log(t+beta));
f0 = firstTerm + secondTerm;
end