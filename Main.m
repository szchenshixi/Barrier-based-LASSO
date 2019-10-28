globalParas;

global g_beta0;
global g_m0; % The coefficient associated with the logarithmic barriers.
t = 1.5*g_beta0 + 0.3*ones(size(g_beta0));
[beta, t, iter, dualityGap, obj] = NewtonsMethod(g_beta0, t, g_m0);
objective(beta)