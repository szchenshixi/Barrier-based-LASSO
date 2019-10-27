globalParas;

global g_beta0;
global g_m0; % The coefficient associated with the logarithmic barriers.
t = 2*g_beta0 + ones(size(g_beta0));
[beta, t, iter] = NewtonsMethod(g_beta0, t, g_m0);