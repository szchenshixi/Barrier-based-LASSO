globalParas;

global g_beta0;
global g_m0; % The coefficient associated with the logarithmic barriers.
global g_p;
t = 20*ones(g_p,1);
[beta, t, iter] = NewtonsMethod(g_beta0, t, g_m0);