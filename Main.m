globalParas;

global g_beta0;
global g_m0; % The coefficient associated with the logarithmic barriers.
t = 1*g_beta0 + 20*ones(size(g_beta0));  % Just a random guess of the initial fesible t.
[beta, t, iter, dualityGap, obj] = NewtonsMethod(g_beta0, t, g_m0);
objective(beta)
index = 1:iter;
subplot(211)
plot(index, log(dualityGap)/log(10));
title("Duality gap versus Newton's iteration");
xlabel("iteration")
ylabel("gap (in log10)")
subplot(212)
plot(index, obj);
title("Objective function value versus Newton's iteration");
xlabel("iteration")
ylabel("Objective function value")