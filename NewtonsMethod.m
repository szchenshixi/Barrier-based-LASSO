% This function solves the original problem by using interior point method
% and introduces a series of logarithmic barriers to substitudes the
% inequality constraints. m is the parameter associated with the
% logarithmic barriers.
function [beta, t, iter] = NewtonsMethod(beta, t, m0)
global g_mu;
global g_p;
m = m0;
iter = 0;
while(dualityGap(beta) > 1e-3)
    iter = iter + 1;
    step = stepComput(beta, t, m);
    beta = beta + step(1:g_p);
    t = t + step(g_p+1,end);
    m = g_mu * m;
end
end

function grad = gradientComput(beta, t, m)
global g_X;
global g_y;
global g_lambda;
global g_p;
temp1 = 2*beta./(t.^2-beta.^2);
gradientBeta = 2*m*g_X'*(g_X*beta-g_y) + temp1;
temp2 = -2*t./(t.^2-beta.^2);
gradientT = m*g_lambda*ones(g_p,1) + temp2;
grad = [gradientBeta; gradientT];
end

function Hessian = hessianComput(beta, t, m)
global g_X;
D1 = diag(2*(t.^2 + beta.^2)./(t.^2-beta.^2));
D2 = diag(-4*t.*beta./(t.^2-beta.^2));
Hessian = [2*m*g_X'*g_X + D1 D2;
           D2 D1];
end

function step = stepComput(beta, t, m)
global g_a;
global g_b;
global g_p;
gradient = gradientComput(beta, t, m);
Hessian = hessianComput(beta, t, m);
direction = -Hessian^-1*gradient;
step = direction;   % Initial step is equal to the direction vector
while(objective(beta + step(1:g_p)) > objective(beta) + g_a*gradient(1:g_p)'*step(1:g_p))
    step = g_b * step;
end
end

function f0 = objective(beta)
global g_X;   % nxp
global g_y;   % nx1
global g_lambda;  %scalar
firstTerm = (g_y-g_X*beta)' * (g_y-g_X*beta);
secondTerm = g_lambda * ones(size(beta))' * abs(beta);
f0 = firstTerm + secondTerm;
end