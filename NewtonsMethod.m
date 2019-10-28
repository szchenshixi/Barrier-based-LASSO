% This function solves the original problem by using interior point method
% and introduces a series of logarithmic barriers to substitudes the
% inequality constraints. m is the parameter associated with the
% logarithmic barriers.
function [beta, t, iter, gap, obj] = NewtonsMethod(beta, t, m0)
global g_mu;
global g_p;
global g_epsilon;
m = m0;
iter = 0;
exitCheck = inf;
gap = [];
obj = [];
while(exitCheck/2 > g_epsilon)
    iter = iter + 1;
    [step, exitCheck] = stepComput(beta, t, m);
    beta = beta + step(1:g_p);
    t = t + step(g_p+1:end);
    m = m * g_mu;
    gap = [gap dualityGap(beta)];
    obj = [obj objective(beta)];
end
end

function grad = gradientComput(beta, t, m)
global g_X;
global g_y;
global g_lambda;
global g_p;
temp1 = 2*beta./(t.^2-beta.^2);
gradientBeta = 2*m*(g_X')*(g_X*beta-g_y) + temp1;
temp2 = -2*t./(t.^2-beta.^2);
gradientT = m*g_lambda*ones(g_p,1) + temp2;
grad = [gradientBeta; gradientT];
end

function Hessian = hessianComput(beta, t, m)
global g_X;
D1_num = 2*(t.^2 + beta.^2);
D1_deno = (t.^2-beta.^2).^2;
D1 = diag(D1_num./D1_deno);
D2_num = -4*t.*beta;
D2_deno = (t.^2-beta.^2).^2;
D2 = diag(D2_num./D2_deno);
Hessian = [2*m*(g_X')*g_X + D1 D2;
           D2 D1];
end

function [step, exitCheck] = stepComput(beta, t, m)
global g_a;
global g_b;
global g_p;
gradient = gradientComput(beta, t, m);
Hessian = hessianComput(beta, t, m);
direction = -inv(Hessian)*gradient;
step = 2 * direction;   % Initial step constructed from the direction vector
exitCheck = (1/m^3) * (gradient') * inv(Hessian) * gradient;   % The square of Newton decrement
% counter = 0;
while(augmentedObjective(beta + step(1 : g_p),                      ...
                         t + step(g_p+1 : 2*g_p),                   ...
                         m) > augmentedObjective(beta, t, m) +      ...
                         g_a*(gradient')*step                       ...
                         || min(t + step(g_p+1 : 2*g_p)) < 0        ...
     )
    step = g_b * step;
%     counter = counter + 1;
%     if(counter > 1e2)
%         break;
%     end
end
end