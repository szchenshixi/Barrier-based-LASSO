% Symbols used in this project is consistent with the paper
% "An Interior-Point Method for Large-Scale 1-Regularized Least Squares"
% though with minor adjustments as required by the assignment.

% Initial settings as given in the assignment that should not be modified.
randn('seed',1); 
% rng(1,'v4')
global g_beta0;
g_beta0 = zeros(10,1); 
g_beta0(3) = 1; 
g_beta0(5) = 7; 
g_beta0(10) = 3; 
global g_n g_p;
g_n = 100; g_p = 10;
global g_X; % nxp
global g_y; % nx1
global g_lambda;  % scalar
g_X = randn(g_n,g_p); 
g_y= g_X*g_beta0 + 0.1*randn(g_n,1); 
g_lambda = 0.2;

% User defined parameters
% Line search
global g_a;
global g_b;
g_a = 0.1; % Slope of the reference line
g_b = 0.9;   % The iterative decrement coefficient

% Newton's method
global g_mu;
global g_m0; % The coefficient associated with the logarithmic barriers.
g_mu = 1.1;
g_m0 = 1; % The initial weight of the original objective function


