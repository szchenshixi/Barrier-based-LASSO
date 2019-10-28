randn('seed',1); 
n=50; 
P=randn(3*n,n); 
M=P'*P; 
x0_c = 2/sqrt(n)*ones(n, 1);  % Current x0
x0_n = 2/sqrt(n)*ones(n, 1);  % Next x0
% x0_c = 1/n*ones(n, 1);  % Current x0
% x0_n = 1/n*ones(n, 1);  % Next x0
m = 1;  % The coefficient associated with the logarithmic barriers.

counter = 0;
minEigOpt = [];
minEigReal = min(eig(M));
while(1)
    cvx_begin
    variables x0_n(n);
    obj = (m*(x0_n)'*M*x0_n) - log(2*(x0_c)'*(x0_n-x0_c) + (x0_c)'*x0_c - 1);
    minimize obj;
    cvx_end
    
    m = m * 5;
    x0_c = x0_n;
    minEigOpt_ = (x0_c)'*M*(x0_c);
    if(abs(minEigOpt_ - minEigReal) < 1e-3)
        break;
    end
    minEigOpt = [minEigOpt minEigOpt_];
end
index = 1:length(minEigOpt);
plot(index, minEigOpt);
title("Minimal eigen value by cvx optimization")
xlabel("iteration")
ylabel("approximated min eigen value")