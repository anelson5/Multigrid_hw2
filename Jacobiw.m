function [x,residual] = Jacobiw(omega,x0,b,iter,h)
% [x,noit,err] = Jacobi(A,b,x0,eps,maxit);
%
% Purpose: Solve Ax=b using Jacobi iteration, with initial guess x0
%           to tolerance eps in less that maxit iterations.
%
% Output : x    Solution
%          noit Number of iterations
%          err  L2 of difference between last two iterations
%          Q    Splitting matrix

N = length(b);
w = omega;  
d = 2*ones(size(x0));
N = length(x0); 
s = -1*ones(length(x0),1); 
A = spdiags([s d s], -1:1, N,N); 
A = 1/(h^2)*A; 
% Create splitting matrix
Q = zeros(N);
for i=1:N
  Q(i,i) = A(i,i);
end;

% Compute reduction factor
R = inv(Q)*(Q-omega*A);
rho = max(abs(eig(R)));

% Set initial values
x= x0;
noit = iter;
err = 1.0;

while (noit > 0)
    y = (Q-omega*A)*x + omega*b;
    xn = Q\y;
    err = norm(xn-x);
    errvec(noit+1) = err;
    noit = noit-1
    x = xn;

end;

residual = A*x-b; 
return
