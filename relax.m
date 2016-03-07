function [v_new,residual] = relax(omega, v, rhs, n, h)
%Relax - Performs weighted Jacobi relaxation with a given omega value.
%Takes the v_approx, f_approx on the given grid and does n sweeps with i. 
%       Weighted Jacobi: R_j = D^(-1)(L+U), c = D^(-1)b
%           v_m+1 = [(1-w)I + wR_j]v_m + wc
%
w = omega;  
d = 2*ones(size(v));
N = length(v); 
s = -1*ones(length(v),1); 
A = spdiags([s d s], -1:1, N,N); 
A = 1/(h^2)*A; 
%Extract the lower, upper and diagonal part of A, to form R_j
L = tril(A); 
U = triu(A); 
D = diag(diag(A)); 

if iscolumn(rhs) == 0
    rhs = rhs'; 
end
if iscolumn(v)==0
    v = v'; 
end

%R_j formed using the Jacobi matrix
R_j = D\(L+U); 

c = D\rhs; 
I = eye(size(R_j));

%Weighted Jacobi matrix formed 
R_w = (1-w)*I + w*R_j;

for i = 1:n

    v = R_w*v + w*c; 
end
v_new = v 
residual = rhs - A*v_new; 

