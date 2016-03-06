function [v_new] = relax(omega, approx, rhs, n, A )
%Relax - Performs weighted Jacobi relaxation with a given omega value.
%Takes the v_approx, f_approx on the given grid and does n sweeps with i. 
%       Weighted Jacobi: R_j = D^(-1)(L+U), c = D^(-1)b
%           v_m+1 = [(1-w)I + wR_j]v_m + wc
%
w = omega; 
v = approx; 
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
v_new = v; 

