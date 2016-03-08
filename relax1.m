function [v_new,residual] = relax1(omega, v, rhs, n, h)
%Relax - Performs weighted Jacobi relaxation with a given omega value.
%Takes the v_approx, f_approx on the given grid and does n sweeps with i. 
%       Weighted Jacobi: R_j = D^(-1)(L+U), c = D^(-1)b
%           v_m+1 = [(1-w)I + wR_j]v_m + wc
%
N = length(v); 
w = omega; 
e = ones(N,1); 
D = 1/(h^2)*spdiags(2*e,0,N,N); 
L = 1/(h^2)*spdiags(e,-1,N,N); 
U = L'; 
if iscolumn(rhs) == 0
    rhs = rhs'; 
end
if iscolumn(v)==0
    v = v'; 
end

%R_j formed using the Jacobi matrix



%Weighted Jacobi matrix formed 
for i = 1:n

    v = (1-w)*v + w*(D\((L+U)*v + rhs)); 
end
v_new = v ;
residual = rhs - (D-L-U)*v_new; 


