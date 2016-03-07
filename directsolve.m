function [v_new] = directsolve(v,rhs,h)
%Direct solve Summary of this function goes here
%   Detailed explanation goes here
d = 2*ones(size(v));
N = length(v); 
s = -1*ones(length(v),1); 
A = spdiags([s d s], -1:1, N,N); 
A = 1/(h^2)*A; 

v_new = A\rhs; 


end

