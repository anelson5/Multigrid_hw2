function [v] = vcycle(h,f,v, nu1, nu2)
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here
% l - level numbers 
% h - step size
% f_array - cell array of rhs
% v_arry - cell array of v
% nui - relaxation parameter. 
w = 2/3; 



[v_new, residual] = relax(2/3, v, f, nu1, h);
v = v_new; 
if length(v) ==1
%if levels == 1 
    %direct solve
    
    v = directsolve(v,f,h);  
    %go to next relax
    
else
    f_coarse = restrict('fw',residual); 
    vcoarse = zeros(size(f_coarse)); 
    [v_cycle] = vcycle(2*h,f_coarse, vcoarse, nu1, nu2); 
end

if length(v) ~= 1
    v = v + interpolate(v_cycle); 
end
[v_new, residual] = relax(2/3, v, f, nu1, h);
v= v_new; 
 
end

