function [v] = vcycle(h,f,v,l, nu1, nu2)
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here
% l - level numbers 
% h - step size
% f_array - cell array of rhs
% v_arry - cell array of v
% nui - relaxation parameter. 
levels = l; 
w = 2/3; 
levels

if levels == 1 
    %direct solve
    v = directsolve(v,f,h); 
    %go to next relax
    
else
    [v_relax, residual] = relax1(w, v, f, nu1, h);
    f_coarse = restrict('fw',residual); 
    vcoarse = zeros(size(f_coarse)); 
    coarser = 2*h; 
    [v_cycle] = vcycle(coarser,f_coarse, vcoarse,levels-1, nu1, nu2); 
    v_new = v_relax+ interpolate(v_cycle);
    [v, residual] = relax1(w, v_new, f, nu1, h);
 
end

