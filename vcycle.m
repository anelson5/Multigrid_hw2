function [v_array,f_array] = vcycle(h,f_array,v_array,l, nu1, nu2)
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
f = f_array{levels};
v = v_array{levels};

[v_new, residual] = relax(2/3, v, f, nu1, h);
v = v_new - v; 
if levels == 1 
    %direct solve
    v = directsolve(v,f,h); 
    %go to next relax
    
else
    f_coarse = restrict('fw',residual); 
    f_array{levels-1} = f_coarse; 
    [v_array,f_array] = vcycle(h/2,f_array, v_array,levels-1, nu1, nu2); 
    f_array
    v_array
end

if levels ~= 1
    vold = v_array{levels-1}; 
    v = v + interpolate(v_array{levels-1}); 
end
[v_new, residual] = relax(2/3, v, f, nu1, h);
v= v_new - v; 
v_array{levels} = v; 
 
end

