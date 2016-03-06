%%Script for Homework 2 problem 6. V Cycle warm up problem for the problem
%%-u''(x) = f(x), 0<x<1, u(0)=0, u(1)=0. 
%
%   Procedure: At each level, going down
%       (1) Apply weighted Jacobi to current approximation (v_1) with
%       initial guess 0 (EXCEPT on first level v^k)
%       (2) compute f = I_restrict*u
%   Procedure: At coarsest grid
%       (3) Solve exactly A_coarsest u = = f
%   Procedure: At each level, going up
%       (4) Correct v_level = v_level + I_interpolate*r_previouslevel
%       (5) Apply weighted Jacobi on level with initial guess v_level
%   
%   Since we're looking at a one dimensional problem by the Galerkin
%   property, we know that A_coarser = I_restrict*A_finer*I_interpolate

w = 2/3; 
n=8; 
h = 1/n
x = 0:h:1; 
size(x);

main = 2*ones(n-1,1); 
off = -1*ones(n-1,1); 
Ah = spdiags([off main off],-1:1, n-1,n-1);
Ah = 1/h^2*Ah; 

levels = 4; 

for i = 1:levels
    if i ==1 
        v = x.^2;
        v = v(2:end-1); 
        v_array{i} = v; 
        f_array{i} = zeros(length(x)-2,1); 
    else
        h = h*2;
        x = 0:h:1; 
    v_array{i} = zeros(length(x)-2,1); 
    f_array{i} = zeros(length(x)-2,1); 
    end
end
v_array = fliplr(v_array)
f_array = fliplr(f_array);
%Gives the solution at the coarsest grid
[v_array,f_array,Ah_coarsest] = coarsen(Ah,f_array,v_array,levels); 
[v_array, f_array] = refine(Ah,f_array,v_array,2);