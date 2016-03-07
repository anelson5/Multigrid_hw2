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

clear
clc
nu1 = 2; 
nu2 = 2; 
w = 2/3; 
n=32; 
h=1/n; 
x = 0:h:1;
%Need to have vh be length n-1
main = 2*ones(n-2,1); 
off = -1*ones(n-2,1); 
Ah = spdiags([off main off],-1:1, n-2,n-2);
Ah = 1/h^2*Ah; 

levels = 2; 

for i = 1:levels
    if i ==1
        x = x(2:end-1); 
        v= -x.^2;
        if iscolumn(v) == 0
            v_array{i} = v'; 
        end
        
        f_array{i} = ones(length(x),1); 
    else
        h = h*2;
        x = 0:h:1; 
    v_array{i} = zeros(length(x)-2,1); 
    f_array{i} = ones(length(x)-2,1); 
    end
end
f_array = fliplr(f_array);
v_array = fliplr(v_array);

[v_array,f_array] = vcycle(h,f_array,v_array,levels, nu1, nu2);

%Gives the solution at the coarsest grid
n=32; 
h=1/n; 
x = 0:h:1;
true = -0.5*x.^2 + 0.5*x; 
plot(x,true); 
hold on; 
this = v_array{levels};
this = this'; 
y = [0 this 0]; 
plot(x,y,'rx')