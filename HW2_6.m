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

errorvect = zeros(1,1); 
nu1 = 2; 
nu2 = 1; 
w = 2/3; 
n=32; 
n1 = n; 
h=1/n; 
x = 0:h:1;
%Need to have vh be length n-1


%homogeneous f(x)
b = 0;
%inhomogeneous f(x) = 1
%b = 1; 
switch b
    case 0 
        x1 = x(2:end-1); 
        v = rand(size(x1)); 
        f = zeros(length(x)-2,1);
        true = zeros(size(x)); 
        figure(b+1)
    case 1
        x1 = x(2:end-1); 
        v =rand(size(x1)); 
        f = ones(length(x)-2,1);
        true = -0.5*x.^2 + 0.5*x;
        figure(b)
end

%Initialize guess for v and rhs

v0 = [0 v 0]; 
error1 = max(abs((v0 - true))); 
%threshold = error1/10^5
threshold = 1e-5; 
counter = 1; 
plot(x,v0); 
hold on; 
while error1 > threshold
    errorvect(counter) = error1;
    [v] = vcycle(h,f,v, nu1, nu2);
%Gives the solution at the coarsest grid
h=1/n1; 
x = 0:h:1;
%true = zeros(size(x)); 
plot(x,true); 
hold on;
v = v'; 
y = [0 v 0]; 
plot(x,y,'rx')
error1 = max(abs(true-y));
counter = counter + 1;
end
error1

figure(2)
semilogy(errorvect,'m')
figure(3)
plot(x,y,'rx-', x, true); 
hold on; 