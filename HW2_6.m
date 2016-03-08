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
nu1 = 1; 
nu2 = 1; 
w = 2/3; 
n=128; 
h=1/n; 
x = 0:h:1;
%Need to have vh be length n-1

levels = 5; 

%Initialize guess for v and rhs
x1 = x(2:end-1); 
v =ones(size(x1)); 
f = ones(length(x)-2,1);
%v = ones(length(x)-2,1); 
[v] = vcycle(h,f,v,levels, nu1, nu2);

%Gives the solution at the coarsest grid
n=128; 
h=1/n; 
x = 0:h:1;
true = -0.5*x.^2 + 0.5*x;
%true = zeros(size(x)); 
plot(x,true); 
hold on;
v = v'; 
size(v)
size(true)
y = [0 v 0]; 
size(y)
plot(x,y,'rx')
error = norm(true-y)