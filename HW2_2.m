%%Homework 2: Problem 2, solving -u'' = exp(x) with two sets of boundary
%%conditions
% Set A: u(0) = 1, u(1) = 0
% Set B: u'(0) = 1, u(1) = 0
% Use Jacobi Method
clear;
h=0.1; 
l=1; 
norms1 = zeros(3,3); 
n = 'B'; 

switch n
    case 'A'  
        %Set A
        while h>= 0.00625
            x = 0:h:1;
            N = length(x);
            %initial guess
            u0 = ones(N,1);
            size(u0)
            u(1) = 1;
       
            b = exp(x);
            b(ceil(end/2)) = b(ceil(end/2)) + h;
            boundary_data = zeros(N-2,1); 
            boundary_data(1) = 1; 
            %When u(1) = h
            boundary_data(end)=h; 
            b = b(2:end-1); 
            figure(1);
            hold on;
            u_true = 2+(exp(1) - 2)*x - exp(x);
            plot(x,u_true)
            main = 2*ones(N-2,1); 
            s = -1*ones(N-3,1); 
            A = (diag(main) + diag(s,-1) + diag(s,1));
            b = b';
            b=h^2*b+boundary_data;
            u = A\b;
            u = u'; 
            
            u = [1 u 0]; 
            
            plot(x,u,'r')
            error = u-u_true;

            norm1 = 0;
            %1 norm
            for i = 1:length(error)
                norm1 = norm1 + abs(error(i)) ;
            end
            norm1 = norm1*h;
            norms1(l,1) = norm1;
            
            %inf norm
            infnorm = max(abs(error));
            norms1(l,3) = infnorm;
            norm2 = 0;
            %two norm
            for i = 1:length(error)
                norm2 = norm2 + abs(error(i))^2;
            end
            norm2 = (h*norm2)^(1/2);
            norms1(l,2) = norm2;
            l = l+1;
            
            h = h/2;
            h
        end
        
        % Order of convergence
        format short
        fprintf('$h = 0.1$ & %.4f & %.4f & %.4f\\\\ \n',-log(norms1(2,1)/norms1(1,1))/log(2),-log(norms1(2,2)/norms1(1,2))/log(2),-log(norms1(2,3)/norms1(1,3))/log(2));
        fprintf('$h = 0.05$ & %.4f & %.4f & %.4f\\\\ \n',-log(norms1(3,1)/norms1(2,1))/log(2),-log(norms1(3,2)/norms1(2,2))/log(2),-log(norms1(3,3)/norms1(2,3))/log(2));
        fprintf('$h = 0.025$ & %.4f & %.4f & %.4f\\\\ \n',-log(norms1(4,1)/norms1(3,1))/log(2),-log(norms1(4,2)/norms1(3,2))/log(2),-log(norms1(4,3)/norms1(3,3))/log(2));
        fprintf('$h = 0.0125$ & %.4f & %.4f & %.4f\\\\ \n',-log(norms1(5,1)/norms1(4,1))/log(2),-log(norms1(5,2)/norms1(4,2))/log(2),-log(norms1(5,3)/norms1(4,3))/log(2));
        
        p_1 = -log(norms1(2,1)/norms1(1,1))/log(2);
        p_1 = -log(norms1(3,1)/norms1(2,1))/log(2);
        p_1 = -log(norms1(4,1)/norms1(3,1))/log(2);
        p_1 = -log(norms1(5,1)/norms1(4,1))/log(2);
        p_2 = -log(norms1(2,2)/norms1(1,2))/log(2);
        p_2 = -log(norms1(3,2)/norms1(2,2))/log(2);
        p_2 = -log(norms1(4,1)/norms1(3,1))/log(2);
        p_2 = -log(norms1(5,1)/norms1(4,1))/log(2);
        p_i = -log(norms1(2,3)/norms1(1,3))/log(2);
        p_i = -log(norms1(3,3)/norms1(2,3))/log(2);
        p_i = -log(norms1(4,1)/norms1(3,1))/log(2);
        p_i = -log(norms1(5,1)/norms1(4,1))/log(2);
    
    case 'B'
        %Set B: Discretize u'(0) = 1 two ways: (u1-u0)/h = 1 or (?3u0 + 4u1 ? u2)/(2h) = 1
        while h>= 0.00625
            x = 0:h:1;
            N = length(x);
            %initial guess
            u0 = ones(N,1);
            size(u0)
            u(1) = 1;
       
            b = exp(x);
            boundary_data = zeros(N-2,1); 
            %For first order one sided
            %boundary_data(1) = -h; 
            %For first order one sided edited
            %boundary_data(1) = -h - h^2; 
            %for second order one sided
            %boundary_data(1) = -2/3*h;
            %For second order one sided edited
            boundary_data(1) = -2/3*(h+h^2);
            b = b(2:end-1); 
            figure(1);
            hold on;
            u_true = 2*x+exp(1)-2 - exp(x);
            plot(x,u_true)
            main = 2*ones(N-2,1); 
            s = -1*ones(N-3,1); 
            A = (diag(main) + diag(s,-1) + diag(s,1));
            %%% With one sided first order, 
            %A(1,1) = 1; 
            %%% With one sided seoncd order; 
            A(1,1) = 2-4/3; 
            A(1,2) = -(1-1/3); 
            
            b = b';
            b=h^2*b+boundary_data;
            u = A\b;
            u = u'; 
            %One sided first order
            a = u(1)-h; 
            u = [a u 0]; 
            
            plot(x,u,'r')
            error = u-u_true;

            norm1 = 0;
            %1 norm
            for i = 1:length(error)
                norm1 = norm1 + abs(error(i)) ;
            end
            norm1 = norm1*h;
            norms1(l,1) = norm1;
            
            %inf norm
            infnorm = max(abs(error));
            norms1(l,3) = infnorm;
            norm2 = 0;
            %two norm
            for i = 1:length(error)
                norm2 = norm2 + abs(error(i))^2;
            end
            norm2 = (h*norm2)^(1/2);
            norms1(l,2) = norm2;
            l = l+1;
            
            h = h/2;
            h
        end
        
        % Order of convergence
        format short
        fprintf('$h = 0.1$ & %.4f & %.4f & %.4f\\\\ \n',-log(norms1(2,1)/norms1(1,1))/log(2),-log(norms1(2,2)/norms1(1,2))/log(2),-log(norms1(2,3)/norms1(1,3))/log(2));
        fprintf('$h = 0.05$ & %.4f & %.4f & %.4f\\\\ \n',-log(norms1(3,1)/norms1(2,1))/log(2),-log(norms1(3,2)/norms1(2,2))/log(2),-log(norms1(3,3)/norms1(2,3))/log(2));
        fprintf('$h = 0.025$ & %.4f & %.4f & %.4f\\\\ \n',-log(norms1(4,1)/norms1(3,1))/log(2),-log(norms1(4,2)/norms1(3,2))/log(2),-log(norms1(4,3)/norms1(3,3))/log(2));
        fprintf('$h = 0.0125$ & %.4f & %.4f & %.4f\\\\ \n',-log(norms1(5,1)/norms1(4,1))/log(2),-log(norms1(5,2)/norms1(4,2))/log(2),-log(norms1(5,3)/norms1(4,3))/log(2));
        
        p_1 = -log(norms1(2,1)/norms1(1,1))/log(2);
        p_1 = -log(norms1(3,1)/norms1(2,1))/log(2);
        p_1 = -log(norms1(4,1)/norms1(3,1))/log(2);
        p_1 = -log(norms1(5,1)/norms1(4,1))/log(2);
        p_2 = -log(norms1(2,2)/norms1(1,2))/log(2);
        p_2 = -log(norms1(3,2)/norms1(2,2))/log(2);
        p_2 = -log(norms1(4,1)/norms1(3,1))/log(2);
        p_2 = -log(norms1(5,1)/norms1(4,1))/log(2);
        p_i = -log(norms1(2,3)/norms1(1,3))/log(2);
        p_i = -log(norms1(3,3)/norms1(2,3))/log(2);
        p_i = -log(norms1(4,1)/norms1(3,1))/log(2);
        p_i = -log(norms1(5,1)/norms1(4,1))/log(2);


end