function [v_array,f_array,Ah] = coarsen(Ah,f_array,v_array,l)
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here
% l - level numbers 
% f - right hand side
% Ah - discretized operator 
levels = l; 
w = 2/3; 
N = 2; %relaxation parameter, how many times we relax
if levels == 1 
    %direct solve step 
    v_array_temp = v_array; 
    f_array_temp = f_array; 
    f = f_array{levels};
    if iscolumn(f) ==0
        f = f'; 
    end
    Ah_temp = Ah; 
    v = Ah_temp\f;
    v_array_temp{levels} = v; 
    f_array_temp{levels} =f; 
else 
    % Form interpolation matrix and restriction matrix to make A^(2h) =
    % I_h^(2h)A^hI_(2h)^h
    v = v_array{levels};
    f = f_array{levels};
    if iscolumn(f) ==0
        f = f';
    end
    if iscolumn(v)==0 
        v = v'; 
    end
    levels
    n = length(v) % is length n-1
    N = (n-1)/2; 
    d = 2*ones(N,1);
    d1 = ones(N+1,1);
    
    %create matrix with twos on diagonal size n/2-1 by n/2-1
    two_matrix = spdiags(d, 0, N+1,N);
    %create matrix with oness on diagonal size n/2 by n/2-1
    ones_matrix = spdiags([d1 d1],-1:0, N+1,N);
   % ones_matrix = full(ones_matrix);
    rows = size(two_matrix,1) + size(ones_matrix,1); %adding the number of rows of the two matrices
   % two_matrix = full(two_matrix);
    row_interweave = reshape([ones_matrix(:) two_matrix(:) ]',rows, []);
    row_interweave = row_interweave(1:end-1,:);
    I_interp = row_interweave; 
    I_restrict = row_interweave'; 
    levels
    size(I_restrict)
    size(I_interp)
    size(Ah)
    Ah_coarse = I_restrict*Ah*I_interp; 
    
    %relax v_1 times with initial guess
    if levels == l
        %Use different initial guess for wGS
        v_new = relax(w, v, f,N,Ah);
      
        r = f - Ah*v_new; 
        e_new = v - v_new; 
    else 
        v_new = relax(w, zeros(size(v)), f,N,Ah);
        r = f - Ah*v_new; 
       % e_new = v - v_new; 
    end
    v_array{levels} = v_new; 
    [r_new] = restrict('fw',r); 
    f_array{levels} = r_new; 
    fprintf('size of matrices Ah, v,f  before the recursive call \n')
    size(Ah_coarse)
    size(v_array)
    size(f_array)
    [v_array_temp,f_array_temp,Ah_temp] = coarsen(Ah_coarse, v_array, f_array, levels-1); 
end
v_array = v_array_temp; 
f_array = f_array_temp;
Ah = Ah_temp;
end

