function [v_array, f_array,Ah_fine] = refine(Ah,f_array,v_array,l)
%UNTITLED12 Summary of this function goes here
%   Detailed explanation goes here
levels = l; 
w = 2/3; 
N = 2; %relaxation parameter, how many times we relax
if levels == 4
    f_array_temp = f_array; 
    v_array_temp = v_array; 
    f = f_array{levels};
    v = v_array{levels};
    v_array_temp{levels} = relax(w,v, f,N,Ah);
    Ah_fine_temp = Ah; 
else 
    % Form interpolation matrix and restriction matrix to make A^(2h) =
    % I_h^(2h)A^hI_(2h)^h
    v = v_array{levels};
    f = f_array{levels};
    n = length(v); % is length n/2-1
    d = 2*ones(n,1);
    d1 = ones(n+1,1);
    
    %create matrix with twos on diagonal size n/2-1 by n/2-1
    two_matrix = spdiags(d, 0, n+1,n);
    %create matrix with twos on diagonal size n/2 by n/2-1
    ones_matrix = spdiags([d1 d1],-1:0, n+1,n);
    %ones_matrix = full(ones_matrix)
    rows = size(two_matrix,1) + size(ones_matrix,1); %adding the number of rows of the two matrices
    %two_matrix = full(two_matrix)
    row_interweave = reshape([ones_matrix(:) two_matrix(:) ]',rows, []);
    row_interweave = row_interweave(1:end-1,:);
    
    I_restrict = row_interweave; 
    I_interp = row_interweave'; 
    Ah_fine = I_interp*Ah*I_restrict; 
    

    v_array{levels} = v_array{levels} + I_interp*v_array{levels-1}; 
    v_new = relax(w,v_array{levels}, f,N,Ah_fine);
    f_array{levels} = r_new; 
    v_array{levels} = v_new; 
    [v_array_temp,f_array_temp,Ah_fine_temp] = refine(Ah_fine, v_array, f_array, levels+1);        
end
 v_array = v_array_temp; 
 f_array = f_array_temp; 
 Ah_fine = Ah_fine_temp; 

end

