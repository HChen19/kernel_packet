%% compute basis function phi(x) and coefficient band matrix A
function [A, phi, phi_fun] = compute_basis_solver_new(x, k, x_phi, rho)
%==========================================================================
%initialization
n = size(x,2);%number of observations
temp_index = 1:size(x_phi,1);%temp_index = [1,2,3,...,size(x,1)]

row_A = zeros((n-k+1)*k,1);
col_A = zeros((n-k+1)*k,1);
values_A = zeros((n-k+1)*k,1);

row_phi = zeros((n-k+1)*k,1);                                             
col_phi = zeros((n-k+1)*k,1);
values_phi = zeros((n-k+1)*k,1);

if length(x_phi) <= length(x)
    row_phi_fun = zeros((n-k+1)*k,1);
    col_phi_fun = zeros((n-k+1)*k,1);
    values_phi_fun = zeros((n-k+1)*k,1);
else
    row_phi_fun = zeros((n-k+1)*(length(x_phi)-k+1),1);
    col_phi_fun = zeros((n-k+1)*(length(x_phi)-k+1),1);
    values_phi_fun = zeros((n-k+1)*(length(x_phi)-k+1),1);
end

ii = 1;
jj = 1;
%==========================================================================



%==========================================================================
%intermediate basis phi_[(k+1)/2 : (n-(k-1)/2)]
for i = (k+1)/2 : (n-(k-1)/2)
    %%compute V_temp
    V_temp = [x(i-(k-1)/2 : i+(k-1)/2) - (x(i-(k-1)/2)+x(i-(k-1)/2))/2] ...
        .^([0:(k-3)/2]');
    
    %%compute A
    null_temp =  ...
        [V_temp ...
        .*exp((x(i-(k-1)/2 : i+(k-1)/2) - x(i-(k-1)/2)).*sqrt(k-2)./rho); ...
        V_temp ...
        .*exp((-x(i-(k-1)/2 : i+(k-1)/2) + x(i-(k-1)/2)).*sqrt(k-2)./rho); ...
        [1,zeros(1,(k-1))]] ...                                           
        \[zeros((k-1),1);1];
    
    
    row_A(ii : ii+k-1) = i-(k-1)/2 : i+(k-1)/2;
    col_A(ii : ii+k-1) = ones(k,1)*i;
    values_A(ii : ii+k-1) = null_temp;
    
    %compute phi
    %distance = abs(x(i-(k-1)/2 : i+(k-1)/2)' - x(i-(k-1)/2 : i+(k-1)/2));
    x_mid = x(i-(k-1)/2 : i+(k-1)/2);
    phi_temp = matern_halfint(x_mid', x_mid', (k-2)/2, 1, rho) ...
        *null_temp;  
                                                   
    row_phi(ii : ii+k-1) = i-(k-1)/2 : i+(k-1)/2;
    col_phi(ii : ii+k-1) = ones(k,1)*i;
    values_phi(ii : ii+k-1) = phi_temp;
    
    ii = ii + k;
    
    %%compute phi_fun
    x_logic = [x_phi <= x(i+(k-1)/2) & x_phi >= x(i-(k-1)/2)];
    if isempty(x_phi(x_logic))
        continue
    end
    %distance_phi_fun = abs(x_phi(x_logic) - x(i-(k-1)/2 : i+(k-1)/2));
    %distance = [|x-x1|,|x-x2|,...,|x-xk|]
    x_index = temp_index(x_logic);
    phi_fun_temp = matern_halfint(x_phi(x_logic), x_mid', (k-2)/2, 1, rho) ...
        *null_temp;
    
    row_phi_fun(jj : jj+length(x_index)-1) = x_index;
    col_phi_fun(jj : jj+length(x_index)-1) = ones(length(x_index),1)*i;
    values_phi_fun(jj : jj+length(x_index)-1) = phi_fun_temp;
    

    jj = jj + length(x_index);
end%end for loop over intermediate basis

A = sparse(row_A, col_A, values_A, n, n);
phi = sparse(row_phi, col_phi, values_phi, n, n);
phi_fun = sparse(row_phi_fun(1:jj-1), col_phi_fun(1:jj-1), ...
    values_phi_fun(1:jj-1), length(x_phi), n);
%==========================================================================



%==========================================================================
%boundary basis phi_[1 : (k-1)/2] & x[n-(k-3)/2 : n]
for i = 1 : (k-1)/2
    %%compute V_temp_left & right
    V_temp_left = [x(1 : i+(k-1)/2) - (x(1)+x(i+(k-1)/2))/2] ...
        .^([0:(k-3)/2]');
    V_temp_right = [x(n-k+1+i : n) - (x(n-k+1+i)+x(n))/2] ...
        .^([0:(k-3)/2]');
    
    %%compute A
    if i == 1
        %left-most A[:,1]
        A(1 : i+(k-1)/2,i) = [V_temp_left ...
            .*exp((x(1 : i+(k-1)/2) - x(1)).*sqrt(k-2)./rho); ...
            [1,zeros(1,i+(k-3)/2)]] ...
            \[zeros(i+(k-3)/2,1);1];
        %right boundary A[:,n-(k-1)/2+1]
        A(n-k+1+i : n, n-(k-1)/2+i) = ...
            [V_temp_right(1 : (k-1)/2-i, :) ...
            .*exp((x(n-k+1+i : n) - x(n-k+1+i)).*sqrt(k-2)./rho); ...
            V_temp_right ...
            .*exp((-x(n-k+1+i : n) + x(n-k+1+i)).*sqrt(k-2)./rho); ...
            [1,zeros(1,k-1-i)]] ...
            \[zeros(k-1-i,1);1];
    elseif i == (k-1)/2
        %right-most basis A[:,n]
        A(n-k+1+i : n, n-(k-1)/2+i) = [V_temp_right ...
            .*exp((-x(n-k+1+i : n) + x(n-k+1+i)).*sqrt(k-2)./rho); ...
            [1,zeros(1,k-1-i)]] ...
            \[zeros(k-1-i,1);1];
        %left boundary A[:,(k-1)/2]
        A(1 : i+(k-1)/2,i) = [V_temp_left ...
            .*exp((x(1 : i+(k-1)/2) - x(1)).*sqrt(k-2)./rho); ...
            V_temp_left(1 : i-1, :) ...
            .*exp((-x(1 : i+(k-1)/2) + x(1)).*sqrt(k-2)./rho); ...
            [1,zeros(1,i+(k-3)/2)]] ...
            \[zeros(i+(k-3)/2,1);1];
    else
        %left boundary A[:,i]
        A(1 : i+(k-1)/2,i) = [V_temp_left ...
            .*exp((x(1 : i+(k-1)/2) - x(1)).*sqrt(k-2)./rho); ...
            V_temp_left(1 : i-1, :) ...
            .*exp((-x(1 : i+(k-1)/2) + x(1)).*sqrt(k-2)./rho); ...
            [1,zeros(1,i+(k-3)/2)]] ...
            \[zeros(i+(k-3)/2,1);1];
        %right boundary A[:,n-(k-1)/2+i]
        A(n-k+1+i : n, n-(k-1)/2+i) = ...
            [V_temp_right(1 : (k-1)/2-i, :) ...
            .*exp((x(n-k+1+i : n) - x(n-k+1+i)).*sqrt(k-2)./rho); ...
            V_temp_right ...
            .*exp((-x(n-k+1+i : n) + x(n-k+1+i)).*sqrt(k-2)./rho); ...
            [1,zeros(1,k-1-i)]] ...
            \[zeros(k-1-i,1);1];
    end%end if statement when computing A in boundary basis
    
    %%compute phi                                                         
    %left boundary basis                                                  
    %distance_left = abs(x(1 : i+(k-1)/2)' - x(1 : i+(k-1)/2));
    x_left = x(1 : i+(k-1)/2);
    x_index_left = [1 : i+(k-1)/2];                                       
    phi(x_index_left, i) = ...                                            
        matern_halfint(x_left', x_left', (k-2)/2, 1, rho) ...          
        *A(1 : i+(k-1)/2, i);                                             
                                                                          
    %right boundary basis                                                 
    %distance_right = abs(x(n-k+1+i : n)' - x(n-k+1+i : n));  
    x_right = x(n-k+1+i : n);
    x_index_right = [n-k+1+i : n];                                        
    phi(x_index_right, n-(k-1)/2+i) = ...                                 
        matern_halfint(x_right', x_right', (k-2)/2, 1, rho) ...         
        *A(n-k+1+i : n, n-(k-1)/2+i);
    
    %%compute phi_fun
    %left boundary basis
    x_logic_left = [x_phi <= x(i+(k-1)/2)];
    if ~isempty(x_phi(x_logic_left))
        %distance_left_phi_fun = abs(x_phi(x_logic_left) - x(1 : i+(k-1)/2));
        x_index_left_phi_fun = temp_index(x_logic_left);
        phi_fun(x_index_left_phi_fun, i) = ...
            matern_halfint(x_phi(x_logic_left), x_left', (k-2)/2, 1, rho) ...
            *A(1 : i+(k-1)/2, i);
    end
    
    %right boundary basis
    x_logic_right = [x_phi >= x(n-k+1+i)];
    if ~isempty(x_phi(x_logic_right))
        %distance_right_phi_fun = abs(x_phi(x_logic_right) - x(n-k+1+i : n));
        x_index_right_phi_fun = temp_index(x_logic_right);
        phi_fun(x_index_right_phi_fun, n-(k-1)/2+i) = ...
            matern_halfint(x_phi(x_logic_right), x_right', (k-2)/2, 1, rho) ...
            *A(n-k+1+i : n, n-(k-1)/2+i);
    end
end%end for loop over boundary basis
%==========================================================================

end%end function