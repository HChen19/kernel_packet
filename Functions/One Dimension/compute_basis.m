%% Compute Basis in 1D
%compute basis function phi(x) and coefficient band matrix A
function [A, varargout] = compute_basis(x, k, method, rho, x_phi)

%==========================================================================
%==========================================================================
%input: x is 1D observations, a 1*n row vector
%       k is the smoothness of Matern kernel, denotes(k-2)/2 Matern
%
%optional input: default method = 'null',default rho = sqrt(k-2),default x_phi = x'
%       method is a string = 'null' or 'solver', the way of computing A
%       rho is lengthscale parameter in Matern kernel
%       x_phi is a column vector variable of of function phi_i(x)
%
%output: A is a coefficient band matrix by solving the linear systems
%        phi = [phi_1(x),...,phi_n(x)], is n*n matrix
%
%optional output:
%        phi_fun=[phi_1(x_phi),...,phi_n(x_phi)], is n_var*n matrix
%
%where sum_{j=1}^{k} x[j]^l*exp{delta*x[j]}*A[j,i] = 0,
%      delta=+-1, l=0,...,(k-3)/2
%      n is the number of observations, n_var is the dimension of x_phi
%      phi_i(x_phi) is a column vector with the same dimension of x_phi
%      phi_i(x_phi) = sum_{j=1}^{k}K(x_phi,x_j)*A[j,i]
%                   = [K(x,x_1),...,K(x,x_k)]*A[1:k,i]
%==========================================================================
%==========================================================================

%check if input sapcing is too small
if min(diff(x)) < 1.2500e-05
    error('Spacing between points in input x is too small, try a larger spacing input.');
    return
end

%check if input include x_phi
if nargin < 5%~exist('x_phi','var')
    % third parameter does not exist, so default it to something
    x_phi = x';
end
%check if input include rho
if nargin < 4%~exist('rho','var')
    % third parameter does not exist, so default it to something
    rho = sqrt(k-2);
end
%check if input include method
if nargin < 3%~exist('method','var')
    % third parameter does not exist, so default it to something
    method = 'null';
end
if nargin < 2
    error('Not enough input arguments.')
end

if nargout < 2 || nargout > 3
    error('Incorrect output arguments.');
end

%sort the elements in x and x_phi
x = sort(x);%sort the elements of x in ascending order
x_phi = sort(x_phi);%sort the elements of x_phi in ascending order

%check if x_phi' == x
if isequal(x_phi',x) && contains(method, 'null')
    [A, varargout{1}] = compute_basis_null(x, k, rho);
elseif isequal(x_phi',x) && contains(method, 'solver')
    [A, varargout{1}] = compute_basis_solver(x, k, rho);
elseif ~isequal(x_phi',x) && contains(method, 'null')
    [A, varargout{1:2}] = compute_basis_null_new(x, k, x_phi, rho);
else
    [A, varargout{1:2}] = compute_basis_solver_new(x, k, x_phi, rho);
end%end if statement

end%end function


%%
%compute basis function phi(x) and coefficient band matrix A
function [A, phi, phi_fun] = compute_basis_null_new(x, k, x_phi, rho)
%==========================================================================
%initialization
n = size(x,2);%number of observations
temp_index = [1:size(x_phi,1)];%temp_index = [1,2,3,...,size(x,1)]

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
        null([V_temp ...
        .*exp(x(i-(k-1)/2 : i+(k-1)/2) - x(i-(k-1)/2)); ...
        V_temp ...
        .*exp(-x(i-(k-1)/2 : i+(k-1)/2) + x(i-(k-1)/2))]);
    
    row_A(ii : ii+k-1) = i-(k-1)/2 : i+(k-1)/2;
    col_A(ii : ii+k-1) = ones(k,1)*i;
    values_A(ii : ii+k-1) = null_temp;
    
    %compute phi
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
        A(1 : i+(k-1)/2,i) = null(V_temp_left ...
            .*exp(x(1 : i+(k-1)/2) - x(1)));
        %right boundary A[:,n-(k-1)/2+1]
        A(n-k+1+i : n, n-(k-1)/2+i) = ...
            null([V_temp_right(1 : (k-1)/2-i, :) ...
            .*exp(x(n-k+1+i : n) - x(n-k+1+i)); ...
            V_temp_right ...
            .*exp(-x(n-k+1+i : n) + x(n-k+1+i))]);
    elseif i == (k-1)/2
        %right-most basis A[:,n]
        A(n-k+1+i : n, n-(k-1)/2+i) = null(V_temp_right ...
            .*exp(-x(n-k+1+i : n) + x(n-k+1+i)));
        %left boundary A[:,(k-1)/2]
        A(1 : i+(k-1)/2,i) = null([V_temp_left ...
            .*exp(x(1 : i+(k-1)/2) - x(1)); ...
            V_temp_left(1 : i-1, :) ...
            .*exp(-x(1 : i+(k-1)/2) + x(1))]);
    else
        %left boundary A[:,i]
        A(1 : i+(k-1)/2,i) = null([V_temp_left ...
            .*exp(x(1 : i+(k-1)/2) - x(1)); ...
            V_temp_left(1 : i-1, :) ...
            .*exp(-x(1 : i+(k-1)/2) + x(1))]);
        %right boundary A[:,n-(k-1)/2+i]
        A(n-k+1+i : n, n-(k-1)/2+i) = ...
            null([V_temp_right(1 : (k-1)/2-i, :) ...
            .*exp(x(n-k+1+i : n) - x(n-k+1+i)); ...
            V_temp_right ...
            .*exp(-x(n-k+1+i : n) + x(n-k+1+i))]);
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

%%
%compute basis function phi(x) and coefficient band matrix A
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
        .*exp(x(i-(k-1)/2 : i+(k-1)/2) - x(i-(k-1)/2)); ...
        V_temp ...
        .*exp(-x(i-(k-1)/2 : i+(k-1)/2) + x(i-(k-1)/2)); ...
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
            .*exp(x(1 : i+(k-1)/2) - x(1)); ...
            [1,zeros(1,i+(k-3)/2)]] ...
            \[zeros(i+(k-3)/2,1);1];
        %right boundary A[:,n-(k-1)/2+1]
        A(n-k+1+i : n, n-(k-1)/2+i) = ...
            [V_temp_right(1 : (k-1)/2-i, :) ...
            .*exp(x(n-k+1+i : n) - x(n-k+1+i)); ...
            V_temp_right ...
            .*exp(-x(n-k+1+i : n) + x(n-k+1+i)); ...
            [1,zeros(1,k-1-i)]] ...
            \[zeros(k-1-i,1);1];
    elseif i == (k-1)/2
        %right-most basis A[:,n]
        A(n-k+1+i : n, n-(k-1)/2+i) = [V_temp_right ...
            .*exp(-x(n-k+1+i : n) + x(n-k+1+i)); ...
            [1,zeros(1,k-1-i)]] ...
            \[zeros(k-1-i,1);1];
        %left boundary A[:,(k-1)/2]
        A(1 : i+(k-1)/2,i) = [V_temp_left ...
            .*exp(x(1 : i+(k-1)/2) - x(1)); ...
            V_temp_left(1 : i-1, :) ...
            .*exp(-x(1 : i+(k-1)/2) + x(1)); ...
            [1,zeros(1,i+(k-3)/2)]] ...
            \[zeros(i+(k-3)/2,1);1];
    else
        %left boundary A[:,i]
        A(1 : i+(k-1)/2,i) = [V_temp_left ...
            .*exp(x(1 : i+(k-1)/2) - x(1)); ...
            V_temp_left(1 : i-1, :) ...
            .*exp(-x(1 : i+(k-1)/2) + x(1)); ...
            [1,zeros(1,i+(k-3)/2)]] ...
            \[zeros(i+(k-3)/2,1);1];
        %right boundary A[:,n-(k-1)/2+i]
        A(n-k+1+i : n, n-(k-1)/2+i) = ...
            [V_temp_right(1 : (k-1)/2-i, :) ...
            .*exp(x(n-k+1+i : n) - x(n-k+1+i)); ...
            V_temp_right ...
            .*exp(-x(n-k+1+i : n) + x(n-k+1+i)); ...
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
%%
%compute basis function phi(x) and coefficient band matrix A
function [A, phi] = compute_basis_null(x, k, rho)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%input: x is 1D observations, a row vector                                %
%       k is the smoothness of Matern kernel, denotes(k-2)/2 Matern       %
%                                                                         %
%output: A is a coefficient band matrix by solving the linear systems     %
%        phi_fun=[phi_1(x'),...,phi_n(x')], is n*n matrix                 %
%                                                                         %
%where sum_{j=1}^{k} x[j]^l*exp{delta*x[j]}*A[j,i] = 0,                   %
%      delta=+-1, l=0,...,(k-3)/2                                         %
%      n is the number of observations                                    %
%      phi_i(x') is a column vector with the same dimension of x'         %
%      phi_i(x') = sum_{j=1}^{k}K(x',x_j)*A[j,i]                          %
%                   = [K(x',x_1),...,K(x',x_k)]*A[1:k,i]                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%initialization                                                           %
n = length(x);%number of observations                                     %
                                                                          %
row_A = zeros((n-k+1)*k,1);                                               %
col_A = zeros((n-k+1)*k,1);                                               %
values_A = zeros((n-k+1)*k,1);                                            %
                                                                          %
row_phi = zeros((n-k+1)*k,1);                                             %
col_phi = zeros((n-k+1)*k,1);                                             %
values_phi = zeros((n-k+1)*k,1);                                          %
                                                                          %
ii = 1;                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%intermediate basis phi_[(k+1)/2 : (n-(k-1)/2)]                           %
for i = (k+1)/2 : (n-(k-1)/2)                                             %
    %%compute V_temp                                                      %
    V_temp = [x(i-(k-1)/2 : i+(k-1)/2) - ...                              %
        (x(i-(k-1)/2)+x(i-(k-1)/2))/2] ...                                %
        .^([0:(k-3)/2]');                                                 %
                                                                          %    
    %compute A                                                            %      
    null_temp = null([V_temp ...                                          %
        .*exp(x(i-(k-1)/2 : i+(k-1)/2) - x(i-(k-1)/2)); ...               %
        V_temp ...                                                        %
        .*exp(-x(i-(k-1)/2 : i+(k-1)/2) + x(i-(k-1)/2))]);                %
                                                                          %
    row_A(ii : ii+k-1) = i-(k-1)/2 : i+(k-1)/2;                           %
    col_A(ii : ii+k-1) = ones(k,1)*i;                                     %
    values_A(ii : ii+k-1) = null_temp;                                    %
                                                                          %
    %compute phi                                                          %
    x_mid = x(i-(k-1)/2 : i+(k-1)/2);                                     %
    phi_temp = matern_halfint(x_mid', x_mid', (k-2)/2, 1, rho) ...        %
        *null_temp;                                                       %  
                                                                          %
    row_phi(ii : ii+k-1) = i-(k-1)/2 : i+(k-1)/2;                         %
    col_phi(ii : ii+k-1) = ones(k,1)*i;                                   %
    values_phi(ii : ii+k-1) = phi_temp;                                   %
                                                                          %   
    ii = ii + k;                                                          %
end%end for loop over intermediate basis                                  %
                                                                          %
A = sparse(row_A, col_A, values_A, n, n);                                 %
phi = sparse(row_phi, col_phi, values_phi, n, n);                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%boundary basis phi_[1 : (k-1)/2] & phi_[n-(k-3)/2 : n]                   %
for i = 1 : (k-1)/2                                                       %
    %%compute V_temp_left & right                                         %
    V_temp_left = [x(1 : i+(k-1)/2) - (x(1)+x(i+(k-1)/2))/2] ...          %
        .^([0:(k-3)/2]');                                                 %
    V_temp_right = [x(n-k+1+i : n) - (x(n-k+1+i)+x(n))/2] ...             %
        .^([0:(k-3)/2]');                                                 %
                                                                          %
    %%compute A                                                           %
    if i == 1                                                             %
        %left-most A[:,1]                                                 %
        A(1 : i+(k-1)/2,i) = null(V_temp_left ...                         %
            .*exp(x(1 : i+(k-1)/2) - x(1)));                              %
        %right boundary A[:,n-(k-1)/2+1]                                  %
        A(n-k+1+i : n, n-(k-1)/2+i) = ...                                 % 
            null([V_temp_right(1 : (k-1)/2-i, :) ...                      %
            .*exp(x(n-k+1+i : n) - x(n-k+1+i)); ...                       %
            V_temp_right ...                                              %
            .*exp(-x(n-k+1+i : n) + x(n-k+1+i))]);                        %
    elseif i == (k-1)/2                                                   %
        %right-most basis A[:,n]                                          %
        A(n-k+1+i : n, n-(k-1)/2+i) = null(V_temp_right ...               %
            .*exp(-x(n-k+1+i : n) + x(n-k+1+i)));                         %
        %left boundary A[:,(k-1)/2]                                       %
        A(1 : i+(k-1)/2,i) = null([V_temp_left ...                        %
            .*exp(x(1 : i+(k-1)/2) - x(1)); ...                           %
            V_temp_left(1 : i-1, :) ...                                   %
            .*exp(-x(1 : i+(k-1)/2) + x(1))]);                            %
    else                                                                  %
        %left boundary A[:,i]                                             %
        A(1 : i+(k-1)/2,i) = null([V_temp_left ...                        %
            .*exp(x(1 : i+(k-1)/2) - x(1)); ...                           %
            V_temp_left(1 : i-1, :) ...                                   %
            .*exp(-x(1 : i+(k-1)/2) + x(1))]);                            %
        %right boundary A[:,n-(k-1)/2+i]                                  %
        A(n-k+1+i : n, n-(k-1)/2+i) = ...                                 %
            null([V_temp_right(1 : (k-1)/2-i, :) ...                      %
            .*exp(x(n-k+1+i : n) - x(n-k+1+i)); ...                       %
            V_temp_right ...                                              %
            .*exp(-x(n-k+1+i : n) + x(n-k+1+i))]);                        %
    end%end if statement                                                  %
                                                                          %
    %%compute phi                                                         %
    %left boundary basis                                                  %           %
    x_index_left = 1 : i+(k-1)/2;                                         %
    x_left = x(x_index_left);                                             %
    phi(x_index_left, i) = ...                                            %
        matern_halfint(x_left', x_left', (k-2)/2, 1, rho) ...             %
        *A(1 : i+(k-1)/2, i);                                             %
                                                                          %
    %right boundary basis                                                 %
    x_index_right = n-k+1+i : n;                                          %
    x_right = x(x_index_right);                                           %
    phi(x_index_right, n-(k-1)/2+i) = ...                                 %
        matern_halfint(x_right', x_right', (k-2)/2, 1, rho) ...           %
        *A(n-k+1+i : n, n-(k-1)/2+i);                                     %
end%end for loop over boundary basis                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end%end function
%%
%compute basis function phi(x) and coefficient band matrix A
function [A, phi] = compute_basis_solver(x, k, rho)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%input: x is 1D observations, a row vector                                %
%       k is the smoothness of Matern kernel, denotes(k-2)/2 Matern       %
%                                                                         %
%output: A is a coefficient band matrix by solving the linear systems     %
%        phi_fun=[phi_1(x'),...,phi_n(x')], is n*n matrix                 %
%                                                                         %
%where sum_{j=1}^{k} x[j]^l*exp{delta*x[j]}*A[j,i] = 0,                   %
%      delta=+-1, l=0,...,(k-3)/2                                         %
%      n is the number of observations, n_var is the dimension of x_phi   %
%      phi_i(x') is a column vector with the same dimension of x'         %
%      phi_i(x') = sum_{j=1}^{k}K(x',x_j)*A[j,i]                          %
%                   = [K(x',x_1),...,K(x',x_k)]*A[1:k,i]                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%initialization                                                           %
n = length(x);%number of observations                                     %
                                                                          %
row_A = zeros((n-k+1)*k,1);                                               %
col_A = zeros((n-k+1)*k,1);                                               %
values_A = zeros((n-k+1)*k,1);                                            %
                                                                          %
row_phi = zeros((n-k+1)*k,1);                                             %
col_phi = zeros((n-k+1)*k,1);                                             %
values_phi = zeros((n-k+1)*k,1);                                          %
                                                                          %
ii = 1;                                                                   % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%intermediate basis phi_[(k+1)/2 : (n-(k-1)/2)]                           %
for i = (k+1)/2 : (n-(k-1)/2)                                             %
    %%compute V_temp                                                      %
    V_temp = ...                                                          %
        [x(i-(k-1)/2 : i+(k-1)/2) - (x(i-(k-1)/2)+x(i-(k-1)/2))/2] ...    %
        .^([0:(k-3)/2]');                                                 %
                                                                          %
    %%compute A by solving linear equations                               %
    null_temp =  ...                                                      %
        [V_temp ...                                                       %
        .*exp(x(i-(k-1)/2 : i+(k-1)/2) - x(i-(k-1)/2)); ...               %
        V_temp ...                                                        %
        .*exp(-x(i-(k-1)/2 : i+(k-1)/2) + x(i-(k-1)/2)); ...              %
        [1,zeros(1,(k-1))]] ...                                           %
        \[zeros((k-1),1);1];                                              %
                                                                          %
    row_A(ii : ii+k-1) = i-(k-1)/2 : i+(k-1)/2;                           %
    col_A(ii : ii+k-1) = ones(k,1)*i;                                     %
    values_A(ii : ii+k-1) = null_temp;                                    %
                                                                          %
    %%compute phi_fun                                                     %
    x_mid = x(i-(k-1)/2 : i+(k-1)/2);                                     %
    phi_temp = matern_halfint(x_mid', x_mid', (k-2)/2, 1, rho) ...        %
        *null_temp;                                                       %
                                                                          %
    row_phi(ii : ii+k-1) = i-(k-1)/2 : i+(k-1)/2;                         %
    col_phi(ii : ii+k-1) = ones(k,1)*i;                                   %
    values_phi(ii : ii+k-1) = phi_temp;                                   %
                                                                          %   
    ii = ii + k;                                                          % 
end%end for loop over intermediate basis                                  %
                                                                          %
A = sparse(row_A, col_A, values_A, n, n);                                 %
phi = sparse(row_phi, col_phi, values_phi, n, n);                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%boundary basis phi_[1 : (k-1)/2] & phi_[n-(k-3)/2 : n]                   %
for i = 1 : (k-1)/2                                                       %
    %%compute V_temp_left & right                                         %
    V_temp_left = [x(1 : i+(k-1)/2) - (x(1)+x(i+(k-1)/2))/2] ...          %
        .^([0:(k-3)/2]');                                                 %
    V_temp_right = [x(n-k+1+i : n) - (x(n-k+1+i)+x(n))/2] ...             %
        .^([0:(k-3)/2]');                                                 %
                                                                          %
    %%compute A                                                           %
    if i == 1                                                             %
        %left-most A[:,1]                                                 %
        A(1 : i+(k-1)/2,i) = [V_temp_left ...                             %
            .*exp(x(1 : i+(k-1)/2) - x(1)); ...                           %
            [1,zeros(1,i+(k-3)/2)]] ...                                   %
            \[zeros(i+(k-3)/2,1);1];                                      %
        %right boundary A[:,n-(k-1)/2+1]                                  %
        A(n-k+1+i : n, n-(k-1)/2+i) = ...                                 % 
            [V_temp_right(1 : (k-1)/2-i, :) ...                           %
            .*exp(x(n-k+1+i : n) - x(n-k+1+i)); ...                       %
            V_temp_right ...                                              %
            .*exp(-x(n-k+1+i : n) + x(n-k+1+i)); ...                      %
            [1,zeros(1,k-1-i)]] ...                                       %
            \[zeros(k-1-i,1);1];                                          %
    elseif i == (k-1)/2                                                   %
        %right-most basis A[:,n]                                          %
        A(n-k+1+i : n, n-(k-1)/2+i) = [V_temp_right ...                   %
            .*exp(-x(n-k+1+i : n) + x(n-k+1+i)); ...                      %
            [1,zeros(1,k-1-i)]] ...                                       %
            \[zeros(k-1-i,1);1];                                          %
        %left boundary A[:,(k-1)/2]                                       %
        A(1 : i+(k-1)/2,i) = [V_temp_left ...                             %
            .*exp(x(1 : i+(k-1)/2) - x(1)); ...                           %
            V_temp_left(1 : i-1, :) ...                                   %
            .*exp(-x(1 : i+(k-1)/2) + x(1)); ...                          %
            [1,zeros(1,i+(k-3)/2)]] ...                                   %
            \[zeros(i+(k-3)/2,1);1];                                      %
    else                                                                  %
        %left boundary A[:,i]                                             %
        A(1 : i+(k-1)/2,i) = [V_temp_left ...                             %
            .*exp(x(1 : i+(k-1)/2) - x(1)); ...                           %
            V_temp_left(1 : i-1, :) ...                                   %
            .*exp(-x(1 : i+(k-1)/2) + x(1)); ...                          %
            [1,zeros(1,i+(k-3)/2)]] ...                                   %
            \[zeros(i+(k-3)/2,1);1];                                      %
        %right boundary A[:,n-(k-1)/2+i]                                  %
        A(n-k+1+i : n, n-(k-1)/2+i) = ...                                 %
            [V_temp_right(1 : (k-1)/2-i, :) ...                           %
            .*exp(x(n-k+1+i : n) - x(n-k+1+i)); ...                       %
            V_temp_right ...                                              %
            .*exp(-x(n-k+1+i : n) + x(n-k+1+i)); ...                      %
            [1,zeros(1,k-1-i)]] ...                                       %
            \[zeros(k-1-i,1);1];                                          %
    end%end if statement                                                  %
                                                                          %
    %%compute phi_fun                                                     %
    %left boundary basis                                                  %
    %distance_left = abs(x(1 : i+(k-1)/2)' - x(1 : i+(k-1)/2));            %
    x_index_left = 1 : i+(k-1)/2;                                       %
    x_left = x(x_index_left);
    phi(x_index_left, i) = ...                                            %
        matern_halfint(x_left', x_left', (k-2)/2, 1, rho) ...                %
        *A(1 : i+(k-1)/2, i);                                             %
                                                                          %
    %right boundary basis                                                 %
    %distance_right = abs(x(n-k+1+i : n)' - x(n-k+1+i : n));               %
    x_index_right = n-k+1+i : n;                                           %
    x_right = x(x_index_right);
    phi(x_index_right, n-(k-1)/2+i) = ...                                 %
        matern_halfint(x_right', x_right', (k-2)/2, 1, rho) ...               %
        *A(n-k+1+i : n, n-(k-1)/2+i);                                     %
end%end for loop over boundary basis                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end%end function
