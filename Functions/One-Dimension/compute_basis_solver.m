%% compute basis function phi(x) and coefficient band matrix A
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
        [x(i-(k-1)/2 : i+(k-1)/2) - (x(i-(k-1)/2)+x(i-(k-1)/2))/2] ...    
        .^([0:(k-3)/2]');                                                 %
                                                                          %
    %%compute A by solving linear equations                               %
    null_temp =  ...                                                      %
        [V_temp ...                                                       %
        .*exp((x(i-(k-1)/2 : i+(k-1)/2) - x(i-(k-1)/2)) ...%
        .*sqrt(k-2)./rho); ...               %
        V_temp ...                                                        %
        .*exp((-x(i-(k-1)/2 : i+(k-1)/2) + x(i-(k-1)/2)) ...%
        .*sqrt(k-2)./rho); ...              %
        [1,zeros(1,(k-1))]] ...                                           %
        \[zeros((k-1),1);1];                                              %
                                                                          %
    row_A(ii : ii+k-1) = i-(k-1)/2 : i+(k-1)/2;                           %
    col_A(ii : ii+k-1) = ones(k,1)*i;                                     %
    values_A(ii : ii+k-1) = null_temp;                                    %
                                                                          %
    %%compute phi_fun                                                     %
    x_mid = x(i-(k-1)/2 : i+(k-1)/2);                                     %
    phi_temp = matern_halfint_1d(x_mid', x_mid', (k-2)/2, 1, rho) ...        %
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
            .*exp((x(1 : i+(k-1)/2) - x(1)).*sqrt(k-2)./rho); ...         %
            [1,zeros(1,i+(k-3)/2)]] ...                                   %
            \[zeros(i+(k-3)/2,1);1];                                      %
        %right boundary A[:,n-(k-1)/2+1]                                  %
        A(n-k+1+i : n, n-(k-1)/2+i) = ...                                 % 
            [V_temp_right(1 : (k-1)/2-i, :) ...                           %
            .*exp((x(n-k+1+i : n) - x(n-k+1+i)).*sqrt(k-2)./rho); ...     %
            V_temp_right ...                                              %
            .*exp((-x(n-k+1+i : n) + x(n-k+1+i)).*sqrt(k-2)./rho); ...    %
            [1,zeros(1,k-1-i)]] ...                                       %
            \[zeros(k-1-i,1);1];                                          %
    elseif i == (k-1)/2                                                   %
        %right-most basis A[:,n]                                          %
        A(n-k+1+i : n, n-(k-1)/2+i) = [V_temp_right ...                   %
            .*exp((-x(n-k+1+i : n) + x(n-k+1+i)).*sqrt(k-2)./rho); ...    %
            [1,zeros(1,k-1-i)]] ...                                       %
            \[zeros(k-1-i,1);1];                                          %
        %left boundary A[:,(k-1)/2]                                       %
        A(1 : i+(k-1)/2,i) = [V_temp_left ...                             %
            .*exp((x(1 : i+(k-1)/2) - x(1)).*sqrt(k-2)./rho); ...         %
            V_temp_left(1 : i-1, :) ...                                   %
            .*exp((-x(1 : i+(k-1)/2) + x(1)).*sqrt(k-2)./rho); ...        %
            [1,zeros(1,i+(k-3)/2)]] ...                                   %
            \[zeros(i+(k-3)/2,1);1];                                      %
    else                                                                  %
        %left boundary A[:,i]                                             %
        A(1 : i+(k-1)/2,i) = [V_temp_left ...                             %
            .*exp((x(1 : i+(k-1)/2) - x(1)).*sqrt(k-2)./rho); ...         %
            V_temp_left(1 : i-1, :) ...                                   %
            .*exp((-x(1 : i+(k-1)/2) + x(1)).*sqrt(k-2)./rho); ...        %
            [1,zeros(1,i+(k-3)/2)]] ...                                   %
            \[zeros(i+(k-3)/2,1);1];                                      %
        %right boundary A[:,n-(k-1)/2+i]                                  %
        A(n-k+1+i : n, n-(k-1)/2+i) = ...                                 %
            [V_temp_right(1 : (k-1)/2-i, :) ...                           %
            .*exp((x(n-k+1+i : n) - x(n-k+1+i)).*sqrt(k-2)./rho); ...                       %
            V_temp_right ...                                              %
            .*exp((-x(n-k+1+i : n) + x(n-k+1+i)).*sqrt(k-2)./rho); ...                      %
            [1,zeros(1,k-1-i)]] ...                                       %
            \[zeros(k-1-i,1);1];                                          %
    end%end if statement                                                  %
                                                                          %
    %%compute phi_fun                                                     %
    %left boundary basis                                                  %
    x_index_left = 1 : i+(k-1)/2;                                         %
    x_left = x(x_index_left);                                             %
    phi(x_index_left, i) = ...                                            %
        matern_halfint_1d(x_left', x_left', (k-2)/2, 1, rho) ...             %
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
