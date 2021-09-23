%% Compute log-determinant of band matrix
function [logdet_M] = bdlogdet(M, bandwidth)

%==========================================================================
%==========================================================================
%input: M is a band matrix
%       bandwidth is bandwidth of band matrix M
%
%output: logdet_M = log|M|
%==========================================================================
%==========================================================================
logdet_M = 0;
width = bandwidth+1;
q = floor(length(M)/width);
q_mod = length(M) - width*q;
index_A = 1 : width;
A_hat = M(index_A,index_A);

for i = 1:q
    [L,U,P] = lu(A_hat);
    logdet_M = logdet_M + log(det(U)) + log(det(P'));
    inv_U = inv(U);
    
    if i < q
        index_Ahat = 1+width*i : width*(i+1);
        index_B_row = 2+width*(i-1) : width*i;
        index_B_col = 1+width*i : width*(i+1)-1;
        index_C_row = index_B_col;
        index_C_col = index_B_row;
    else
        index_Ahat = length(M)-q_mod+1 : length(M);
        index_B_row = 2+width*(i-1) : width*i;
        index_B_col = length(M)-q_mod+1 : length(M);
        index_C_row = index_B_col;
        index_C_col = index_B_row;
    end
    
    B = M(index_B_row,index_B_col);
    C = M(index_C_row,index_C_col);
    X = C*inv_U(end-size(C,2)+1:end,end-size(C,2)+1:end);
    G = (P'*L)\[zeros(1,size(B,2));B];
    G2 = G(2:end,:);
    if i < q
        A_hat = M(index_Ahat,index_Ahat) - [X*G2,zeros(width-1,1);zeros(1,width)];
    else
        A_hat = M(index_Ahat,index_Ahat) - [X*G2];
    end
end

[L,U,P] = lu(A_hat);
%det_M = det_M* det(U)*det(P');
logdet_M = logdet_M + log(det(U)) + log(det(P'));

end%end bdlogdet function