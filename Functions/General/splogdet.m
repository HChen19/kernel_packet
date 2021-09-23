%% Compute log-determinant of covariance matrix under sparse grid design in high dimension
function [splogdet_M] = splogdet(k, sg_logdet, rho)

%==========================================================================
%==========================================================================
%input: sg is the sparse grid construction
%       k is the smoothness of Matern kernel, denotes(k-2)/2 Matern
%       sg_logdet is the sparse grid design structure starting from d
%
%optional input: default rho = sqrt(k-2)
%       rho is lengthscale parameter in Matern kernel
%
%output: splogdet_M = log|M|, M is Matern kernel matrix over sparse grid
%   design 
%==========================================================================
%==========================================================================

sg = sg_logdet;

d = sg.d;
eta = sg.eta;
design_fun = sg.design_fun;

X_tot = sg.X_tot;
ind_X_grid = sg.ind_X_grid;
dims = sg.dims;

%check if input include rho
if nargin < 3%~exist('rho')
    rho = sqrt(k-2)*ones(1,d);
end
if nargin < 2
    error('Not enough input arguments.')
end

%check if rho is d-dimensional
if length(rho) ~= d
    rho = rho(1)*ones(1,d);
end
%check if d <= eta
if d > eta
    error('Input argument d should be less than or equal to eta')
end

splogdet_M = 0;

for j_sum = d:eta%loop over j_sum
    [n_partition, j_arrow] = nsumk(d,j_sum);
    for partition_i = 1:n_partition %loop over partitions of j_sum
        
        logdet_K = zeros(1,d);
        logdet_K_sub = zeros(1,d);
        
        ind = ind_X_grid{j_sum-d+1, partition_i};
        X = X_tot(ind,:);
        dim_kron = dims{j_sum-d+1, partition_i};
        
        for dim_i = 1:d %loop over dimension
            
            X_ij = unique(X(:,dim_i))';
            X_sub_ij = design_fun([j_arrow(partition_i,dim_i)-1, sg.interval(1),sg.interval(2)]);
            dim_kron_sub(dim_i) = length(X_sub_ij);
            
            if dim_kron(dim_i) >= k
                [A, Phi] = compute_basis_solver(X_ij, k, rho(dim_i));
                logdet_A = bdlogdet(A, (k-1)/2);
                logdet_Phi = bdlogdet(Phi, (k-1)/2);
                logdet_K(dim_i) = logdet_Phi - logdet_A;
            else
                K = matern_halfint(X_ij', X_ij',(k-2)/2, 1, rho(dim_i));
                logdet_K(dim_i) = log(det(K));
            end
            
            if dim_kron_sub(dim_i) >= k
                [A_sub, Phi_sub] = compute_basis_solver(X_sub_ij, k, rho(dim_i));
                logdet_A_sub = bdlogdet(A_sub, (k-1)/2);
                logdet_Phi_sub = bdlogdet(Phi_sub, (k-1)/2);
                logdet_K_sub(dim_i) = logdet_Phi_sub - logdet_A_sub;
            elseif dim_kron_sub(dim_i) == 0
                logdet_K_sub(dim_i) = 0;
            else
                K_sub = matern_halfint(X_sub_ij', X_sub_ij', (k-2)/2, 1, rho(dim_i));
                logdet_K_sub(dim_i) = log(det(K_sub));
            end
           
        end%loop for loop over dimension
        
        dims_diff = dim_kron - dim_kron_sub;
        dims_prod = prod(dims_diff);
        dims_mvprod = dims_prod./dims_diff;
        splogdet_M = splogdet_M +  ...
            sum((logdet_K - logdet_K_sub).*dims_mvprod);

    end%end for loop over partitions of j_sum
end%end for loop over j_sum

end%end splogdet function