function [varargout]= sg_w(sg, k, multiplier, rho)

%==========================================================================
%==========================================================================
%input: sg is the sparse grid design structure
%       k is the smoothness of Matern kernel, denotes(k-2)/2 Matern
%       multiplier is a function whose value matrix multiplies kron product
%
%optional input: default rho = sqrt(k-2),default x_left = 0,default x_right = 1
%       rho is lengthscale parameter in Matern kernel
%
%output: w_value = kron(K_inv{1}, K_inv{2}, ..., K_inv{d})*m, 
%        w_value is a N*p matrix of nonzero value of w, where N is the 
%   number of sparse grid points, p = size(m, 2),     
%   where K_inv{i} is inverse of covariance matrix in dimension i,
%   m = multiplier(sg.X_set)
%
%optional output:       
%        w_sptensor is a 1*p cell of sparse tensor, its nonzero index is 
%   find(w_sptensor{i}), nonzero value is w_sptensor{i}(find(w_sptensor{i})) 
%==========================================================================
%==========================================================================

X_tot = sg.X_tot;
subs_X_tot = sg.subs_X_tot;
ind_X_grid = sg.ind_X_grid;
dims = sg.dims;
d = sg.d;
eta = sg.eta;

%check if input include rho
if nargin < 4%~exist('rho')
    rho = sqrt(k-2)*ones(1,d);
end
if nargin < 3
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

%Initialization
subs_w = [];%subscripts of the nonzeros in sparse tensor w
vals_w = [];%values of the nonzeros in sparse tensor w

for j_sum = max(d,eta-d+1):eta
    [n_partition, j_arrow] = nsumk(d,j_sum);
    a = (-1)^(eta-j_sum) * nchoosek(d-1, eta-j_sum);
    for partition_i = 1:n_partition %loop over partitions of eta(differnet j_arrow for the same eta)
        ind = ind_X_grid{j_sum-max(d,eta-d+1)+1, partition_i};
        X = X_tot(ind,:);
        m = multiplier(X);
        dim_kron = dims{j_sum-max(d,eta-d+1)+1, partition_i};
        dim_total = [size(m,2),dim_kron];
        dims_flip = fliplr(dim_total);
        
        m_ts = reshape(m, dims_flip);
        for dim_i = 1:d%loop over dimension
            X_ij = unique(X(:,dim_i))'; 
            if dim_kron(dim_i) >= k
                [A, Phi] = compute_basis_solver(X_ij, k, rho(dim_i));
                m_ts = refold(A* (Phi\unfold(m_ts, dim_i+1, dims_flip)), ...
                    dim_i+1, dim_total);
            elseif length(m_ts) > 1
                [Q,~] = chol(matern_halfint_1d(X_ij', X_ij', ...
                    (k-2)/2, 1, rho(dim_i)));
                K_inv =Q\(Q')^(-1);
%                 K_inv = inv(matern_halfint(X_ij', X_ij', ...
%                     (k-2)/2, 1, rho(dim_i)));
                m_ts = refold(K_inv*unfold(m_ts, dim_i+1, dims_flip), dim_i+1, dim_total);
            end%end if statement
        end%end for loop over dimension
       kron_rs = reshape(m_ts,[],size(m,2));
        
%         m_ts = tensor(reshape(m, dims_flip));
%         for dim_i = 1:d %loop over dimension
%             X_ij = unique(X(:,dim_i))';       
%             if dim_kron(dim_i) >= k
%                 [A, Phi] = compute_basis(X_ij, k, 'solver', rho(dim_i));
%                 m_mat = tenmat(m_ts,d-dim_i+1);
%                 m_mat = A*(Phi\double(m_mat));
%                 m_ts = reshape(tensor(m_mat),[dims_flip(d-dim_i+1),dims_flip(setdiff(1:d+1,d-dim_i+1))]);
%                 if dim_i < d
%                     m_ts = permute(m_ts,[2:d-dim_i+1,1,d-dim_i+2:d+1]);
%                 end
%                 
%                 %K_inv = ( (Phi')\(A') )';
%                 %m_ts = ttm(m_ts, inv(Phi), d-dim_i+1);
%                 %m_ts = ttm(m_ts, A, d-dim_i+1);
%                 %K_inv{dim_i} = A{dim_i}*inv(Phi{dim_i});
%             else
%                 K_inv = inv(matern_halfint(X_ij', X_ij', ...
%                     (k-2)/2, 1, rho(dim_i)));
%                 m_ts = ttm(m_ts,K_inv,d-dim_i+1);
%             end
%         end%loop for loop over dimension
%         kron_rs = double(tenmat(m_ts,d+1)');
        
        subs_w_temp = subs_X_tot(ind,:);
        w_temp = a*kron_rs;
        
        subs_w = [subs_w;subs_w_temp];
        vals_w = [vals_w;w_temp];

    end%end loop over partitions of eta
end%end loop over jsum

w_sptensor = cell(1, size(m,2));
%w_value = [];
for i = 1: size(m,2)
    w_sptensor{i} = sptensor(subs_w,vals_w(:,i));
    num_nnz = nnz(w_sptensor{i});
    w_value(1:num_nnz,i) = w_sptensor{i}.vals;
end

varargout{1} = w_value;

if nargout >= 2
    varargout{2} = w_sptensor;
end
if nargout > 2
    error('Incorrect output arguments.');
end

end%end sg_w function