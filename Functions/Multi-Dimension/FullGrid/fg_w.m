%% Compute w for Lattice Design
function [w_value]= fg_w(fg, k, multiplier, rho)

%==========================================================================
%==========================================================================
%input: fg is the full grid design structure
%       k is the smoothness of Matern kernel, denotes(k-2)/2 Matern
%       multiplier is a function whose value matrix multiplies kron product
%
%optional input: default rho = sqrt(k-2),default x_left=0,default x_right=1
%       rho is lengthscale parameter in Matern kernel
%       [x_left, x_right]^d is the interval our input sampled in
%
%output: w_value = kron(K_inv{1}, K_inv{2}, ..., K_inv{d})*m
%        w_value is a N*p matrix of value of w, where N is the number of 
%   full grid points,N=(2^(eta-1)-1)^d, p=size(m, 2),
%   where K_inv{i} is inverse of covariance matrix in dimension i
%   m = multiplier(fg.X_set)      
%
%==========================================================================
%==========================================================================

X_1d = fg.X_grid;
m = multiplier(fg.X_set);
dims = [size(m,2), fg.dims];
dims_flip = fliplr(dims);
d = fg.d;
eta = fg.eta;

%check if input include rho
if nargin < 4%~exist('rho')
    rho = sqrt(k-2);
end
if nargin < 3
    error('Not enough input arguments.')
end

%check if d <= eta
if d > eta
    error('Input argument d should be less than or equal to eta')
end

vt = reshape(m, dims_flip);
for i = 1:d%loop over dimension
    if length(X_1d) >= k
        [A, Phi] = compute_basis(X_1d, k, 'solver', rho);
        vt = refold(A* (Phi\unfold(vt, i+1, dims_flip)), ...
            i+1, dims);
    else
        K_inv = inv(matern_halfint(X_1d', X_1d', (k-2)/2, 1, rho));
        vt = refold(K_inv*unfold(vt, i+1, dims_flip), i+1, dims);
    end%end if statement
end%end for loop over dimension
w_value = reshape(vt,[],dims(1));

end%end fg_w function