%% Compute y_pred in High Dimensional GPR, y_pred = posterior mean
function [y_pred]= compute_pred(xnew, w_value, X_set, k, rho, mu_fun)

%==========================================================================
%==========================================================================
%input: xnew is input we want to predict, m*d matrix
%       w_value = inv(K)*(y-mu), where K is covariance Matern matrix
%       X_set is a N*d matrix for input coordinates of sparse grid 
%   points
%       k is the smoothness of Matern kernel, denotes(k-2)/2 Matern
%
%optional input: default mu_fun = @(x) zeros(size(x,1),1), is zero function
%                default rho = sqrt(k-2), lengthscale of Matern kernel
%   
%output: y_pred is prediction vlaue on xnew, m*1 matrix
%==========================================================================
%==========================================================================

%check if input include mu_fun
if nargin < 6%~exist('mu_fun')
    mu_fun = @(x) zeros(size(x,1),size(w_value,2));
end
if nargin < 5%~exist('rho')
    rho = sqrt(k-2)*ones(1,size(X_set,2));
end
if nargin < 4
    error('Not enough input arguments.')
end
%check if rho is d-dimensional
if length(rho) ~= size(X_set,2)
    rho = rho(1)*ones(1,size(X_set,2));
end

y_pred = mu_fun(xnew) + matern_halfint(xnew, X_set, (k-2)/2, 1, rho)*w_value;
%find(w) extract nonzeros of w by the ascending order of dimension

end%end compute_pred function