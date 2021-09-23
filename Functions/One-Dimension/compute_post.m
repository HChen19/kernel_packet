%% Compute Posterior Mean and Covariance in One Dimensional GPR 
function [varargout] = compute_post(xnew, x, y, k, rho, mu_fun)

% GP Regression
%==========================================================================
% Assumptions:
% y_i = f(x_i) + epsilon_i, i=1,...,n. epsilon_i = 0 with noise off
% f ~ GP(mean, cov), mean and cov are functions
% epsilon_i ~ N (0, sigma^2), i.i.d.
%
% Input:
% x = [x1,x2,...,xn] belongs to R^n, n*1 matrix
% y = [y1,y2,...,yn]' belongs to R^n, n*1 matrix
% x_new = [x_new1,...,x_new_m] belongs to R^m
%
% Posterior:
% f|y ~ GP(mean_new, cov_new)
% mean_new(x_i) = mean(x_i) + cov_{x_i,x}* (cov_{x,x}+sigma^2*I_n)^(-1)* 
%                (y-mean(x))
% cov_new(x_i,x_j) = cov(x_i,x_j) - cov_{x_i,x}* (cov_{x,x}+sigma^2*I_n)^(-1) 
%                * cov_{x,x_j}'
%==========================================================================
if nargin < 6
    mu_fun = @(x) zeros(length(x),1);%mean function
end
if nargin < 5
    rho = sqrt(k-2);
end
if nargin < 4
    error('Not enough input arguments.');
end

xnew = sort(xnew);

mean_x = mu_fun(x);
mean_xnew = mu_fun(xnew);

[A, Phi_x, Phi_xnew] = compute_basis_solver_new(x, k, xnew, rho);

%compute mean_new
w = (Phi_x) \ (y - mean_x);
mean_new = mean_xnew + Phi_xnew*w;

varargout{1} = mean_new;

if nargout >= 2
    %compute cov_new
    cov_xnew = matern_halfint_1d(xnew, xnew, (k-2)/2, 1, rho);
    w_cov = (A'*Phi_x) \ (Phi_xnew');
    cov_new = cov_xnew - Phi_xnew*w_cov;
    varargout{2} = cov_new;
end

end%end compute_post function
