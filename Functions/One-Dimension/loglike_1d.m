%% Compute log-likelihood of Parameter 'rho' for One Dimension
function [varargout]= loglike_1d(rho, k, X, y_fun, f)

%==========================================================================
%==========================================================================
%input: rho is lengthscale parameter in Matern kernel
%       k is the smoothness of Matern kernel, denotes(k-2)/2 Matern
%       X = [x1;x2;...;xN] is the vector of input, a N*1 matrix 
%       y_fun is true function of observations, y_fun(X)= [y1;y2;...;yN]
%
%optional input: f = [f1(.), f2(.), ..., fp(.)]
%   where mean = beta1*f1(.) + ... + betap*fp(.)
%
%output: loglike = log-likelihood(beta_hat,sigma2_hat,rho)
%
%optional output: 
%        beta_hat is updated parameter beta in MLE
%        sigma2_hat is updated parameter sigma2 in MLE
%==========================================================================
%==========================================================================

%check if input include f
if nargin < 5%~exist('f')
   f = @(x)ones(size(x,1),1);
    %By default we consider p = 1, F =[f(x1);f(x2);...;f(xN)] = [1;1;...;1]
    %mu_fun = beta*F = [beta;beta;...;beta], which is a N*1 matrix 
end
if nargin < 4
    error('Not enough input arguments.')
end


%update beta
[A, Phi] = compute_basis_solver(X', k, rho);
F = f(X);
w_beta = Phi'\(A'*F);
y = y_fun(X);
beta_hat = (w_beta'*F)\(w_beta')*y;

%update sigma2
w_sigma2 = Phi'\(A'*(y-F*beta_hat));
N = size(X,1);
sigma2_hat = w_sigma2'*(y-F*beta_hat)./N;

%compute log|K|
logdet_A = bdlogdet(A, (k-1)/2);
logdet_Phi = bdlogdet(Phi, (k-1)/2);
logdet_K = logdet_Phi - logdet_A;

%compute log-likelihood L(beta_hat,sigma2_hat,rho)
loglike = -N/2*log(sigma2_hat) - logdet_K/2 - 1/2*N;%(y-F*beta_hat)'*w_sigma2/sigma2_hat;
loglike = real(loglike);

varargout{1} = loglike;

if nargout >= 2
    varargout{2} = beta_hat;
end
if nargout >= 3
    varargout{3} = sigma2_hat;
end
if nargout > 3
    error('Incorrect output arguments.');
end

end%end loglike_1d function