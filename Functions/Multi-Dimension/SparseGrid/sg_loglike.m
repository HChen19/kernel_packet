%% Compute log-likelihood of Parameter 'rho' for High Dimensional Sparse Grid Design
function [varargout]= sg_loglike(rho, k, sg, sg_logdet, y_fun, f)

%==========================================================================
%==========================================================================
%input: rho is lengthscale parameter in Matern kernel
%       sg is the sparse grid design structure starting from max(d, eta-d+1)
%       sg_logdet is the sparse grid design structure starting from d
%       y_fun is true function of observations, y_fun(X)= [y1;y2;...;yN]
%
%optional input: 
%       f = [f1(.), f2(.), ..., fp(.)]
%   where mean = beta1*f1(.) + ... + betap*fp(.)
%
%output: loglike = log-likelihood L(beta_hat,sigma2_hat,rho)
%
%optional output: 
%        beta_hat is updated parameter beta in MLE
%        sigma2_hat is updated parameter sigma2 in MLE
%==========================================================================
%==========================================================================

%check if input include x_left and x_right
if nargin < 6%~exist('f')
   f = @(x)ones(size(x,1),1);
    %Here we consider p = 1, F =[f(x1);f(x2);...;f(xN)] = [1;1;...;1]
    %mu_fun = beta*F = [beta;beta;...;beta], which is a N*1 matrix 
end
if nargin < 5
    error('Not enough input arguments.')
end

X_set = sg.X_set;

%update beta
multiplier_beta = f;
[w_beta] = sg_w(sg, k, multiplier_beta, rho);
F = f(X_set);
y = y_fun(X_set);
beta_hat = (w_beta'*F)\(w_beta')*y;

%update sigma2
multiplier_sigma2 = @(x)y_fun(x)-f(x)*beta_hat;
[w_sigma2] = sg_w(sg, k, multiplier_sigma2, rho);
N = size(X_set,1);
sigma2_hat = w_sigma2'*(y-F*beta_hat)./N;
%sigma2_hat = 1;

%compute log|K|
log_K = splogdet(k, sg_logdet, rho);

%compute log-likelihood L(beta_hat,sigma2_hat,rho)
loglike = -N/2*log(sigma2_hat) - log_K/2 - (1/2)*N;%(y-F*beta_hat)'*w_sigma2/sigma2_hat;
loglike = real(loglike);

varargout{1} = loglike;

if nargout >= 2
    varargout{2} = beta_hat;
end
if nargout >= 3
    varargout{3} = sigma2_hat;
end
if nargout >= 4
    varargout{4} = log_K;
end
if nargout >= 5
    error('Incorrect output arguments.');
end

end%end  sg_loglike function
