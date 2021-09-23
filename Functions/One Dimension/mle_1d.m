%% Compute updated parameter theta_hat, log-like L_hat and L_init in one dimensional MLE
function [theta_hat,L_hat,L_init] = mle_1d(rho_init, k, X, y_fun, f)

%==========================================================================
%==========================================================================
%input: rho_init is initial lengthscale parameter in Matern kernel
%       k is the smoothness of Matern kernel, denotes(k-2)/2 Matern
%       X = [x1;x2;...;xN] is the vector of input, a N*1 matrix 
%       y_fun is true function of observations, y_fun(X)= [y1;y2;...;yN]n
%       f = [f1(.), f2(.), ..., fp(.)]
%   where mean = beta1*f1(.) + ... + betap*fp(.)
%
%optional input: 
%       f = [f1(.), f2(.), ..., fp(.)]
%   where mean = beta1*f1(.) + ... + betap*fp(.)
%
%output: 
%       theta_hat = [beta_hat, sigma2_hat, rho_hat], use theta_hat.beta to
%   extract beta_hat
%       L_hat = log-likelihood L(beta_hat,sigma2_hat,rho_hat)
%       L_init = log-likelihood L(beta_hat,sigma2_hat,rho_init)
%
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

[L_init, beta_hat, sigma2_hat] = loglike_1d(rho_init, k, X, y_fun, f);
loglike_rho = @(rho) -loglike_1d(rho, k, X, y_fun, f);


%update rho
options = optimoptions('fmincon','Algorithm','sqp','Display','off');
rho_hat = fmincon(loglike_rho, rho_init,[],[],[],[],5e-2,50,[],options);

L_hat = - loglike_rho(rho_hat);

theta_hat.beta = beta_hat;
theta_hat.sigma2 = sigma2_hat;
theta_hat.rho = rho_hat;

end%end MLE_1d function