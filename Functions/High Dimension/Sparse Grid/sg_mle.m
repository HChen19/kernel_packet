%% Compute updated parameter theta_hat, log-like L_hat and L_init in high dimensional MLE under sparse grid design
function [theta_hat,L_hat,L_init] = sg_mle(rho_init, k, sg, sg_logdet, y_fun, f)

%==========================================================================
%==========================================================================
%input: rho_init is initial lengthscale parameter in Matern kernel
%       k is the smoothness of Matern kernel, denotes(k-2)/2 Matern
%       sg is the sparse grid design structure starting from max(d, eta-d+1)
%       sg_logdet is the sparse grid design structure starting from d
%       y_fun is true function of observations, y_fun(X)= [y1;y2;...;yN]
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

d = sg.d;
%check if input include f
if nargin < 6%~exist('f')
   f = @(x)ones(size(x,1),1);
    %Here we consider p = 1, F =[f(x1);f(x2);...;f(xN)] = [1;1;...;1]
    %mu_fun = beta*F = [beta;beta;...;beta], which is a N*1 matrix 
end
if nargin < 5
    error('Not enough input arguments.')
end

%[L_init, beta_hat, sigma2_hat] = sg_loglike(rho_init, k, sg, sg_logdet, y_fun, f);
loglike_rho = @(rho) -sg_loglike(rho, k, sg, sg_logdet, y_fun, f);
L_init = - loglike_rho(rho_init);

%update rho
options = optimoptions('fmincon','Algorithm','sqp','Display','off');
rho_hat = fmincon(loglike_rho, rho_init,[],[],[],[],5e-2,5,[],options);
% options = optimoptions('fminunc', 'Display', 'off');
% rho_hat = fminunc(loglike_rho, rho_init,options);
[L_hat, beta_hat, sigma2_hat] = sg_loglike(rho_hat, k, sg, sg_logdet, y_fun, f);
%L_hat = - loglike_rho(rho_hat);

theta_hat.beta = beta_hat;
theta_hat.sigma2 = sigma2_hat;
theta_hat.rho = rho_hat;

end%end sg_mle function