function [C] = matern_halfint_1d(X1, X2, nu, sigma, rho)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%input: X1: n1*1 matrix of locations, n1 is the # of X1
%       X2: n2*1 matrix of locations, n2 is the # of X2
%       nu: roughness, controls the sharpness of ridges in the covariance 
%           function, which ultimately affects the roughness (smoothness) 
%           of realizations
%       sigma: amplitude, controls the scaling of the output along the 
%              y-axis. This parameter is just a scalar multiplier and is 
%              therefore usually left out of implementations of the Matèrn 
%              function
%       rho: lengthscale, complements the amplitude by scaling realizations 
%            on the x-axis. Larger values push points closer together along 
%            this axis
%see the wikipedia page of Matern covariance function:
% https://en.wikipedia.org/wiki/Mat%C3%A9rn_covariance_function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 5
    rho = sqrt(2*nu);
end
if nargin < 4
    sigma = 1;
end
if nargin < 2
    error('Not enough input arguments, requires at least 2 inputs.');
end

d = pdist2(X1,X2);%d = Euclidean distance
p = nu -1/2;
sum = zeros(size(d));
if p == 0%1/2 Matern kernel
    C = sigma^2.*exp(-d./rho);
elseif p == 1%3/2 Matern kernel
    C = sigma^2.*(1+sqrt(3).*d./rho).*exp(-sqrt(3).*d./rho);
elseif p == 2%5/2 Matern kernel
    C = sigma^2.*(1+sqrt(5).*d./rho+5.*d.^2./(3*rho^2)) ...
        .*exp(-sqrt(5).*d./rho);
else
    for i = 0:p
        sum = sum + factorial(p+i)./(factorial(i).*factorial(p-i)) ...
            .*(2.*sqrt(2*p+1).*d./rho).^(p-i);
    end
    %term = factorial(p:2*p)./(factorial(0:p).*factorial(p:-1:0)) ...
    %   .*(2*sqrt(2*p+1).*d./rho).^(p:-1:0);
    C = sigma^2.*exp(-sqrt(2*p+1).*d./rho).*factorial(p)./factorial(2*p) ...
    .*sum; 
end