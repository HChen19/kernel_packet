function [C] = matern_halfint(X1, X2, nu, sigma, rho)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%input: X1: n1*d matrix of locations, d is the dimension, n1 is the # of X1
%       X2: n2*d matrix of locations, d is the dimension, n2 is the # of X2
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
%
%output: C: n1*n2 covariance matrix 
%
%see the wikipedia page of Matern covariance function:
% https://en.wikipedia.org/wiki/Mat%C3%A9rn_covariance_function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n1 = size(X1,1);
n2 = size(X2,1);
d = size(X1,2);

if nargin < 5
    rho = sqrt(2*nu)*ones(1,d);
end
if nargin < 4
    sigma = 1;
end
if nargin < 3
    error('Not enough input arguments, requires at least 2 inputs.');
end

if length(rho) ~= d
    rho = rho(1)*ones(1,d);
end

rho = rho';
X1_scaled = X1./repmat(rho,1,n1)';
X2_scaled = X2./repmat(rho,1,n2)';
F(1,:,:) = X2_scaled';
diff_val = abs(repmat(X1_scaled,[1,1,n2])-repmat(F,[n1,1,1]));
power = squeeze(sum(diff_val,2));

p = nu -1/2;
sum_i = zeros(size(d));
if p == 0%1/2 Matern kernel
    C = sigma^2.*exp(-power);
elseif p == 1%3/2 Matern kernel
    C = sigma^2.*squeeze(prod(1+sqrt(3)*diff_val,2)) .* exp(-sqrt(3).*power);
elseif p == 2%5/2 Matern kernel
    C = sigma^2.*squeeze(prod(1+sqrt(5)*diff_val+5*diff_val.^2./3,2)) ...
        .*exp(-sqrt(5).*power);
else
    for i = 0:p
        sum_i = sum_i + factorial(p+i)./(factorial(i).*factorial(p-i)) ...
            .*(2.*sqrt(2*p+1).*diff_val).^(p-i);
    end
    %term = factorial(p:2*p)./(factorial(0:p).*factorial(p:-1:0)) ...
    %   .*(2*sqrt(2*p+1).*d./rho).^(p:-1:0);
    C = sigma^2 .* exp(-sqrt(2*p+1).*power) .* factorial(p)./factorial(2*p) ...
    .*squeeze(prod(sum_i,2)); 
end

if size(C,1) ~= n1
    C = C';
end

end% end matern_halfint function