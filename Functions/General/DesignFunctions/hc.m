%% Return a vector with hyperbolic cross points(bisection)
function [hc_j] = hc(j,x_1,x_n)

% hc.vals = the point coordinates of the set hc_j
% hc.subs = the subscripts of vector hc_j
%we suppose hc_{i,j} is the same for different dimensin i(1<=i<=d, d is the dimension)
% j = the number of points in set hc_{i,j}
% [x_1,x_n] is the interval of the dimesion i
%
% Here we use grid design with hyperbolic cross points (bisection)
if j == 0
    hc_j = [];
else
    hc_j = ([(2^j-1):-1:1]*x_1 + [1:(2^j-1)]*x_n) * 0.5^j;
end

end%end hyperbolic cross (hc) function