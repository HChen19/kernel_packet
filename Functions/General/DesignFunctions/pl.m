%% Return a vector with Plumlee's design points
function [pl_j] = pl(j,x_1,x_n)

% pl.vals = the point coordinates of the set pl_j
% pl.subs = the subscripts of vector pl_j
%we suppose pl_{i,j} is the same for different dimensin i(1<=i<=d, d is the dimension)
% j = the number of points in set pl_{i,j}
% [x_1,x_n] is the interval of the dimesion i
%
% Here we use Plumplee's sparse grid design
X = zeros(16,2);
X(1,:) = [0.5, 0.5];
X(2,:) = [.125,.875];
X(3,:) = [.25,.75];
X(4,:) = [0,1];
X(5,:) = [.375,.625];
X(6,:) = [.25-1/16,.75+1/16];
X(7,:) = [1-1/16,0+1/16];
X(8,:) = [.5-1/16,.5+1/16];
X(9,:) = [0.3125,   0.6875];
X(10,:)= [.25+1/32,.75-1/32];
X(11,:) = [0+1/32,1-1/32];
X(12,:) = [.875+1/32,.125-1/32];
X(13,:) = [.875-1/32,.125+1/32];
X(14,:) = [.375+1/32,.625-1/32];
X(15,:) = [.5-1/32,.5+1/32];
X(16,:) = [.25-1/32,.75+1/32];

if j == 0
    pl_j = [];
elseif j <= 16
    pl_j = x_1 + unique(X(1:j,:))' .* (x_n-x_1);
else
    error('Out of range subscript')
end

end%end Pluemplee's (pl) function