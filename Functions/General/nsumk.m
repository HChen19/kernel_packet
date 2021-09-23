%% nsumk(n,k) return a partition matirx where n positive integers summing to k
function [m,x] = nsumk(n,k)

% NSUMK Number and listing of positive integer n-tuples summing to k
%    M = NSUMK(N,K) where N and K are positive integers returns M=nchoosek(K-1,N-1)
%    This is the number of ordered N-tuples of positive integers summing to K
%
%    [M,X] = NSUMK(N,K) produces a matrix X with
%    nchoosek(K-1,N-1) rows and n columns. Each row comprises
%    positive integers summing to k. The ordering of rows follows the
%    same convention as NCHOOSEK, which is undocumented but with some
%    reliability appears to be lexicographic. The reverse of this presumed ordering
%    is a natural way to list coefficients of polynomials in N variables of degree K.
%    As per nchoosek, this syntax is only practical for situations where N is
%    less than about 15.
% 
%  EXAMPLES:   m = nsumk(3,6), m=10 is the number of partition ways   
%              [~,x] = nsumk(3,6) returns a 10 x 3 matrix x in which rows
%              sum to 6
 
if n <= k && isscalar(n) && isscalar(k) && nargout<=1
    m = nchoosek(k-1,n-1);
elseif isscalar(n) && isscalar(k) && nargout==2
    m = nchoosek(k-1,n-1);
    dividers = [zeros(m,1),nchoosek((1:(k-1))',n-1),ones(m,1)*(k)];
    x = diff(dividers,1,2);
else
    error('nsumk anticipates scalar k and n or n is greater than k');
    
end
end