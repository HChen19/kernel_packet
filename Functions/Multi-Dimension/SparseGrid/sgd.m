function [sg] = sgd(d, eta, design_fun, x_left, x_right, j_sum_start)

%==========================================================================
%==========================================================================
%input: d is the dimension of input, d >= 2
%       eta is level of sparse grid construction
%       design_fun is the design that used in one dimension (Here we assume
%   the sparse grid design is same in each dimension)
%
%optional input: default x_left = 0, default x_right = 1
%                default j_sum_start = 'start=max'        
%       [x_left, x_right]^d is the interval our input sampled in
%
%output: sg is a structure of sparse grid design
%        sg.X_tot: a M*d matrix, the total points of sparse grid design 
%   with repetition in this smolyak iteration, where M is the number of 
%   points without using unique
%        sg.X_set: a N*d matrix, the set of points of sparse grid desgin 
%   without repetition, where N is the number of points in the sparse grid
%   construction
%        sg.subs_X_tot: a M*d matrix, subscriptions of sg.X_tot
%        sg.ind_X_grid: a (eta-max(d,eta-d+1)+1)*(n_partition) cell array,use  
%   sg.subs_X_tot(sg. ind_X_grid{j_sum-j_sum_start+1, partition_i},1:d)
%   to extract subscriptions of X_grid in each smolyak iteration
%        sg.dims: a (eta-j_sum_start+1)*(n_partition) cell array, use
%   sg.dims{j_sum-j_sum_start+1, partition_i} to extract the number of
%   points projected in each dimension in each smolyak iteration
%==========================================================================
%==========================================================================

%check if input include 'j_sum_start' variable
if nargin < 6%~exist('j_sum_start','var')
    % sisth parameter does not exist, so default it to something
    j_sum_start = 'start=max';
end
%check if input include x_left and x_right
if nargin < 5%~exist('x_left','x_right')
    x_left = 0;
    x_right = 1;
end
if nargin < 3
    error('Not enough input arguments.')
end

%check if d <= eta
if d > eta
    error('Input argument d should be less than or equal to eta')
end

%Initialization
X_tot = [];%total points of sparse grid design with repetition
ind_X_grid = {};%indices of points in this smolyak iteration
dims = {};%dimensions in each iteration

ii = 1;

if strcmp(j_sum_start, 'start=max')
    j_sum_start = max(d,eta-d+1);
elseif strcmp(j_sum_start, 'start=d')
    j_sum_start = d;
end

for j_sum = j_sum_start:eta
    [n_partition, j_arrow] = nsumk(d,j_sum);
    for partition_i = 1:n_partition %loop over partitions of eta(differnet j_arrow for the same eta)
        X = cell(1,d);%sparse grid points
        dim_kron = zeros(1,d);%dim_kron(i) is the number of observation points in dimension i
        
        for dim_i = 1:d %loop over dimension
            [X{dim_i}] = design_fun([j_arrow(partition_i,dim_i), x_left,x_right]);
            dim_kron(dim_i) = length(X{dim_i});
        end%loop for loop over dimension
        [X{1:d}] = ndgrid(X{:});
        X = reshape(cat(d,X{:}),[],d);
        X = sortrows(X);
        %X is ()*d matirx, coordinates of points in j_sum for loop
        
        X_tot = [X_tot;X];%set of all the points including same points
        %X_grid{j_sum-max(d,eta-d+1)+1, partition_i} = X;
        ind_X_grid{j_sum-j_sum_start+1, partition_i} = ...
            ii : ii + size(X,1) - 1;
        ii = ii + size(X,1);
            %2^(eta-1)*(X-x_left)./(x_right-x_left);
        dims{j_sum-j_sum_start+1, partition_i} = dim_kron;
        
    end%end loop over partitions of eta
end%end loop over jsum

[X_set,~,~] = unique(X_tot,'rows','sorted');
[~,~,ic] = unique(X_tot);
subs_X_tot = reshape(ic,[],d);

sg.X_tot = X_tot;
sg.X_set = X_set;
sg.subs_X_tot = subs_X_tot;
sg.ind_X_grid = ind_X_grid;
sg.dims = dims;

sg.d = d;
sg.eta = eta;
sg.design_fun = design_fun;
sg.interval = [x_left, x_right];

end%end sgd function