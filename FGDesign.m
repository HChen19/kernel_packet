%% Example with Full Grid Design Functions lt_w.m and compute_pred.m
d = 2; %dimension
eta = 4; %level of sparse grid construction, which eta>= d
k = 5; %smootheness of Matern kernel
rho = sqrt(k-2);%lengthscale of Matern kernel
x_left = 0;
x_right = 1;

y_fun = @(x) sum( (x) .^2, 2)/4000 - prod(cos( (x) ./sqrt(1:d)), 2) + 1; %Griewank function, mean = 0.5, interval = []
mu_fun = @(x) zeros(size(x,1),1);
multiplier = @(x)y_fun(x)-mu_fun(x);

X_1d = x_left + (x_right-x_left).*[1:2^(eta-1)-1]./2^(eta-1);
%X_1d is coordinate of points in one dimension, lentgh(X_1d)=(2^(eta-1)-1)^d
X_grid = cell(1,d);
[X_grid{1:d}] = ndgrid(X_1d);
X_grid = reshape(cat(d,X_grid{:}),[],d);
fg.X_set = sortrows(X_grid);%(2^(eta-1)-1)^d * d matrix, points of full grid
fg.X_grid = X_1d;
fg.dims = ones(1,d)*(2^(eta-1)-1);
fg.d = d;
fg.eta = eta;

% compute w
tic
[w_fg] = fg_w(fg, k, multiplier, rho);
toc1 = toc;
N_FG = (2^(eta-1)-1)^d;%N is the number of full grid input obervations
fprintf('elapsed time of computing w_fg for d=%1.0f, eta=%1.0f and N=%1.0f is: %.8f seconds. \n',d,eta,N_FG,toc1')
% uncomment the following line to get the plot of X_set_unique
%scatter(X_lt(:,1),X_lt(:,2),'o','MarkerFaceColor', 'b');

% compute y_pred
m = 1e3;
xnew = x_left + rand(m,d).*(x_right-x_left);%m inputs in d dimension, m*d matrix
y_true = y_fun(xnew);
tic;
y_pred_fg = compute_pred(xnew, w_fg, fg.X_set, k, rho);
toc2 = toc;
mse_fg = mean((y_true-y_pred_fg).^2);
fprintf('elapsed time of predicting y_pred_sg for d=%1.0f, eta=%1.0f and N=%1.0f, mse=%0.8f is: %.8f seconds. \n',d,eta,N_FG, mse_fg, toc1+toc2')

figure;
scatter3(xnew(:,1),xnew(:,2),y_true(:,1))
hold on;
scatter3(xnew(:,1),xnew(:,2),y_pred_fg(:,1))
xlabel('X1')
ylabel('X2')
zlabel('mean(x,2).*sin(mean(x,2))')