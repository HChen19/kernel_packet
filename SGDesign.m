%% Example with Sparse Grid Design Functions sg_w.m and compute_pred.m
d = 2; %dimension
eta = 10; %level of sparse grid construction, which satisfies eta>= d
k = 5; %smootheness of Matern kernel
rho = sqrt(k-2);%lengthscale of Matern kernel
x_left = 0;
x_right = eta - 1;

y_fun = @(x) sum( (x) .^2, 2)/4000 - prod(cos( (x) ./sqrt(1:d)), 2) + 1; %Griewank function, mean = 0.5, interval = []
mu_fun = @(x) zeros(size(x,1),1);
multiplier = @(x)y_fun(x)-mu_fun(x);

design_fun = @(x)hc(x(1),x(2),x(3));
%design_fun = @(x)pl(x(1),x(2),x(3));

[sg] = sgd(d, eta, design_fun, x_left, x_right);
X_tot = sg.X_tot;
X_set = sg.X_set;
subs_X_tot = sg.subs_X_tot;
ind_X_grid = sg.ind_X_grid;
dims = sg.dims;
d = sg.d;
eta = sg.eta;
% uncomment the following lines to get the plot of X_set
% figure;
% scatter(sg.X_set(:,1),sg.X_set(:,2),'o','MarkerFaceColor', 'b');

% compute w
tic;
[w_sg] = sg_w(sg, k, multiplier, rho);
toc1=toc;
N_SG = size(w_sg,1);%or N = size(X_set,1), N is the number of input obervations
fprintf('elapsed time of computing w_sg for d=%1.0f, eta=%1.0f and N=%1.0f is: %.8f seconds. \n',d,eta,N_SG,toc1')

% compute y_pred
m = 1e3;
xnew = x_left + rand(m,d).*(x_right-x_left);%m inputs in d dimension, m*d matrix
y_true = y_fun(xnew);
tic;
y_pred_sg = compute_pred(xnew, w_sg, X_set, k, rho);
toc2 = toc;
mse_sg = mean((y_true-y_pred_sg).^2);
fprintf('elapsed time of predicting y_pred_sg for d=%1.0f, eta=%1.0f and N=%1.0f, mse=%0.8f is: %.8f seconds. \n',d,eta,N_SG, mse_sg, toc1+toc2')

% plot
figure;
scatter3(xnew(:,1),xnew(:,2),y_true(:,1))
hold on;
scatter3(xnew(:,1),xnew(:,2),y_pred_sg(:,1))
xlabel('X1')
ylabel('X2')
zlabel('Griewank function')
