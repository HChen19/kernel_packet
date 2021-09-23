%% Example with mle_1d.m
d = 1;%dimension of input
k = 5;%smoothness of Matern kernel, denotes(k-2)/2 Matern
rho_init = sqrt(k-2);%initial value of rho
X = linspace(0, 1, 1e3)';%input

y_fun = @(x) sum( (x) .^2, 2)/4000 - prod(cos( (x) ./sqrt(1:size(x,2))), 2) + 1; %Griewank function. true function
f = @(x)ones(size(x,1),1);

[theta_hat,L_hat,L_init] = mle_1d(rho_init, k, X, y_fun, f);
%% Example with sg_mle.m
d = 2;%dimension of input
eta = 9;%level of sparse grid construction, which satisfies eta>= d 
k = 5;%smoothness of Matern kernel, denotes(k-2)/2 Matern
rho_init = sqrt(k-2);%initial value of rho
f = @(x)ones(size(x,1),1);

y_fun = @(x) sum( (x) .^2, 2)/4000 - prod(cos( (x) ./sqrt(1:size(x,2))), 2) + 1; %Griewank function. true function

x_left = 0;
x_right = 1;

design_fun = @(x)hc(x(1),x(2),x(3));
%design_fun = @(x)pl(x(1),x(2),x(3));

[sg] = sgd(d, eta, design_fun, x_left, x_right,'start=max');
[sg_logdet] = sgd(d, eta, design_fun, x_left, x_right,'start=d');

% [splogdet_M] = splogdet(d, eta, k, sg_logdet, rho_init);
% euclid_dist = pdist2(sg.X_set, sg.X_set);
% M = matern_halfint(euclid_dist, (k-2)/2, 1, rho_init);
% logdet_M = log(det(M));
tic
[theta_hat,L_hat,L_init] = sg_mle(rho_init, k, sg, sg_logdet, y_fun, f);
toc