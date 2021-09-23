%% Example with compute_basis.m, return band matirx A and Phi

x_1 = 0;%left-endpoint
x_n = 1;%right-endpoint
n = 11;%the number of observations
k = 5;%smoothness of Matern kernel, denotes(k-2)/2 Matern

rho = sqrt(k-2);%lengthscale of Matern kernel
x_input = linspace(x_1, x_n, n);%input observation row vector
x_Phi = linspace(x_1, x_n, 1e3)';

y_fun = @(x) sum( (x) .^2, 2)/4000 - prod(cos( (x) ./sqrt(1:size(x,2))), 2) + 1; %Griewank function. true function
y_true = y_fun(x_input');%true response

tic;
[A, Phi, Phi_fun] = compute_basis(x_input, k, 'null', rho, x_Phi);
toc_basis = toc;

fprintf('elapsed time of computing basis in 1D for %1.0f input is: %.8f seconds. \n',n,toc_basis')
plot(x_Phi', abs(Phi_fun),'LineWidth', 1);
title(sprintf('Kernel Packet Basis of degree k=%1.0f, Matern-%1.0f/2', k, k-2),'FontSize', 20)
%% Example with compute_post.m, return posterior mean and covariance

x_1 = 0;%left-endpoint
x_n = 1;%right-endpoint
n = 1e2+1;%the number of observations
k = 5;%smoothness of Matern kernel, denotes(k-2)/2 Matern

rho = sqrt(k-2)+1;%lengthscale of Matern kernel
x_input = linspace(x_1, x_n, n);%input observation row vector

y_fun = @(x) sum( (x) .^2, 2)/4000 - prod(cos( (x) ./sqrt(1:size(x,2))), 2) + 1; %Griewank function
y_input = y_fun(x_input');%true response

m = 1e3;
xnew = sort(rand(m,1));
y_true = y_fun(xnew);

tic;
[mean_new, cov_new] = compute_post(xnew, x_input, y_input, k, rho);
toc_post = toc;
mse_1d = mean((y_true-mean_new).^2);
fprintf('elapsed time of computing posteriorer in 1D for %1.0f input, mse=%.8f is: %.8f seconds. \n',n,mse_1d,toc_post')

figure;
plot(xnew, y_true);
hold on;
plot(xnew, mean_new);
title('ytrue versus ypred of degree k=5 in one dimension','FontSize', 20)
legend('$y_{true}$','$y_{pred}$','fontsize',14,'interpreter','latex')