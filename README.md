# Kernel Packaet
Kernel Packet is a MATLAB package for fast computing in Gaussian Process Modeling. 

## File Description
### Functions:
* [General](https://github.com/HChen19/kernel_packet/tree/main/Functions/General): 
  * [bdlogdet](https://github.com/HChen19/kernel_packet/blob/main/Functions/General/bdlogdet.m) is a function file to compute log determinant of [band matrix](https://en.wikipedia.org/wiki/Band_matrix#:~:text=In%20mathematics%2C%20particularly%20matrix%20theory,more%20diagonals%20on%20either%20side.).
  * [matern_halfint](https://github.com/HChen19/kernel_packet/blob/main/Functions/General/matern_halfint.m) is a function file to compute [Mat&eacute;rn covariance matrix for half integer](https://en.wikipedia.org/wiki/Mat%C3%A9rn_covariance_function#Simplification_for_%CE%BD_half_integer).
  * [nsumk](https://github.com/HChen19/kernel_packet/blob/main/Functions/General/nsumk.m) is a function file to compute a partition matirx where n positive integers summing to k.
  * [splogdet](https://github.com/HChen19/kernel_packet/blob/main/Functions/General/splogdet.m) is a function file to compute log determinant of high dimensional covariance matrix under sparse grid design.
  * [hc](https://github.com/HChen19/kernel_packet/tree/main/Functions/General/DesignFunctions) in [DesignFunctions](https://github.com/HChen19/kernel_packet/tree/main/Functions/General/DesignFunctions) folder is a function file in to compute a vector with hyperbolic cross points (bisection).
  * [pl](https://github.com/HChen19/kernel_packet/blob/main/Functions/General/DesignFunctions/pl.m) in [DesignFunctions](https://github.com/HChen19/kernel_packet/tree/main/Functions/General/DesignFunctions) folder is a function file to returen a vector with Plumlee's design.
 

* [One-dimension]():
  * Functions:
    * [compute_basis](https://github.com/HChen19/kernel_packet/blob/main/compute_basis.m) is a function file to compute kernel packet (which also our basis function) in one dimensional [Kriging (GPR)](https://en.wikipedia.org/wiki/Kriging).
    * [compute_post](https://github.com/HChen19/compact_support/blob/main/compute_post.m) is a function file to to compute posterior mean and covariance in one dimensional Kriging.
    * [loglike_1d](https://github.com/HChen19/compact_support/blob/main/loglike_1d.m) is a function file to compute log likelihood of parameter 'rho' in one dimension.
    * [MLE_1d](https://github.com/HChen19/compact_support/blob/main/MLE_1d.m) is a function file to compute updated parameter theta_hat, updated log likelihood L_hat and initial loglikelihood L_init in one dimensional MLE.

  * Examples:
    * [Example1D](https://github.com/HChen19/compact_support/blob/main/Example1D.m) is an example to use [compute_basis](https://github.com/HChen19/compact_support/blob/main/compute_basis.m) and [compute_post](https://github.com/HChen19/compact_support/blob/main/compute_post.m) function files. 

* High deimensional:
  * Functions:
    * [unfold](https://github.com/HChen19/compact_support/blob/main/unfold.m) is a function file to unfold tensor into matrix.
    * [refold](https://github.com/HChen19/compact_support/blob/main/refold.m) is a function file to refold matrix into tensor.
    * [kron_mat_prod](https://github.com/HChen19/compact_support/blob/main/kron_mat_prod.m) is a function file to compute matrix-matrix multiplication between matrix kron(As{1}, As{2}, ..., As{N}) and matrix v without forming the full kronecker product by using unfold and refold functions.
    * [sg_w](https://github.com/HChen19/compact_support/blob/main/sg_w.m) is a function file to compute multiplication of inverse of covariance matrix and a random matrix in sparse grid design.
    * [lt_w](https://github.com/HChen19/compact_support/blob/main/lt_w.m) is a function file to compute multiplication of inverse of covariance matrix and a random matrix in lattice design. 
    * [compute_pred](https://github.com/HChen19/compact_support/blob/main/compute_pred.m) is a funciton file to compute prediction value of response (which is also the posterior mean) of high dimensional Kriging.
    * [loglike_sg](https://github.com/HChen19/compact_support/blob/main/loglike_sg.m) is a function file to compute log likelihood of parameter 'rho' in high dimensional sparse grid design.
    * [MLE_sg](https://github.com/HChen19/compact_support/blob/main/MLE_sg.m) is a function file to compute updated parameter theta_hat, updated log likelihood L_hat and initial loglikelihood L_init in high dimensional MLE under sparse grid design.
  
  * Examples:
    * [SGDesign](https://github.com/HChen19/compact_support/blob/main/SGDesign.m) is an example to use [sg_w](https://github.com/HChen19/compact_support/blob/main/sg_w.m) and [compute_pred](https://github.com/HChen19/compact_support/blob/main/compute_pred.m) in sparse grid design.
    * [LTDeisgn](https://github.com/HChen19/compact_support/blob/main/LTDesign.m) is an example to compute multiplication of inverse of covariance matrix and response in lattice design (or called full grid design), then compute prediction value of response based on it.
    * [ExampleMLE](https://github.com/HChen19/compact_support/blob/main/ExampleMLE.m) is an example to use [MLE_1d](https://github.com/HChen19/compact_support/blob/main/MLE_1d.m) and [MLE_sg](https://github.com/HChen19/compact_support/blob/main/MLE_sg.m) function files.

* Others:
  * [tensor_toolbox-master](https://github.com/HChen19/compact_support/tree/main/tensor_toolbox-master) is a MATLAB tensor toolbox by Brett W. Bader, Tamara G. Kolda and others, we use [sparse tensor](https://www.tensortoolbox.org/sptensor_doc.html) function 'sptensor' in our file [sg_w](https://github.com/HChen19/compact_support/blob/main/sg_w.m). For more information please see their official website [here](https://www.tensortoolbox.org/).
  * [threeDplot](https://github.com/HChen19/compact_support/blob/main/threeDplot.m) is a file to plot kernel packet (our proposed basis function) with its components in 3D, which can help you understant intuitively our algorithm.

## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

Please make sure to update tests as appropriate.

## License
This work is published under [MIT](https://choosealicense.com/licenses/mit/) License.
