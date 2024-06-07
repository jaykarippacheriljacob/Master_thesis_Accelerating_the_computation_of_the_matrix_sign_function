# Master_thesis_Accelerating_the_computation_of_the_matrix_sign_function
Thesis title: Accelerating the computation of the matrix sign function 

Content: The matrx sign function arisies in computations in lattice QCD. We look at the computation of the action sign(Q)x of the sign function of the matrix Q on a vector x. In our application, Q is the symmetrized Wilson-Dirac operator. This is a Hermitian matrix if the chemical potential is 0, otherwise it is non-Hermitian. Actually, we will always consider the inverse square root function, since sign(Q)x = (Q^2)^{1/2}Qx

The Arnoldi Kroylv subspace approximation is the basis method to approximate sign(Q)x. >There are several ways to accelerate the convergence of this basic scheme: 
  1. Restarts (in the non-hermitian case). This avoids to have too many inne rproducts in the Arnoldi aorthogonalization.
  2. Deflation (explicit and implicit). This make the matrix better coinditioned and thus reduces the number of iterations. Explicit deflation use the snmallest left and right eigenvectors, implicit deflation is present in the thick restart approach of Eiermann and Güttel; see also the funm Matlab code
  3. Polynomial preconditioning. This also makes the matrix better conditioned and this reduces the number of iterations. We have a recent paper on this and numerical results fpr QCD on a parallel machine.
  4. Sketching. This is a randomized approach where we save orthogonalizations and sketch the Arnoldi matrix. The relevant paper is by Güttel and Schweitzer.

The purpose of the thesis is to consider the following combination of the above approaches:
  1. 2 + 1 (as is already done in funmat)
  2. 2 + 3 (building on existing work and code of Gustavo)
  3. 2 + 4 (this is new, but Stefan Güttel just gave a talk on it at a conference in Paris)

Tasks:
  1. Understand and describe the individual methods (1-4)
  2. Describe, formulate algorithmically and discuss the combined methods (2+1, 2+3, 2+4)
  3. Test the combined methods, both in Matlab on small configurations and in C on large configurations and in            parallel (using existing code)
