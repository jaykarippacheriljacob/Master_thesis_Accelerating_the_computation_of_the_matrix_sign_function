A = read_matrix('periodic_L4_b3.55_k0.137n0_1.mat'); % Read the input matrix from a file.

N = size(A, 2); % Size of the matrix
gamma5hat = [speye(6), zeros(6,6); zeros(6,6), -speye(6)];
Gamma5 = kron(speye(N/12),gamma5hat);
A = Gamma5*A;
% Calculate the conjugate transpose of A
A_conjugate_transpose = A';

% Check if A is equal to its conjugate transpose
if isequal(A, A_conjugate_transpose)
    disp('The matrix is Hermitian.');
else
    disp('The matrix is not Hermitian.');
end

% Calculate eigenvalues of matrix A and A^2
spec_A = eig(full(A));          % Eigenvalues of matrix A
spec_A_sqr = eig(full(A*A));    % Eigenvalues of matrix A^2

% Plot the eigenvalues of A on the complex plane
figure;
plot(real(spec_A), imag(spec_A), '*'); 
xlabel('Real Part');            % Label for the x-axis
ylabel('Imaginary Part');       % Label for the y-axis
title('Eigenvalues of Matrix A'); % Title for the plot
grid on;                        % Add grid for better visualization
legend('Eigenvalues of A');     % Add legend
axis equal;                     % Equal scaling on both axes for accurate representation

% Plot the eigenvalues of A^2 on the complex plane
figure;
plot(real(spec_A_sqr), imag(spec_A_sqr), 's'); 
xlabel('Real Part');            % Label for the x-axis
ylabel('Imaginary Part');       % Label for the y-axis
title('Eigenvalues of Matrix A^2'); % Title for the plot
grid on;                        % Add grid for better visualization
legend('Eigenvalues of A^2');   % Add legend
axis equal;                     % Equal scaling on both axes for accurate representation
