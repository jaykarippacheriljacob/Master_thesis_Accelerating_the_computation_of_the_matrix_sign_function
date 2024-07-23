clear;
clc;
close all;

% Define test parameters
rng(2130);
% N = 50; % Size of the matrix
% A = rand(N) + 10.0 * eye(N); % Generate a random N x N matrix
% A = gallery('poisson', 50);
%d = eig(A); % Compute the eigen values of the generated matrix
A = read_matrix('4x4x4x4b6.0000id3n1.mat'); % Read the input matrix from a file.
N = size(A, 2); % Size of the matrix
gamma5hat = [speye(6), zeros(6,6); zeros(6,6), -speye(6)];
Gamma5 = kron(speye(N/12),gamma5hat);
A = Gamma5*A;

% A = A - 0.8 *speye(N);

% d = eig(full(A)); % Compute the eigen values of the generated matrix
b = randn(N, 1); % Generate a random N x 1 vector
% b = ones(N, 1);
s = 801; % Sketch matrix row dimension
trunc = 3; % Truncate orthogonalization to the last '3' vector

% plot(real(d), imag(d), '*'); % plot the real vs imaginary part of the eigen values

A_sqr = full(A) * full(A);
% d = eig(full(A_sqr)); % Compute the eigen values of the generated matrix
% figure;
% plot(real(d), imag(d), '*'); % plot the real vs imaginary part of the eigen values
Ab = full(A) * b;

start = cputime;

% Compute f(A)x directly using the sign function
exact_result = (A*(inv(sqrtm(full(A * A)))))*b;
% exact_result = sqrtm(inv(full(A)))*b;
% Save the value to exact_result.mat file
save('exact_result.mat', 'exact_result');
% Load the value from the file
% loadedData = load('exact_result.mat', 'exact_result');
% exact_result = loadedData.exact_result;  % Extract the value from the structure

finish = cputime;
disp(['Time taken without Quadrature based sketched FOM scheme = ', num2str(finish - start), ' s']);
    
% Initialize the figure and set up the plot
figure;
hold on;
ylabel('log(Relative Error)');
xlabel('m');
title('log(Relative Error) vs. m');
grid on;

for m = 5:10:800
    start = cputime;
    
    % Call the Quadrature based restarted arnoldi function
    % Quadrature_based_sketched_FOM_approximation = Quad_based_sketched_trunc_FOM(A_sqr, Ab, m, s, trunc);
    Quadrature_based_sketched_FOM_approximation = Quad_based_sketched_trunc_FOM(A, b, m, s, trunc);
    
    finish = cputime;
    disp(['Time taken by Quadrature based sketched FOM scheme with s "', num2str(s), '" & m "', num2str(m), '" = ', num2str(finish - start), ' s']);
     
    % Compute relative error
    rel_err = norm(exact_result - Quadrature_based_sketched_FOM_approximation) / norm(exact_result);

    % Plot the point
    plot(m, log10(rel_err), 'bo');
    % plot(m, rel_err, 'bo');
end

hold off;