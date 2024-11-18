%% Set matrix size and load data
% A = read_matrix('../Matrix_A/periodic_L4_b3.55_k0.137n0_1.mat');
A = read_matrix('../Matrix_A/periodic_L8_b3.55_k0.137n0_1.mat');
N = size(A, 2);

gamma5hat = [speye(6), zeros(6,6); zeros(6,6), -speye(6)];
Gamma5 = kron(speye(N/12),gamma5hat);
A = Gamma5*A;

m = 2000; %no.of critical eigenvalues

% Setting up of the gmres solver
tol = 1e-6;
maxit = 200;

% Generate the Ritz values (example for m=10, adjust as needed)
m1 = 100; % Krylov subspace dimension

%% Start the timing
start = cputime;
theta1 = ritz_value(A, m1);
precond1 = @(x) eval_precon_poly_inv(A, x, theta1, m1);

Afun1 = @(x) gmres(A, x, [], tol, maxit, precond1);
% Afun = @(x) gmres(A, x, [], tol, maxit);

[Vm, D1] = eigs(Afun1, N, m, 'smallestabs'); % Computation time for 2000 eigenvalues of 8^4 lattice: 19480.8594 seconds 

% End the timing
finish = cputime;
disp(['Computation time: ', num2str(finish - start), ' seconds']);

% Save eigenvalues and eigenvectors to file
% save('eigen_values.mat', 'D1');
%% Start the timing
start = cputime;
% theta2 = ritz_value(A', m1);
precond2 = @(x) eval_precon_poly_inv(A', x, theta1, m1);

Afun2 = @(x) gmres(A', x, [], tol, maxit, precond1);
% Afun = @(x) gmres(A, x, [], tol, maxit);

[Wm, D2] = eigs(Afun2, N, m+1, 'smallestabs'); % Computation time for 2000 eigenvalues of 8^4 lattice: 62166.7656 seconds

% End the timing
finish = cputime;
disp(['Computation time: ', num2str(finish - start), ' seconds']);

%% Evaluation of the eigen values to see if the left and right
% eigenvectors are ordered in pairs
tol = 1e-8; % Tolerance for comparing floating point values
for i = 1:size(D1, 1)
    
    % Check if D2m(i,i) is the conjugate of D1m(i,i)
    if abs(conj(D1(i,i)) - D2(i,i)) > tol
        %disp(['Conjugate eigenvalue correction needed for eigenvalue ' num2str(i)]);
        
        % Find the conjugate of D1m(i,i) from the elements inside D2m
        for j = i+1:size(D2, 1)
            if abs(conj(D1(i, i)) - D2(j,j)) < tol
                
                % Swap places
                %D2(i, i) = D2(i, i) + D2(j, j);
                %D2(j, j) = D2(i, i) - D2(j, j);
                %D2(i, i) = D2(i, i) - D2(j, j);
                temp = D2(i, i);
                D2(i, i) = D2(j, j);
                D2(j, j) = temp;

                
                % Swap places for corresponding eigenvectors
                temp = Wm(:, j);
                Wm(:, j) = Wm(:, i);
                Wm(:, i) = temp;
                %Wm(:, [i, j]) = Wm(:, [j, i]);
                break;
            end
        end
    end
end

% Remove the last column of Wm
Wm = Wm(:, 1:m);

% D = diag(D1);

% Save eigenvalues and eigenvectors to file
save('critical_eigenvalues_2000_8to4.mat', 'Vm', 'Wm', 'D1');

% % Output results
% finish = cputime;
% disp('Eigenvalues:');
% disp(diag(D));
% disp(['Computation time: ', num2str(finish - start), ' seconds']);

% start = cputime;
% [Vm, D1] = eigs(A, m, 'smallestabs');
% finish = cputime;
% disp('Eigenvalues:');
% disp(diag(D1));
% disp(['Computation time: ', num2str(finish - start), ' seconds']);
