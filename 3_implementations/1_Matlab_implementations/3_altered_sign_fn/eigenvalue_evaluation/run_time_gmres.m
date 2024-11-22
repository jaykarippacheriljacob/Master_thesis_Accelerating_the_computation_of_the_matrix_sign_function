%% Initialize
A = read_matrix('../Matrix_A/periodic_L4_b3.55_k0.137n0_1.mat');
% A = read_matrix('../Matrix_A/periodic_L8_b3.55_k0.137n0_1.mat');
N = size(A, 2);

gamma5hat = [speye(6), zeros(6,6); zeros(6,6), -speye(6)];
Gamma5 = kron(speye(N/12), gamma5hat);
A = Gamma5 * A;

tol = 1e-6;
maxit = 200;

% Define m1 values
m1_values = [32];

% Fixed m value
m = 128;  

% Storage for results
run_times_precond = zeros(size(m1_values));  % Run times with preconditioner

%% Run GMRES without preconditioner for fixed m
start = cputime;
Afun_no_precond = @(x) gmres(A, x, [], tol, maxit);  % No preconditioner
% [~, D_no_precond] = eigs(Afun_no_precond, N, m, 'smallestabs');
finish = cputime;
run_times_no_precond = finish - start;
disp(['GMRES without preconditioner computation time: ', num2str(run_times_no_precond), ' seconds']);

%% Run for each m1 with preconditioner
for idx = 1:length(m1_values)
    m1 = m1_values(idx);
    disp(['Running for m1 = ', num2str(m1)]);
    
    % Precompute Ritz values for current m1
    start = cputime;
    theta1 = ritz_value(A, m1);  % Precompute Ritz values
    precond1 = @(x) eval_precon_poly_inv(A, x, theta1, m1);
    precond_setup_time = cputime - start;
    disp(['Preconditioning setup time for m1 = ', num2str(m1), ': ', num2str(precond_setup_time), ' seconds']);
    
    % Timing start
    start = cputime;
    
    % GMRES with preconditioner
    Afun1 = @(x) gmres(A, x, [], tol, maxit, precond1);
    [Vm, D1] = eigs(Afun1, N, m, 'smallestabs');
    
    % Timing end
    finish = cputime;
    run_times_precond(idx) = finish - start;
    disp(['Computation time for m1 = ', num2str(m1), ': ', num2str(run_times_precond(idx)), ' seconds']);
    
    % Save eigenvalues and eigenvectors for this m1
    save(['critical_eigenvalues_m1_', num2str(m1), '_4to4.mat'], 'Vm', 'D1');
end

%% Save runtime data
save('run_times_preconditioned_m1.mat', 'm1_values', 'run_times_precond');
save('run_time_no_preconditioner.mat', 'run_times_no_precond');

disp('All computations complete and results saved.');
