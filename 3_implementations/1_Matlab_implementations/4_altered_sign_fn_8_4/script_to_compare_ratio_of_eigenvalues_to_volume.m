% Adding paths for accessing the functions
addpath(fullfile(pwd, 'Matrix_A'));

% Hermitian matrix 4^4 lattice
% A = read_matrix('4x4x4x4b6.0000id3n1.mat'); % Read the input matrix from a file.

% Non-Hermitian matrix 4^4 lattice
% A = read_matrix('periodic_L4_b3.55_k0.137n0_1.mat'); % Read the input matrix from a file.

% Non-Hermitian matrix 8^4 lattice
A = read_matrix('periodic_L8_b3.55_k0.137n0_1.mat'); % Read the input matrix from a file.
filename = 'critical_eigenvalues_2000_8to4.mat';
% Load data
data = load(filename);
D = data.D1;

N = size(A, 2); % Size of the matrix

gamma5hat = [speye(6), zeros(6,6); zeros(6,6), -speye(6)];
Gamma5 = kron(speye(N/12),gamma5hat);
A = Gamma5*A;

% Values of m to test
m_values = [2, 4, 8, 16, 32, 64, 128];

% Store results
ratios = zeros(size(m_values));

% Main loop for each m
for i = 1:length(m_values)
    m = m_values(i); % Current number of eigenvalues
    
    % Find the m eigenvalues with smallest absolute values
    % eigenvalues = eigs(A, m, 'smallestabs');
    eigenvalues = diag(D(1:m, 1:m));
    
    % Find the largest eigenvalue among those found
    largest_selected = max(abs(eigenvalues));
    
    % Find the largest real eigenvalue of A
    largest_real = abs(eigs(A, 1, 'largestabs'));
    
    % Compute the ratio
    ratios(i) = largest_selected / largest_real;
end

% Display the results
disp('Ratios for each m:');
disp(table(m_values', ratios', 'VariableNames', {'m', 'Ratio'}));

% Save the results to a LaTeX table
output_filename = 'eigenvalue_ratios.tex';
fid = fopen(output_filename, 'w');
fprintf(fid, '\\begin{tabular}{|c|c|}\n');
fprintf(fid, '\\hline\n');
fprintf(fid, 'm & $\\max |\\lambda_{\\text{defl}}| / \\max |\\lambda_{\\text{all}}|$ \\\\\n');
fprintf(fid, '\\hline\n');
for i = 1:length(m_values)
    fprintf(fid, '%d & %.6f \\\\\n', m_values(i), ratios(i));
end
fprintf(fid, '\\hline\n');
fprintf(fid, '\\end{tabular}\n');
fclose(fid);

% Plot the results
figure;
plot(m_values, ratios, '-o');
xlabel('m (Number of Eigenvalues)');
ylabel('Ratio');
title('Ratio of Largest Selected Eigenvalue to Largest Real Eigenvalue');
grid on;

disp(['Results saved to ', output_filename]);
