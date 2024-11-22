%% Adding paths for accessing the functions

addpath(fullfile(pwd, 'Matrix_A'));

% Hermitian matrtix 4^4 lattice
% A = read_matrix('4x4x4x4b6.0000id3n1.mat'); % Read the input matrix from a file.


% Non-Hermitian matrtix 4^4 lattice
% A = read_matrix('periodic_L4_b3.55_k0.137n0_1.mat'); % Read the input matrix from a file.

% Non-Hermitian matrtix 8^4 lattice
% A = read_matrix('periodic_L8_b3.55_k0.137n0_1.mat'); % Read the input matrix from a file.
filename = 'critical_eigenvalues_2000_8to4.mat';

N = size(A, 2); % Size of the matrix

% gamma5hat = [speye(6), zeros(6,6); zeros(6,6), -speye(6)];
% Gamma5 = kron(speye(N/12),gamma5hat);
% A = Gamma5*A;

largest_eigenvalue = eigs(A, 1, 'largestabs');
% d = eig(full(A)); % Example if Hw is the matrix

% Load data
data = load(filename);
D = data.D1;
d = diag(D);

% Find the indices of zero elements
zeroIndices = find(d == 0);

% Check if the vector contains a zero
if ~isempty(zeroIndices)
    disp('The vector contains zeros at the following positions:');
    disp(zeroIndices);
else
    disp('The vector does not contain any zeros.');
end

save('eigenvalues.mat', 'd', "largest_eigenvalue");
