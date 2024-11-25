function [Vm, Wm, D1] = compute_eigenvectors(A, m)
    % load_and_trim_eigenvectors - Load and trim eigenvector-related variables.
    %
    % Inputs:
    %    m - Number of columns to extract from Vm and Wm and size of submatrix for D1
    %
    % Outputs:
    %    Vm - First m columns of the eigenvector matrix
    %    Wm - First m columns of the weight matrix
    %    D1 - Top m x m submatrix of the diagonal matrix
    
    % Hardcoded filename
    filename = 'critical_eigenvalues_2000_8to4.mat';
    
    % Check if the file exists
    if ~isfile(filename)
        error('File not found: %s', filename);
    end
    
    % Load data
    data = load(filename);
    
    % Ensure all required variables exist
    if ~isfield(data, 'Vm') || ~isfield(data, 'Wm') || ~isfield(data, 'D1')
        error('File does not contain the required variables: Vm, Wm, D1');
    end
    
    % Extract variables
    Vm_full = data.Vm;
    Wm_full = data.Wm;
    D1_full = data.D1;
    
    % Validate dimensions
    if size(Vm_full, 2) < m || size(Wm_full, 2) < m || size(D1_full, 1) < m || size(D1_full, 2) < m
        error('Not enough data in Vm, Wm, or D1 for the specified m.');
    end
    
    % Trim matrices
    Vm = Vm_full(:, 1:m);
    Wm = Wm_full(:, 1:m);
    D1 = D1_full(1:m, 1:m);
end
