% Test script for compute_sign_function_diag function with sparse input

matrix_size = 5; % Size of the matrix
% Generate random diagonal elements
diagonal_elements = rand(matrix_size, 1);

% Create the diagonal matrix
input_matrix = diag(diagonal_elements);

% Display the input matrix
disp('Input Sparse Matrix:');
disp(full(input_matrix)); % Displaying the full matrix for clarity

% Compute sign function applied to the diagonal elements
result_matrix = compute_sign_function_diag(input_matrix);

% Display the resultant sign function matrix
disp('Resultant Sign Function Matrix:');
disp(full(result_matrix)); % Displaying the full matrix for clarity

% Verify the result
expected_result_matrix = sign(input_matrix);

% Display the expected sign function matrix
disp('Expected Sign Function Matrix (computed using sign function on diagonal elements):');
disp(full(expected_result_matrix));

% Check if the computed result matches the expected result
if isequal(full(result_matrix), full(expected_result_matrix))
    disp('Test Passed: Result matches expected sign function matrix.');
else
    disp('Test Failed: Result does not match expected sign function matrix.');
end
