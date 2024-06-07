% Define a test matrix
A = [2 1; 1 3];

% Test 1
% Define the number of critical eigenvalues to extract
m = 1;

% Call the function
try
  [critical_left_eigenvectors, critical_right_eigenvectors,critical_eigenvalues] = compute_eigenvectors(A, m);
  disp('Test Passed: Function executed without errors.');
  display(critical_right_eigenvectors);
  display(critical_left_eigenvectors);
  display(critical_eigenvalues);
catch ME
  disp('Test Failed!');
  disp(ME.message);
end

% Test 2
% Define the number of critical eigenvalues to extract
m = 2;

% Call the function
try
  [critical_left_eigenvectors, critical_right_eigenvectors,critical_eigenvalues] = compute_eigenvectors(A, m);
  disp('Test Passed: Function executed without errors.');
  display(critical_right_eigenvectors);
  display(critical_left_eigenvectors);
  display(critical_eigenvalues);
catch ME
  disp('Test Failed!');
  disp(ME.message);
end
[Rm, D, Lm] = eig(A);
 display(Rm)
 display(Lm)
 display(D)

% Additional tests (optional)

% Test with non-square matrix
disp('--- Testing with non-square matrix ---');
B = [1 2 3; 4 5 6];
try
  [~, ~, ~] = getCriticalEigenvectors(B, m);
catch ME
  disp('Test Passed: Error thrown for non-square matrix.');
end

% Test with different number of critical eigenvalues
disp('--- Testing with different m ---');
m = 2;
try
  [critical_left_eigenvectors, critical_right_eigenvectors,critical_eigenvalues] = compute_eigenvectors(A, m);
  disp('Test Passed: Function extracts specified number of critical eigenvalues.');
  disp(['Extracted critical eigenvalues:', num2str(critical_eigenvalues)]);
catch ME
  disp('Test Failed!');
  disp(ME.message);
end
