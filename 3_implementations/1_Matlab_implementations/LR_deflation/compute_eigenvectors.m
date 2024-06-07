function [critical_left_eigenvectors, critical_right_eigenvectors, critical_eigenvalues] = compute_eigenvectors(A, m)

  % Compute eigenvalues and eigenvectors
  %[Vm, D, Wm] = eig(A);

  % A * Vm = Vm * Dm (eigen values) -> A * vi = vi * di
  % Wm' * A = Dm * Wm' -> A' * Wm = Wm * Dm' -> A' * wi = wi * conj(di)

  % (v1)0.4376 + 0.0000i, (v2)-0.4524 + 0.0000i, (v3)-0.4655 - 0.2961i
  % (w1)0.4376 + 0.0000i, (w2)-0.4524 + 0.0000i, (w3)-0.4655 - 0.2961i

   % (v1)0.3285 + 0.0000i, (v2)-0.5582 - 0.1525i, (v3)-0.5582 + 0.1525i
   % (w1)0.3285 + 0.0000i, (w2)-0.5582 - 0.1525i, (w3)-0.5582 + 0.1525i

  [Vm, D1] = eigs(A,m,'smallestabs');
  display(diag(D1))

  [Wm, D2] = eigs(A',m,'smallestabs');
  display(diag(D2))

  % Task: To keep Vm as it is but to reorder Wm to match Vm, via D2'

  % A = X * D * inv(X) -> A.' = inv(X).' * D.' * X.' -> inv(X).' * X.' = I
  % because inv(X).' = inv(X.')
  % the eigen values does not change but the eigen vectors iinterchange for transpose of a matrix 

  % Ensure matrix is square (error handling)
  if size(A, 1) ~= size(A, 2)
      error('Input matrix A must be square.');
  end

  % Sort eigenvalues (absolute values) in ascending order
  %[eigenvalues_abs, sorted_indices] = sort(abs(diag(D)), 'ascend');

  % Extract critical eigenvalues and corresponding eigenvectors
  %critical_eigenvalues = diag(eigenvalues_abs(1:m)); % Make critical eigenvalues a diagonal matrix
  %critical_right_eigenvectors = V(:, sorted_indices(1:m));
  %critical_left_eigenvectors = W(:, sorted_indices(1:m))';

end
