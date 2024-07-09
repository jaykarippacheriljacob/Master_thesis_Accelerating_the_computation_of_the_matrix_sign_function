function A = generateSparseMatrix(N)
    % Check if N is odd
    if mod(N, 2) == 0
        diag_values = [1:N/2, -(N/2):-1];
    else
        mid = (N+1)/2;
        diag_values = [1:mid-1, 1, -mid+1:-1];
    end

    % Create the sparse matrix with the central diagonal
    A = spdiags(diag_values', 0, N, N);
end
