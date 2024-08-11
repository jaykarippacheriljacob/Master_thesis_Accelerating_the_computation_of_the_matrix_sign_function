function p_A = eval_pre_condi_poly(A, b, theta, m)
    
    N = size(A, 1);
    p_A = zeros(N, 1); % Polynomial for preconditioning of matrix A
    I = eye(N); % Identity matrix of size N
    f = @(t) t.^(-1/2); % Define the function
    dd = zeros(m, m); % Initialize the divided difference table
    
    % Compute the first column of divided differences (function values at points theta)
    for i = 1:m
        dd(i, 1) = f(theta(i));
    end
    
    % Compute the divided differences
    for j = 2:m
        for i = 1:(m-j+1)
            dd(i, j) = (dd(i+1, j-1) - dd(i, j-1)) / (theta(i+j-1) - theta(i));
        end
    end
    
    % Extract the coefficients
    coefficients = dd(1, :);
    
    % Expressing the polynomial using newton basis
    for i=1:m
        if i > 1
                p_term = (A - theta(i-1) * I) * p_term;
        else
            % p_term = I;
            p_term = b;
        end
        p_A = p_A + coefficients(i) * p_term;
    end
end