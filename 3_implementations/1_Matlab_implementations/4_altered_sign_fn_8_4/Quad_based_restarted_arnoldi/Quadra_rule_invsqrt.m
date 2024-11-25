function [h2, N2] = Quadra_rule_invsqrt(A, active_nodes, subdiag, H, tol, N1)
    %% Computes the quadrature rules ql1(A, Vm, b) and ql2(A, Vm, b) 
    % Input: 
    %      A            - N x N, matrix
    %      active_nodes - m x 1, eigenvalues of hessenberg matrix
    %      subdiag      - subdiagonal elements of Hessenberg matrix
    %      H            - m x m, Hessenberg matrix
    %      tol          - tolerance for the error computed to be
    %      N1           - numbers of quadrature nodes
    % Output: 
    %     h2            - coefficients calculated based on the quadrature rule
    %     N2            - numbers of quadrature nodes
    
    % For f(z) = 1/sqrt(z), using Gauss-Jacobi quadrature
    % This implementation is for non - implicit non-Hermitian matrices.
    
    %% Step 1: Set l_ := 8 and l := round(sqrt(2)*l_)
    N1 = 8;
    N2 = floor(sqrt(2) * N1);
    m = size(H, 2);

    %% Step 2: Set accurate := false and refined := false
    quad_err = inf; % Initial error set to infinity for the while loop condition
    first = true; % Flag to indicate the first iteration
    % beta_transform = trace(A) / size(A, 1);
    beta_transform = 1;
    
    % Run the loop till convergence
    while quad_err > tol && N2 < 1000
        %% Step 3: Choose sets (ti~, wi) i = 1,...,l~ and (ti, wi) i =
        %         1,...,l of quadrature nodes / weights
        % and
        %% Step 4: Compute hm~ and hm by quadrature of order l~ and l
        %         respectively
        if first

            weights1 = pi / N1 * ones(1, N1);
            t1 = zeros(1, N1);
            
            for ii = 1:N1
                t1(ii) = cos((2 * ii - 1) / (2 * N1) * pi);
            end
            
            % Evaluate the reciprocal of the nodal polynomial at the
            % quadrature points
            tt = -beta_transform * (1 - t1) ./ (1 + t1);

            % Find rho_v - nodal evaluation
            rho_vec = 0 * tt + 1;
            for j = 1:length(active_nodes)
                rho_vec = rho_vec * subdiag(j) ./ (tt - active_nodes(j));
            end

            ee = zeros(m, 1);
            ee(1) = 1;
            h1 = zeros(size(ee));

            for j = 1:length(t1)
                h1 = h1 + weights1(j) * rho_vec(j) * ((-beta_transform * (1 - t1(j)) * eye(m) - H * (1 + t1(j))) \ ee);
            end
            h1 = -2 * sqrt(beta_transform) / pi * h1;

            first = false; % Flag to indicate the first iteration to be set false
        end

        weights2 = pi / N2 * ones(1, N2);
        t2 = zeros(1, N2);

        for ii = 1:N2
            t2(ii) = cos((2 * ii - 1) / (2 * N2) * pi);
        end

        % Evaluate the reciprocal of the nodal polynomial at the
        % quadrature points
        tt = -beta_transform * (1 - t2) ./ (1 + t2);

        % Find rho_v - nodal evaluation
        rho_vec = 0 * tt + 1;
        for j = 1:length(active_nodes)
            rho_vec = rho_vec * subdiag(j) ./ (tt - active_nodes(j));
        end

        ee = zeros(m, 1);
        ee(1) = 1;
        h2 = zeros(size(ee));

        for j = 1:length(t2)
            h2 = h2 + weights2(j) * rho_vec(j) * ((-beta_transform * (1 - t2(j)) * eye(m) - H * (1 + t2(j))) \ ee);
        end
        h2 = -2 * sqrt(beta_transform) / pi * h2;

        %% Step 5: Check the quadrature error w.r.t tot the set tolerance.
        quad_err = norm(h1 - h2) / norm(h1);
        % quad_err = norm(h1 - h2);
        if quad_err > tol
            N1 = N2;
            N2 = floor(sqrt(2) * N1);
            h1 = h2;
        else
            break;
        end
    end
end