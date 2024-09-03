function [h2] = Quadrature_rule_invsqrt(A, active_nodes, subdiag, thick_replaced, H, m, tol, ell, k)
    % Input: %%%%% Need to be revised
    %      A - N x N, matrix
    %      active_nodes - Interpolation nodes currently active in f
    %      subdiag - subdiagonal elements of Hessenberg matrix
    %      thick_replaced - subspace of prev. cycle (thick restart)
    %      H - m x m, Hessenberg matrix
    %      m - no. of iterations for the kryl ov's subspace
    %      tol - tollerance for the error computed to be.
    %      ell - thick restart parameter
    %      k - current iteration number of the restarted arnoldi
    % Output: 
    %     h2 - coefficients calculated based on the quadrature rule
    
    % For f(z) = 1/sqrt(z), using Gauss-Jacobi quadrature
    % This implementation is for non - implicit non-Hermitian matrices.
    
    % Step 1: Set l_ := 8 and l := round(sqrt(2)*l_)
    N1 = 8;
    N2 = floor(sqrt(2) * N1);

    % Step 2: Set accurate := false and refined := false
    quad_err = inf; % Initial error set to infinity for the while loop condition
    first = true; % Flag to indicate the first iteration
    % beta_transform = trace(A) / size(A, 1);
    beta_transform = 1;
    
    % Run the loop till convergence
    while quad_err > tol && N2 < 1000
        % Step 3: Choose sets (ti~, wi) i = 1,...,l~ and (ti, wi) i =
        %         1,...,l of quadrature nodes / weights
        % and
        % Step 4: Compute hm~ and hm by quadrature of order l~ and l
        %         respectively
        if first
            % using Gauss-Jacobi quadrature
            weights1 = pi / N1 * ones(1, N1);
            t1 = zeros(1, N1);
            
            for ii = 1:N1
                t1(ii) = cos((2 * ii - 1) / (2 * N1) * pi);
            end
            
            % Evaluate the reciprocal of the nodal polynomial at the
            % quadrature points
            tt = -beta_transform * (1 - t1) ./ (1 + t1);

            % Find rho_vec - nodal evaluation
            % rho_v => prod(subdiag) ./ ( (tt-active_nodes(1))*(tt-active_nodes(2))*... )
            rho_vec = evalnodal(tt, active_nodes(1:end-length(thick_replaced{k-1})), subdiag(1:end-length(thick_replaced{k-1}))).';
            rho_vec_replaced = evalnodal(tt, thick_replaced{k-1}, subdiag(end-length(thick_replaced{k-1})+1:end)).';
            rho_vec = rho_vec .* rho_vec_replaced;

            ee = zeros(ell+m, 1);
            ee(ell+1) = 1;
            h1 = zeros(size(ee));

            for j = 1:length(t1)
                h1 = h1 + weights1(j) * rho_vec(j) * ((-beta_transform * (1 - t1(j)) * eye(ell+m) - H * (1 + t1(j))) \ ee);
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

        % Find rho_vec - nodal evaluation
        rho_vec2 = evalnodal(tt, active_nodes(1:end-length(thick_replaced{k-1})), subdiag(1:end-length(thick_replaced{k-1}))).';
        rho_vec_replaced2 = evalnodal(tt, thick_replaced{k-1}, subdiag(end-length(thick_replaced{k-1})+1:end)).';
        rho_vec2 = rho_vec2 .* rho_vec_replaced2;

        ee = zeros(ell+m, 1);
        ee(ell+1) = 1;
        h2 = zeros(size(ee));

        for j = 1:length(t2)
            h2 = h2 + weights2(j) * rho_vec2(j) * ((-beta_transform * (1 - t2(j)) * eye(ell+m) - H * (1 + t2(j))) \ ee);
        end
        h2 = -2 * sqrt(beta_transform) / pi * h2;

        % Step 5: Check the quadrature error w.r.t tot the set tolerance.
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