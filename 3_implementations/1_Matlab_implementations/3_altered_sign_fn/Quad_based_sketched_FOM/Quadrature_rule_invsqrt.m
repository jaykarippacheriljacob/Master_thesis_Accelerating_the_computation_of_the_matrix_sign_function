function [c, z, l, h] = Quadrature_rule_invsqrt(V, AV, v, tol)
    % Input: 
    %      V - N x m, m arnoldi vectors
    %      AV - N x m,  Matrix
    %      v - N x 1, Sketched vector(Sb) of b
    %      tol - tollerance for the error computed to be.
    % Output: 
    %      c - quadrature nodes
    %      z - quadrature weights
    %      l - no.of quadrature nodes
    %      h - hm calculated based on the quadrature rule
    
    % Step 1: Set l_ := 8 and l := round(sqrt(2)*l_)
    l1 = 8;
    l2 = floor(sqrt(2) * l1);
    % Step 2: Set accurate := false and refined := false
    quad_err = inf; % Initial error set to infinity for the while loop condition
    first = true; % Flag to indicate the first iteration

    while quad_err > tol && l2 < 1000
        % Step 3: Choose sets (ti~, wi) i = 1,...,l~ and (ti, wi) i =
        %         1,...,l of quadrature nodes / weights

        c1 = pi / l1 * ones(1, l1); % quadrature nodes for l~
        c2 = pi / l2 * ones(1, l2); % quadrature nodes for l
        z1 = cos((2 * (1:l1) - 1)/(2 * l1) * pi); % quadrature weights for l~
        z2 = cos((2 * (1:l2) - 1)/(2 * l2) * pi); % quadrature weights for l
        
        % Step 4: Compute hm~ and hm by quadrature of order l~ and l
        %         respectively
        if first
            h1 = 0;
            for j = 1:l1
                h1 = h1 + c1(j) * ((V' * (-(1 - z1(j)) * V - (1 + z1(j)) * AV))\(V' * v));
            end
            h1 = -2/pi * h1;
            first = false;
        end
    
        h2 = 0;
        for j = 1:l2
            h2 = h2 + c2(j) * ((V' * (-(1 - z2(j)) * V - (1 + z2(j)) * AV))\(V' * v));
        end
        h2 = -2/pi*h2;
        
        % Step 5: Check the quadrature error w.r.t tot the set tolerance.
        quad_err = norm(h1 - h2);
        if quad_err <= tol || floor(sqrt(2) * l2) >= 1000
            l = l2;
            c = c2;
            z = z2;
            h = h2;
            disp(['Number of quadrature nodes: ', num2str(l)]);
            break;
        else
            l1 = l2;
            l2 = floor(sqrt(2) * l2);
            h1 = h2;
        end
    end
end