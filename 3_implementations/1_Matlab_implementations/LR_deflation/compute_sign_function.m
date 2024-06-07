% function f_Tm = compute_sign_function(S0, m)
%     % Use the given diagonal matrix as the initial matrix S0
% 
%     % Iterate the recurrence relation
%     S = S0;
%     for i = 1:m
%         S = 0.5 * (S + inv(S));
%     end
% 
%     % Extract the upper triangular part of S_m
%     f_Tm = triu(S);
% end
function f_Tm = compute_sign_function(S0)
    % Use the given diagonal matrix as the initial matrix S0

    % Iterate the recurrence relation
    S = S0;
    tol = 1e-8; % tolerance for convergence
    iterations = 0;
    while true
        S1 = 0.5 * (S + (S \ eye(size(S))));
        if norm(S1 - S, 'fro') < tol
            break;
        else
            S = S1;
            iterations = iterations + 1;
        end
    end
    % Display the number of iterations
    display(iterations);

    % Extract the upper triangular part of S_m
    f_Tm = triu(S);
end
