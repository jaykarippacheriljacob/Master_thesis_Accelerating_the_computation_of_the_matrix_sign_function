function rho_vec = evalnodal( tt, active_nodes, subdiag )
    %% Evaluatioin of the node
    % Input: 
    %      tt           - quadrature points
    %      active_nodes - Interpolation nodes currently active in f
    %      subdiag      - subdiagonal elements of Hessenberg matrix
    % Output: 
    %      rho_vec      - nodal evaluation
    
    % Find rho_v - nodal evaluation
    % rho_v => prod(subdiag) ./ ( (tt-active_nodes(1))*(tt-active_nodes(2))*... )
    rho_vec = 0 * tt + 1;
    for j = 1:length(active_nodes)
        rho_vec = rho_vec * subdiag(j) ./ (tt - active_nodes(j));
    end
end