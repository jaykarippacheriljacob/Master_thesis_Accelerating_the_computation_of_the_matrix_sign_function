function x_perp = orthogonal_projection(S_m, x)
    % Compute orthogonal projection
    I = eye(size(S_m, 1)); % Identity matrix
    x_perp = (I - S_m * pinv(S_m) ) * x;
end
