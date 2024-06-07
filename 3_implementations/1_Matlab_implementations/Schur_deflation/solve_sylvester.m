function Y = solve_sylvester(Tm, Hk, f_Tm, f_Hk)
    % Function to solve the Sylvester equation for Y
    n = size(Hk, 1);
    I = eye(n);
    A = kron(Tm, I) - kron(I, Hk);
    B = kron(f_Tm, ones(n,1)) - kron(ones(n,1), f_Hk);
    Y = A \ B;
end