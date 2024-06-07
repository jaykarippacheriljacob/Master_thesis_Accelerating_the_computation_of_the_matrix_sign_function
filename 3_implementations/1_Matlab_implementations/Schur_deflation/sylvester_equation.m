function Y = sylvester_equation(Tm, Hk, X, f_Tm, f_Hk)
  % This function implements the Sylvester equation.
  % Tm, Hk are square matrices, Y and X are vectors, and f is a function.
  tol = 1e-8;
  Tm_inv = Tm\ eye(size(Tm));
  Y0 = eye(size(X));
  display(Y0)
  while true
      %Y = Tm_inv * (f_Tm * X -  f_Hk * X + Y0 * Hk);
    Y = Tm_inv * (f_Tm * X -  f_Hk * X + Hk * Y0);
    if (norm(Y - Y0) < tol)
        break;
    else
        Y0 = Y;
    end
  end
end
