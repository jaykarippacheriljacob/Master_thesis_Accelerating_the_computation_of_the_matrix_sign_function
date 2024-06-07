function Y = sylvester_equation(A, Sm, Tm, Hk, Vk, f_Tm, f_Hk)
  % This function implements the Sylvester equation.
  % Tm, Hk are square matrices, Y and X are vectors, and f is a function.
  
  % Set tolerance level
  tol = 1e-8;

  % Setting X = Sm'*A*Vk
  X = Sm'*(A*Vk);

  % C = f_Tm * X -  X * f_Hk
  C = (f_Tm * X) -  (X * f_Hk);
  Tm_inv = inv(Tm);

  Y0 = zeros(size(Tm, 1), size(Hk, 1));
  iter = 0;
  while true
      
    %Y = Tm_inv * (C + Y0 * Hk);
    Y = Tm_inv * (C + (Y0 * Hk));
    rel_err = (norm(Y - Y0) / norm(Y));
    % disp(rel_err);
    if (rel_err < tol)
        break;
    elseif iter > 100
        % disp(rel_err);
        disp('No.of iteration > 100 exiting');
        break;
    else
        Y0 = Y;
        iter = iter + 1;
    end
  end
end
