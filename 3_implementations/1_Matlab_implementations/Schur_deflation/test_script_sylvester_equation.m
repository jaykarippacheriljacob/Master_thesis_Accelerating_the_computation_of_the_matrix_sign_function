% Define matrices and vectors for testing
Tm = [6, -1; -1, 2];
Hk = [1, 0; 8, -1];
X = [19; 5];

% Define functions f_Tm and f_Hk
f_Tm = eye(size(Tm));
f_Hk = eye(size(Hk));

% Call sylvester_equation function
Y = sylvester_equation(Tm, Hk, X, f_Tm, f_Hk);

% Display result
disp('Solution Y:');
disp(Y);