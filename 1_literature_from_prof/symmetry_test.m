load("4x4x4x4b6.0000id3n1.mat");
n = size(D,1);

%gamma5 = [speye(2), zeros(2,2); zeros(2,2), -speye(2)]

gamma5hat = [speye(6), zeros(6,6); zeros(6,6), -speye(6)];
Gamma5 = kron(speye(n/12),gamma5hat);

Q = Gamma5*D;

norm(Q-Q','fro')
