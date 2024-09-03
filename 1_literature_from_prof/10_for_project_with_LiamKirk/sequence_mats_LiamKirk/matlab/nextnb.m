% function [NB] = nextnb(d,T,L,BC,order,bc_dirs)
%
%  x = {x_0, x_1, ..., x_d}
%  
%  n_x = x_0 + x_1 L_0 + x_2 L_0 L_1 + ... + x_d L_0 ... L_d-1
%
% Constructs the nearest-neighbour matrix for d+1 dimensional
% lattice of size V = T*L^d
%
% NB(nx,i),   i = 1..2(d+1) : x_mu+1 -> i=mu, x_mu-1 -> i=mu+d+1
%
% order == 1 -> L_0 = T  (default)
% order ~= 1 -> L_d = T
%
function [NB] = nextnb(d,T,L,BC,order,bc_dirs)

if nargin < 4,
    disp('Too few arguments! Usage:');
    help nextnb;
    return;
end
if nargin < 5,
    order=1;
end
if nargin < 6,
    bc_dirs=1;
end

if length(L)~=d,  error('Dimension mismatch!'); end

if BC<0 || BC>1, error('Unknown BC!'); end

V = T*prod(L);
NB = zeros(V,2*(d+1));
if BC,
    NB(V+1,:)=V+1;
end
if order==1
    Lx(1)=T;
    Lx(2:d+1)=L;
    mu_index=1:d+1;
else
    Lx(1:d)=L(end:-1:1);
    Lx(d+1)=T;
    mu_index=d+1:-1:1;
end

for n_x=0:V-1,
    r_x = n_x;
    Lmu = V;
    for mu=d:-1:0,
        % L^mu
        Lmu = floor(Lmu/Lx(mu+1));
        % x_mu
        x_mu = floor(r_x/Lmu);
        % n_x - x_mu*Lmu
        r_x = mod(r_x,Lmu);
        
        % forward
        x_muP=x_mu+1;
        while (x_muP >= Lx(mu+1)), x_muP=x_muP-Lx(mu+1); end
	    if BC && ~isempty(find(mu_index(mu+1)==bc_dirs, 1)) && x_muP==0,
            NB(n_x+1,mu_index(mu+1)) = V+1;                         % free BC
        else
            NB(n_x+1,mu_index(mu+1)) = 1 + n_x + (x_muP-x_mu)*Lmu;     % periodic BC
        end
        
        % backward
        x_muP=x_mu-1;
        while (x_muP < 0), x_muP=x_muP+Lx(mu+1); end
	    if BC && ~isempty(find(mu_index(mu+1)==bc_dirs, 1)) && x_mu==0,
            NB(n_x+1,mu_index(mu+1)+d+1) = V+1;                         % free BC
        else
            NB(n_x+1,mu_index(mu+1)+d+1) = 1 + n_x + (x_muP-x_mu)*Lmu;     % periodic BC
        end
    end
end

[v,type] = sint(V);
NB = eval([type '(NB)']);

% nnb = zeros(T,L);
% nx = zeros(T,L);
% dp=2;
% for x1=0:T-1,
%     for x2=0:L-1,
%         nx(x1+1,x2+1)=1+x1+x2*T*L^(dp-2);
%     end
% end
% nx
% for mu=[1 dp 1+d+1 dp+d+1],
%     mu
%     for x1=0:T-1,
%         for x2=0:L-1,
%             nnb(x1+1,x2+1)=NB(1+x1+x2*T*L^(dp-2),mu);
%         end
%     end
%     nnb
% end
