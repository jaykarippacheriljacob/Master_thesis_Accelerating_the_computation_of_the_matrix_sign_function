function [gamma] = gamma_matrices_4d(rep)
%
%   Definition of the gamma matrices in 4 dim
%
%   GAMMA_MATRICES_4D(REP)
%
%   Bjoern Leder, 2008
%
% See also: SPARSE, KRON

if nargin==0 || isempty(rep)
    rep=2;
end

% identity
ID = eye(2);

% sigma matrices
sig1=[0 1;1 0];
sig2=[0 -i;i 0];
sig3=[1 0;0 -1];

rep_str={'SF','DDHMC','QOP'};
nrep=length(rep_str);
for k=1:nrep,
   id=sprintf('dirac_operator:gamma_matrices_%s',rep_str{k});
   try 
       q=warning('query',id);
       s(k)=strcmp(q.state,'off');
   catch
       s(k)=0;
   end
end
id=sprintf('dirac_operator:gamma_matrices_%s',rep_str{rep});
rep_old=find(s==1);
if ~isempty(rep_old) && rep_old~=rep
   wstr=sprintf('Changed gamma-matrix convention from %s to %s',rep_str{rep_old},rep_str{rep});
   warning('dirac_operator:gamma_matrices',wstr);
   id_old=sprintf('dirac_operator:gamma_matrices-%s',rep_str{rep_old});
   warning('on',id_old);
else
   wstr=sprintf('Gamma-matrix convention: %s',rep_str{rep});
   warning(id,wstr);
end
warning('off',id)  ;

% reps
if rep==1,
    % \gamma_0 diagonal
    error('dirac-operator:gamma-matrices','Not yet implemented!');
elseif rep==2,
    % \gamma_5 diagonal
    gamma{1}=sparse(kron(sig1,-ID));
    gamma{2}=sparse(blkdiag(-i*sig1,i*sig1)*(-gamma{1}));
    gamma{3}=sparse(blkdiag(-i*sig2,i*sig2)*(-gamma{1}));
    gamma{4}=sparse(blkdiag(-i*sig3,i*sig3)*(-gamma{1}));
    gamma{5}=gamma{1}*gamma{2}*gamma{3}*gamma{4};
    gamma{6}=i/2*(gamma{1}*gamma{2}-gamma{2}*gamma{1});
    gamma{7}=i/2*(gamma{1}*gamma{3}-gamma{3}*gamma{1});
    gamma{8}=i/2*(gamma{1}*gamma{4}-gamma{4}*gamma{1});
    gamma{9}=i/2*(gamma{2}*gamma{3}-gamma{3}*gamma{2});
    gamma{10}=i/2*(gamma{2}*gamma{4}-gamma{4}*gamma{2});
    gamma{11}=i/2*(gamma{3}*gamma{4}-gamma{4}*gamma{3});
elseif rep==3,
    % \gamma_5 diagonal
    gamma{1}=sparse(kron(sig1,ID));
    gamma{2}=sparse(blkdiag(i*sig1,-i*sig1)*(gamma{1}));
    gamma{3}=sparse(blkdiag(-i*sig2,i*sig2)*(gamma{1}));
    gamma{4}=sparse(blkdiag(i*sig3,-i*sig3)*(gamma{1}));
    gamma{5}=gamma{1}*gamma{2}*gamma{3}*gamma{4};
    gamma{6}=i/2*(gamma{1}*gamma{2}-gamma{2}*gamma{1});
    gamma{7}=i/2*(gamma{1}*gamma{3}-gamma{3}*gamma{1});
    gamma{8}=i/2*(gamma{1}*gamma{4}-gamma{4}*gamma{1});
    gamma{9}=i/2*(gamma{2}*gamma{3}-gamma{3}*gamma{2});
    gamma{10}=i/2*(gamma{2}*gamma{4}-gamma{4}*gamma{2});
    gamma{11}=i/2*(gamma{3}*gamma{4}-gamma{4}*gamma{3});
end
