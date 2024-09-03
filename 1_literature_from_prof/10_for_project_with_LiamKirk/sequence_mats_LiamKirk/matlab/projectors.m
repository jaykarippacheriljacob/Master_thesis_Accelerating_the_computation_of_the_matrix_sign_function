function [P] = projectors(dim,la)
%
%   Projectors P_\pm(\mu) = 1/2 * ( 1 \pm \gamma_\mu )
%
%   [P] = PROJECTORS(DIM,LA)
%
%   Works for DIM=2,4. Result is cell array of length DIM.
%   Find  P_-(0) in P{1}, P_-(1...DIM-1) in P{2:DIM}
%   and  P_+(0) in P{1+DIM}, P_+(1...DIM-1) in P{DIM+2:2*DIM}.
%
%   If LA given: 
%                   P_+(\mu) -> CONJ(LA) * P_+(\mu)
%                   P_-(\mu) ->      LA  * P_-(\mu)
%
%   Bjoern Leder, 2008
%
% See also: GAMMA_MATRICES_4D, GAMMA_MATRICES_2D, CELL, CONJ

if nargin<2
    la=ones(1,dim);
end
gamma=eval(['gamma_matrices_' num2str(dim) 'd']);
ID=speye(dim);

for i=1:dim
    % P_-
    P{i}=la(i)*1/2*(ID-gamma{i});
    % P_+
    P{i+dim}=conj(la(i))*1/2*(ID+gamma{i});
end