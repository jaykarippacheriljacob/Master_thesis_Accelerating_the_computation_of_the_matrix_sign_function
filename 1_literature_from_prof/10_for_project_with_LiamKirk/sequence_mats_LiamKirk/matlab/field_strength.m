function [F] = field_strength(lat)
% field_strength.m Computes the field strength lat.U
%
% The definition of DD-HMC-.../doc/dirac.ps is used.
%
% F: cell array of size 6xlat.volume
%    F_{\mu\nu}(ix) is found in
%    F({[0 1],[0 2],[0 3],[1 2],[1 3],[2 3]},ix) 
%
% Bjoern Leder, 2011
%

NB=nextnb(lat.dirac-1,lat.T,lat.L,(lat.bc>0),2);
planes={[0 1],[0 2],[0 3],[1 2],[1 3],[2 3]};
npls=length(planes);

X=zeros(3*npls,3*lat.volume);
F=mat2cell(X,3*ones(npls,1),3*ones(lat.volume,1));
%F=cell(npls,lat.volume);
% [F{:}]=deal(zeros(3,3));
for n=1:npls,
    mu=planes{n}(1);
    nu=planes{n}(2);
    
    for nx=1:lat.volume,
        Us=plaq_us(lat,nx,mu,nu,NB);
        F{n,nx}=F{n,nx}+Us{1}*Us{2}*Us{3}'*Us{4}';
        nnx=NB(nx,mu+1);
        F{n,nnx}=F{n,nnx}+Us{2}*Us{3}'*Us{4}'*Us{1};
        nnx=NB(nnx,nu+1);
        F{n,nnx}=F{n,nnx}+Us{3}'*Us{4}'*Us{1}*Us{2};
        nnx=NB(nx,nu+1);
        F{n,nnx}=F{n,nnx}+Us{4}'*Us{1}*Us{2}*Us{3}';     
    end
end

for nx=1:lat.volume,
    for n=1:npls,
        F{n,nx}=(F{n,nx}-F{n,nx}')/8;
    end
end
