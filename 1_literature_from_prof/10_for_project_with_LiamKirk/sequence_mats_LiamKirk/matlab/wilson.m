function [D] = wilson(kappa,theta,dirac,L,T,bc,flavour,zf,ds)
%
%   Wilson-Dirac operator.
%   Builds the sparse matrix that represents the Wilson-Dirac operator.
%
%   [D] = WILSON(KAPPA,THETA,DIRAC,L,T,BC,FLAVOUR,ZF)
%
%   gives the free Wilson-Dirac operator with mass M=1/(2*KAPPA)-DIRAC*R,
%   periodicity in space direction with phases EXP(I*THETA./L), FLAVOUR
%   type of quark, in DIRAC dimensions, spatial extension L=[L1 L2 L3], 
%   temporal extension T, temporal boundary conditions
%   BC (0: (anti-)periodic, 1: SF) and ZF (for twisted bc).
%   The matrix D has size NxN, where for BC=0
%   
%       N = L^(DIRAC-1)*T*DIRAC
%
%   and for BC>=1 
%
%       N = L^(DIRAC-1)*(T-1)*DIRAC.
%
%   [D] = WILSON(M,THETA,FLAVOUR)
%
%   assumes a lattice structure LAT (with gauge filed) in variable FLAVOUR
%   which has been imported with IMPORT_LATTICE. The matrix D has size NxN
%   with
%
%       N = L^(DIRAC-1)*T*LAT.DIRAC*LAT.COLOUR
% 
%   or
% 
%       N = L^(DIRAC-1)*(T-1)*LAT.DIRAC*LAT.COLOUR
% 
%   Bjoern Leder, 2008
%
% See also: PROJECTORS, IMPORT_LATTICE, SPARSE, CELL, KRON

if nargin==3
    % previously imported lattice with gauge field in variable dirac
    lat=dirac;
    colour=lat.colour;
    flavour=lat.flavour;
    dirac=lat.dirac;
    bc=lat.bc;
    L=lat.L;
    T=lat.T;
    zf=lat.zf;
    ds=lat.ds;
elseif nargin<9
    disp('Too few arguments. Usage:')
    help(mfilename)
    return
else
    % free Wilson-Dirac operator
    if length(L)~=dirac-1, error('Dimension mismatch.'); end
    lat=[];
    colour=1;
end
if bc==0
    T=T+1;
end
% constants
d=dirac-1;
ndc=dirac*colour;
Vs=prod(L);
VS=Vs*(T-2)+1;
V=Vs*(T-1);
if bc, VB=Vs; else VB=0; end

% projectors
la=[1 exp(1i*theta./L(:)')];
P=projectors(dirac,la);
G=gamma_matrices_4d();

G2=G{2};
G{2}=G{3};
G{3}=G2;

% diagonal part
m0=1/(2*kappa)-dirac;
M=(dirac+m0)*ones(ndc,1);
Mb=(zf+ds*(dirac-1)+m0)*ones(ndc,1);
adiag=(0:ndc-1)';

% nearest neighbour
NB=double(nextnb(d,T-1,L,bc>0,2));

a=cell(V*(2*dirac+1),1);
b=cell(size(a));
v=cell(size(a));

% k=-1;
% p=1;

for nx=1:Vs*(T-1)
    nnx=ndc*(nx-1)+1;
    if ~isempty(lat), unx=dirac*(VB+nx-1)+1; end
    % x=y
    % D(nnx:nnx+d,nnx:nnx+d) = ID;
    nxi=(2*dirac+1)*(nx-1)+1;
    if ~isempty(lat) && isfield(lat,'sw'),
        [a{nxi},b{nxi},v{nxi}]=find(lat.csw*lat.sw{nx}+spdiags(M,0,ndc,ndc));
        a{nxi}=nnx+a{nxi}-1;
        b{nxi}=nnx+b{nxi}-1;
    else
        a{nxi}=nnx+adiag;
        b{nxi}=a{nxi};
        v{nxi}=M;
    end
    % nearest neighbour
    for mu=1:dirac
        % + mu-direction
        if mu>1 || nx<VS || bc==0,
            % D(nnx:nnx+d,NB(nx,mu):NB(nx,mu)+d) = -P{mu} U(x,mu);
            K=-P{mu};
            if mu>1 && (nx>=VS || nx<=Vs) && bc>1
                %K=K+0.5*(ds-1)*G{mu};
                K=K*ds;
            end
            if ~isempty(lat)
                U=double(lat.U{unx+mu-1});
                if mu==1 && isfield(lat,'mu'),
                   U=exp(lat.mu)*U;
                end
                K=kron(K,U);
            end
            inx=nxi+mu;
            [a{inx},b{inx},v{inx}]=find(K);
            a{inx}=nnx+a{inx}-1;
            b{inx}=ndc*(NB(nx,mu)-1)+b{inx};
        end
        % - mu-direction
        if mu>1 || nx>Vs || bc==0,
            % D(nnx:nnx+d,NB(nx,mu+dirac):NB(nx,mu+dirac)+d) = -P{mu+dirac} U(x-mu,mu)^-1;
            K=-P{mu+dirac};
            if mu>1 && (nx>=VS || nx<=Vs) && bc>1
                %K=K-0.5*(ds-1)*G{mu};
                K=K*ds;
            end
            if ~isempty(lat)
                U=double(lat.U{dirac*(VB+NB(nx,mu+dirac)-1)+mu}');
                if mu==1 && isfield(lat,'mu'),
                   U=exp(-lat.mu)*U;
                end
                K=kron(K,U);
            end
            inx=nxi+mu+dirac;
            [a{inx},b{inx},v{inx}]=find(K);
            a{inx}=nnx+a{inx}-1;
            b{inx}=ndc*(NB(nx,mu+dirac)-1)+b{inx};
        end
    end
    % boundary x_0=T-1
    if nx>=VS,
        if bc>1,
            if bc==2,
                K=(-1)^flavour * 1i*G{5}*P{1+dirac};
%                 K=-(p*P{1}-k/2*(speye(dirac)-i*G{1}*G{5}));
            elseif bc==3
                K=-speye(dirac);
            end
            if ~isempty(lat)
                K=kron(K,speye(colour));
            end
            [a{nxi},b{nxi},v{nxi}]=find(K+diag(Mb));
            a{nxi}=nnx+a{nxi}-1;
            b{nxi}=nnx+b{nxi}-1;
        end
    end
    % boundary x_0=1
    if nx<=Vs,
        if bc>1,
            if bc==2,
               K=(-1)^flavour * 1i*G{5}*P{1};
%                 K=-(p*P{dirac+1}-k/2*(speye(dirac)+i*G{1}*G{5}));
            elseif bc==3
                K=-speye(dirac);
            end
            if ~isempty(lat)
                K=kron(K,speye(colour));
            end
            [a{nxi},b{nxi},v{nxi}]=find(K+diag(Mb));
            a{nxi}=nnx+a{nxi}-1;
            b{nxi}=nnx+b{nxi}-1;
        end
    end
end
D=sparse(cat(1,a{:}),cat(1,b{:}),cat(1,v{:}),V*ndc,V*ndc);