function [all_plaq ss_plaq st_plaq] = plaquette(lat)
%[all_plaq ss_plaq st_plaq] = plaquette(lat)
%
%   Compute spatial plaquette.
%
%    x+nu ---- x+nu+mu
%      |         |
%      |         |
%      x  ---- x+mu
%
%   [PLAQ SS_PLAQ ST_PLAQ] = PLAQUETTE(LAT)
% 
%   Bjoern Leder, 2008
%
% See also: IMPORT_LATTICE

nb=nextnb(lat.dirac-1,lat.T,lat.L,(lat.bc>0),2);
ss_plaq=0;
st_plaq=0;
all_plaq=0;
if lat.bc
    Vs=prod(lat.L);
else
    Vs=0;
end
for nx=1:Vs
    for mu=1:lat.dirac-1
        nu=0;
        st_plaq = st_plaq + real(trace(plaq(lat,nx,mu,nu,nb)));
    end
end
for nx=Vs+1:lat.volume
    for mu=2:lat.dirac-1
        for nu=1:mu-1
            p = real(trace(plaq(lat,nx,mu,nu,nb)));
            ss_plaq = ss_plaq + p;
            all_plaq = all_plaq + p;
        end
    end
    if nx<=lat.volume-Vs
        for mu=1:lat.dirac-1
            nu=0;
            p = real(trace(plaq(lat,nx,mu,nu,nb)));
            st_plaq = st_plaq + p;
            all_plaq = all_plaq + p;
        end
    end
end
ss_plaq=ss_plaq/3/(lat.volume-Vs);
st_plaq=st_plaq/3/(lat.volume-Vs);
all_plaq=all_plaq/6/(lat.volume-Vs);
%plaq=(ss_plaq+st_plaq)/2;
