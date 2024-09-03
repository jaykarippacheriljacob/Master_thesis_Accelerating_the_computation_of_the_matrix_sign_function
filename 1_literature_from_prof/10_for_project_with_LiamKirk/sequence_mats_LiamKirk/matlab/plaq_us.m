function [Us] = plaq_us(lat,nx,mu,nu,nb)
%
%   Returns the links arround the plaquette.
%
%          Us_3
%    x+nu --<-- x+nu+mu
%      |         |
% Us_4 v         v Us_2
%      |         |
%      x  --<-- x+mu
%          Us_1
% 
%   Bjoern Leder, 2009
%
% See also: IMPORT_LATTICE

unxmu=4*(nb(nx,mu+1)-1)+1;
unxnu=4*(nb(nx,nu+1)-1)+1;
unx=4*(nx-1)+1;

Us{1}=lat.U{unx+mu};
Us{2}=lat.U{unxmu+nu};
Us{3}=lat.U{unxnu+mu};
Us{4}=lat.U{unx+nu};