function [P] = plaq(lat,nx,mu,nu,nb)
%
%   Returns the plaquette.
%
%          Us_3
%    x+nu -->-- x+nu+mu
%      |         |
% Us_4 ^    P    v Us_2
%      |         |
%      x  --<-- x+mu
%          Us_1
% 
%   Bjoern Leder, 2009
%
% See also: IMPORT_LATTICE

Us=plaq_us(lat,nx,mu,nu,nb);
P = Us{1}*Us{2}*Us{3}'*Us{4}';