function [lat] = sw_term(lat)
% sw_term.m Computes the SW-term for lat.U
%
% It is stored in lat.sw
%
% Bjoern Leder, 2011
%

F=field_strength(lat);
G=gamma_matrices_4d;
S=cell(6);
for n=1:6,
    S{n}=1i/2*G{5+n};
end

X=sparse(12*lat.volume,12);
lat.sw=mat2cell(X,12*ones(lat.volume,1),12);
% lat.sw=cell(lat.volume,1);
% [lat.sw{:}]=deal(sparse(12,12));
for nx=1:lat.volume,
    for n=1:6,
        lat.sw{nx}=lat.sw{nx}+kron(S{n},F{n,nx});
    end
end