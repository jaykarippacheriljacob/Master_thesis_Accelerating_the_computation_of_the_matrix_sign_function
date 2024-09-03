function [D,lat_out] = wilson_dirac_from_cnfg(filename,m0,csw,mu,lat_in,seq_idx,delta_seq)
%
% Read cnfg from file and construct Wilson-Dirac operator
% the spectrum is shifted by m0
% optionally with SW-term and finite mu
%
%  filename    full path to cnfg in DDHMC format
%  m0          shift in the dirac operator (usually negative)
%  csw         parameter multiplying the SW-term
%  mu          finite chemical potential
%              (exp(mu) is multiplied to the temporal links)
%
% Bjoern Leder, 2014
%

addpath('./matlab')

if nargin<2 || isempty(m0),
   m0=0;
end
if nargin<3 || isempty(csw),
   csw=0;
end
if nargin<4 || isempty(mu),
   mu=0;
end

% import cnfg, ie. links
if seq_idx==1
  lat_out = import_ddhmc_lattice(filename,0);
else
  lat_out = lat_in;
end

% at this point we have the gauge links in lat.U
% (see description at the beginning of import_ddhmc_lattice.m),
% and for seq_idx>1 we want to do a variation of the previous lat

if seq_idx>1
  for ix=1:length(lat_out.U)/4
    for mux=0:3
      local_U = lat_out.U{4*(ix-1)+1+mux};
      
      % build perturbation
      local_Pert = eye(3) + 1i*delta_seq*rand(3,3);
      [K,~] = qr(local_Pert);
      local_U_Pert = K*(local_U*K');
      
      % replace the output link with the perturbed version
      lat_out.U{4*(ix-1)+1+mux} = local_U_Pert;
    end
  end
end

pl=plaquette(lat_out);
disp(['average plaquette: ' num2str(pl)]);

% compute SW-term
if csw~=0
   lat_out.csw=csw;
   lat_out = sw_term(lat_out);
end

% set finite mu
if mu~=0
   lat_out.mu=mu;
end

% construct D
kappa=1/(2*m0+8); % m0=1/(2*kappa)-4;
D = wilson(kappa,0,lat_out);
