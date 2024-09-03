function [lat] = import_ddhmc_lattice(fname,bc)
%
%   Imports lattice exported by DD-HMC-1.2.1 code.
%
%   [LAT] = IMPORT_DDHMC_LATTICE(FNAME)
%
%   Reads in binary mode from file FNAME.
%
%   See misc/archive.c in the DD-HMC-1.2.1 for the format specification.
%
%   LAT is a stucture containing the information in the header of the file
%   and the gauge field LAT.U has cell array. Each element is one link variable.
% 
%   LAT.U{4*(ix-1)+1+mu} is the link variable at the point with label ix in
%   direction mu=0,...,3. The label of a point with cartesian coordinates
%   x0, x1, x2, x3 is 
%
%       ix = x3 + x2*lat.L(3) + x1*prod(lat.L(2:3)) + x0*prod(lat.L)
%
%   Bjoern Leder, 2008
%
% See also: IMPORT_LATTICE, STRUCT, CELL

if nargin<2
    error('import_ddhmc_lattice:arguments','Not enough input arguments.')
end

fid = fopen(fname,'r');
if fid<0
    error('import_ddhmc_lattice:fileNotFound',['Unable to open file ' fname '.'])
end

lat=empty_lattice;

lat.bc=bc;
lat.colour=3;
lat.flavour=0; % default
lat.zf=1; % default
lat.ds=1; % default
lat.T=fread(fid,1,'int32');
lat.L=fread(fid,3,'int32');
lat.plaq=fread(fid,1,'double');
lat.dirac=length(lat.L)+1;
lat.volume=prod(lat.L)*lat.T;

tic;
su3_matrix=3*3*2;
u=fread(fid,lat.volume*lat.dirac*su3_matrix,'double');
if length(u)~=lat.volume*lat.dirac*su3_matrix, die(fid); end
eof=isempty(fread(fid,1,'char'));
fclose(fid);
if ~eof
    error('import_ddhmc_lattice:ReadError','Inconsistent lattice dimensions and file size.')
end

% nearest neighbour (periodic BC!)
NB=double(nextnb(lat.dirac-1,lat.T,lat.L,0,2));

os=0;
for nx=1:1:lat.volume
    [x3 x2 x1 x0]=ind2sub([lat.L(3) lat.L(2) lat.L(1) lat.T],nx);
    if mod(x1+x2+x3+x0,2)~=0
        for mu=0:lat.dirac-1
            % + mu-direction
            unx=lat.dirac*(nx-1)+1;
            lat.U{unx+mu}=reshape(u(os+1:2:os+su3_matrix)+1i*u(os+2:2:os+su3_matrix),3,3).';
            os=os+su3_matrix;
            % - mu-direction
            unx=lat.dirac*(NB(nx,mu+1+lat.dirac)-1)+1;
            if x0==0
                lat.U{unx+mu}=reshape(u(os+1:2:os+su3_matrix)-1i*u(os+2:2:os+su3_matrix),3,3);
            else
                lat.U{unx+mu}=reshape(u(os+1:2:os+su3_matrix)+1i*u(os+2:2:os+su3_matrix),3,3).';
            end
            os=os+su3_matrix;
        end
    end
end
toc

function die(fid)
    fclose(fid);
    error('import_ddhmc_lattice:read','Read error.');
