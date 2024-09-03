% sint.m
%
%   function [int,type] = sint(v)
%
% Automatically converts v to the right integer type
% unit8, ..., uint64 or nit8, ..., int64
%

function [int,type] = sint(v)

if nargin < 1,
    disp('Too few arguments! Usage:');
    help uint;
    return;
end


ma = max(v);
while length(ma)>1, ma=max(ma); end
mi = min(v);
while length(mi)>1, mi=min(mi); end
if ma<0 | abs(ma)<abs(mi), ma = abs(mi)-1; end

if mi<0, 
    ma = 2*ma+1;
    type='int';
else
    type='uint';
end

if ma<=255,
    type = [type '8'];
elseif ma<=65535,
    type = [type '16'];
elseif ma<=4294967295,
    type = [type '32'];
elseif ma<=18446744073709551615,
    type = [type '64'];
end

int = eval([type '(v)']);


