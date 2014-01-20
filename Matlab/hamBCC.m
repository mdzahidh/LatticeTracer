function [hama varargout] = hamBCC(n, gridType )
% Generates samples of the HAM dataset in the region [-1,1]^3.
% n specifies the number of samples along the x- axis
% 
% CC usage:
%   f = hamBCC(n, 'CC')
%
% BCC usage:
%  [fa fb] = hamBCC(n, 'BCC')
%    fa and fb are the unshifted and shifted Cartesian grids.
%    Each grid contains n x n x n samples.

switch lower(gridType)
    case 'cc'
        h = 2/(n-1);
    case 'bcc'
        h = 4/(2*n-1);
    otherwise
        error 'Unknown grid type. Specify either CC or BCC'.
end

% Parameters for ham dataset

fm = 6;
amp = 0.25;
beta = 2;


[xa ya za] = meshgrid(-1:h:1);
xb = xa + h/2;
yb = ya + h/2;
zb = za + h/2;

rhoA = sqrt(xa.^2 + ya.^2 + za.^2);
rhoB = sqrt(xb.^2 + yb.^2 + zb.^2);

hama = beta.*rhoA - amp.*cos(2*pi*fm*za./rhoA);
hamb = beta.*rhoB - amp.*cos(2*pi*fm*zb./rhoB);

hama(isnan(hama)) = 0.0;
hamb(isnan(hamb)) = 0.0;

% check output arguments

switch lower(gridType)
    case 'cc'
        if nargout > 1
            error 'Only one output expected for CC'.
        end
    case 'bcc'
        if nargout < 2
            error 'Two output arguments expected for BCC'
        end
        varargout{1} = hamb;
end

return
