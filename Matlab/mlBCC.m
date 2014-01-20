function [mla varargout] = mlBCC(n, gridType )
% Generates samples of the ML dataset in the region [-1,1]^3.
% n specifies the number of samples along the x- axis
% 
% CC usage:
%   f = mlBCC(n, 'CC')
%
% BCC usage:
%  [fa fb] = mlBCC(n, 'BCC')
%    fa and fb are the unshifted and shifted Cartesian grids.
%    Each grid contains n x n x n samples.

% Generate a grid for ML sampled on a BCC grid

switch lower(gridType)
    case 'cc'
        h = 2/(n-1);
    case 'bcc'
        h = 4/(2*n-1);
    otherwise
        error 'Unknown grid type. Specify either CC or BCC'.
end

% Parameters for ML dataset

fm = 6;
alpha = 0.25;

[xa ya za] = meshgrid(-1:h:1);
xb = xa + h/2;
yb = ya + h/2;
zb = za + h/2;

rhoA = cos(2*pi*fm*cos(0.5*pi*sqrt( (xa).^2 + (ya).^2 ) ));
rhoB = cos(2*pi*fm*cos(0.5*pi*sqrt( (xb ).^2 + (yb).^2 ) ));

mla = (1-sin(0.5*pi*(za)) + alpha*(1+rhoA)) ./ (2+2*alpha);
mlb = (1-sin(0.5*pi*(zb)) + alpha*(1+rhoB)) ./ (2+2*alpha);

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
        varargout{1} = mlb;
end

return