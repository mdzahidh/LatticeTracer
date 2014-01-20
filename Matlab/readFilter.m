function [fa fb] = readFilter(name, nX, nY, nZ, varargin)
% Read a BCC filter and create a Cartesian periodic extension in the first
% octant.
%
% Usage:
% [fa fb] = readFilter(name, nX, nY, nZ)
% [fa fb] = readFilter(name, nX, nY, nZ, dir)
%
% where:
% name is one of ''.
% numX, numY and numZ specify the size of the grid to create.
% dir is an optional argument that permutes the coordinates of the filter. 
%  1 - permute x and y, 2 - permute x and z. This is useful for derivative
%  filters.
%
% fa and fb are the two nX x nY x nZ grids that contain the read filter.

temp = struct2cell(load(name));
filt = temp{1,1};
clear temp;

% see if we need to permute any of the coordinates.

if nargin > 4
    if varargin{1} == 1
        disp('Permuting x and y');
        filt = [filt(:,2) filt(:,1) filt(:,3) filt(:,4)];
    elseif varargin{1} == 2
        disp('Permuting x and z');
        filt = [filt(:,3) filt(:,2) filt(:,1) filt(:,4)];
    end
end

% filt consits of 4 columns, the first three contain the coordinates
% and the fourth contains the filter weights

fa = zeros(nY, nX, nZ);
fb = zeros(nY, nX, nZ);

[uX uY uZ] = meshgrid(-floor(nX/2):1:-floor(nX/2)+nX-1, -floor(nY/2):1:-floor(nY/2)+nY-1, -floor(nZ/2):1:-floor(nZ/2)+nZ-1);
uX = 2*uX; uY = 2*uY; uZ = 2*uZ; 
sX = uX + 1;
sY = uY + 1;
sZ = uZ + 1;

% The unshifted and shifted grid points

uPoints = [reshape(uX,numel(uX),1) reshape(uY,numel(uY),1) reshape(uZ,numel(uZ),1)];
sPoints = [reshape(sX,numel(sX),1) reshape(sY,numel(sY),1) reshape(sZ,numel(sZ),1)];

% fill in the unshifted and shifted points

[temp, uIndicesFilt, uIndicesPoints] = intersect(filt(:,1:3), uPoints, 'rows');
clear temp;
fa(uIndicesPoints) = filt(uIndicesFilt, 4);

[temp, sIndicesFilt, sIndicesPoints] = intersect(filt(:,1:3), sPoints, 'rows');
clear temp;
fb(sIndicesPoints) = filt(sIndicesFilt, 4);

fa = ifftshift(fa);
fb = ifftshift(fb);

return