function out = writevudBCC( filename, varargin )
% Barebones writer for CC or BCC volumetric data in vud format
%
% CC usage  :  out = writevudBCC( filename, v)
% BCC usage:  out = writevudBCC( filename, va, vb)

% where
% va and vb are the two 3D interleaved Cartesian grids. They must be of
% the same size.
% filename is the file where the data is to be written. The data is
% written as floats.
%
% out is the total number of elements successfully written

if (nargin == 3)
    va = varargin{1};
    vb = varargin{2};
elseif (nargin == 2)
    va = varargin{1};
else
    error 'Incorrect usage!'
end

% First check to make sure the data is the same size
% and that is contains 3 dimensions
if (nargin == 3)
    if( sum(size(va) == size(vb))  ~= 3 )
        error('Input data is not in the correct format');
    end
elseif (nargin == 2)
    if (ndims(va) ~= 3)
        error('Input data must be 3D');
    end
end

fid = fopen( filename, 'wb');
if( fid < 0 )
    error(['Could not open file ',filename, ' for writing']);
end

% Rearrange the data so that it is in one array
if (nargin == 3)
    vol = zeros([size(va,1),size(va,2),2*size(va,3)], class(va));
    vol(:,:,1:2:end) = va;
    vol(:,:,2:2:end) = vb;
elseif (nargin == 2)
    vol = va;
end

% write the header information

fprintf(fid, '# vu DataFile Version 1.0\n');
fprintf(fid, 'Volume Data\n');
fprintf(fid, 'BINARY\n');

if (nargin == 3)
    fprintf(fid, 'DATASET BCC_POINTS\n');
    disp('Writing BCC_POINTS');
elseif (nargin == 2)
    fprintf(fid, 'DATASET CC_POINTS\n');
    disp('Writing CC_POINTS');
end

fprintf(fid, 'UNIMODAL\n');
fprintf(fid, 'DIMENSIONS %d %d %d 1\n', size(vol,2), size(vol,1), size(vol,3));
fprintf(fid, 'ORIGIN 0 0 0 0\n');
fprintf(fid, 'OFFSET 0 0 0 0\n');
fprintf(fid, 'SCALING 1.6 1.6 1.6 1\n'); %dummy scaling
fprintf(fid, 'SPACING 0 0 0 0\n');
fprintf(fid, 'POINT_DATA %d\n', numel(vol));
fprintf(fid, 'SCALARS data float\n');
fprintf(fid, 'LOOKUP_TABLE default\n');

out = fwrite(fid, permute(vol, [2 1 3]), 'single');
fclose(fid);
end
