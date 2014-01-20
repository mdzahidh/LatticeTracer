function [va varargout] = readvudBCC( filename )
% barebones reader for a vud file
%
% For CC data use:
%  v = readvudBCC( filename )
% 
% For BCC data use:
%  [va vb] = readvudBCC( filename ) 
% where
% filename is the name of a binary vud file.
% va and vb contain the two interleaved Cartesian grids of the
% BCC volumetric data

fid = fopen(filename, 'rb');
if( fid < 0 )
   error(['Cannot open file ', filename] ); 
end

% Parse through the header
fgets(fid);
fgets(fid);
fgets(fid);

% type of the grid
gridType = fgets(fid);
if( ~isempty(findstr( gridType, 'BCC_POINTS')) )
    display('BCC data found');
    mode = 'BCC';
elseif( ~isempty(findstr( gridType, 'CC_POINTS')) )
    display( 'CC data found');
    mode = 'CC';
else
    error('unknown grid type');
end

%Skip another line
fgets(fid);

% The resolution of the volume
res = sscanf(fgets(fid), '%*s %d %d %d');
res = res';

% skip origin, offset, scaling and spacing
fgets(fid);
fgets(fid);
fgets(fid);
fgets(fid);

% Number of elements
num = sscanf(fgets(fid), '%*s %d');

% The kind of data
kind = sscanf(fgets(fid), '%*s %*s %s');

%skip lookup stuff
% Bug fixed, what if the first data stored is 
% encoded as 0A or 0D, in which case the reader
% will read one extra character and screw up!
% So we will read exactly 21 bytes
%fgets(fid);
fread(fid,21,'uint8');

% Now read the data

if( strcmp(kind, 'byte'))
    dataType = 'uint8';
else if( strcmp(kind, 'float'))
        dataType = 'single';
    else
        error( ['Data type ', kind, ' not supported'] );
    end
end

[data, count] = fread(fid, num, dataType);
fclose(fid);

data = reshape(data, res(1), res(2), res(3));
data = permute(data, [2 1 3]);

switch lower(mode)
    case 'cc'
        va = data;
    case 'bcc'
        if nargout < 2
            error('Output arguments not correct for BCC data. Two are expected');
        end
        % va and vb are now the unshifted and shifted versions

        % If the number of slices is not even, add a slice of
        % zeros

        if( ~(mod(res(3),2) == 0) )
            data = padarray(data, [0 0 1], 0, 'post');
        end

        va = data(:,:,1:2:end);
        varargout{1} = data(:,:,2:2:end);
end

return