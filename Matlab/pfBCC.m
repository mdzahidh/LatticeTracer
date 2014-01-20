function out = pfBCC( filein, fileout, gridType )
% Prefilters data for either tricubic B-spline interpolation on CC or
% Quintic box spline interpolation on BCC.
%
% Usage:
% out = pfBCC(filein, fileout, gridType)
%
% where: filein specifies the input file
%             fileout specified the name of the output file 
%            gridType is one of 'CC' or 'BCC'
%
% out returns the number of elements written
% Data is written in single precision floating point vud format.


% Read data, make the filter, perform filtering, and write the data

switch lower(gridType)
    case 'cc'
        va = readvudBCC(filein);
        f = [1/6 2/3 1/6];
        cf = makeFilterCC(tpFilter(f,f,f), size(va,2), size(va,1), size(va,3));
        va = real(ifftn(fftn(va) ./ fftn(cf)));
        out = writevudBCC(fileout, va);
    case 'bcc'
        [va vb] = readvudBCC(filein);
        [filta, filtb] = readFilter('SSQbs', size(va,2), size(va,1), size(va,3));
        [va vb] = ifft3bcc(fft3bcc(va, vb) ./ fft3bcc(filta, filtb) );
        va = real(va);
        vb = real(vb);
        out = writevudBCC(fileout, va,vb);
end

return
