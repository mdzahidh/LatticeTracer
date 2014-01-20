function [ofa ofb] = orthoProjectBCC(fa, fb, filterType, dir)
% Derivative computation via orthogonal projection for the BCC lattice.
%
% Usage:
% [ofa ofb] = orthoProjectBCC(fa, fb, filterType, dir)
% where:
% fa and fb are the unshifted and shifted grids that comprise the BCC
% sampling. They are assumed to contain samples.
% filterType is one of 'LL', 'QL', 'QQ', 'NL' or 'NQ'. When the first
% letter is 'Q' or 'N', the sample values are prefiltered.
% dir is 0 for x-derivatives, 1 for y-derivatives and 2 for z
% derivatives.

fname = '';
pfname = '';
acname = '';

padsize = 10;

switch lower(filterType)
   case 'll'
      fname = 'CdLLx';
      acname = 'ACLbs';
   case 'ql'
      fname = 'CdQLx';
      acname = 'ACLbs';
      pfname = 'SSQbs';
   case 'qq'
      fname = 'CdQQx';
      acname = 'ACQbs';
      pfname = 'SSQbs';
    case 'nl'
      fname = 'CdNLx';
      acname = 'ACLbs';
      pfname = 'SSNbs';
    case 'nq'
      fname = 'CdNQx';
      acname = 'ACQbs';
      pfname = 'SSNbs';
    otherwise
        error('Unknown filter type.')
end

% First prefilter if necessary
if ~strcmp(pfname, '')
    disp('Prefiltering')
    osize = size(fa);
    [prefiltera prefilterb] = readFilter(pfname,osize(2), osize(1), osize(3));
    
    [fa fb] = ifft3bcc( fft3bcc(fa, fb) ./ fft3bcc(prefiltera, prefilterb) );
    fa = real(fa);
    fb = real(fb);
end

% We'll add some padding on the boundary to reduce boundary artefacts.

pfa = padarray(fa, [padsize padsize padsize], 'replicate');
pfb = padarray(fb, [padsize padsize padsize], 'replicate');

sa = size(pfa);

% Read the appropriate derivative filter.

if ((dir ~= 0) && (dir ~= 1) && (dir ~= 2))
    error('dir must be 0, 1 or 2');
end
[dfilta dfiltb] = readFilter(fname, sa(2), sa(1), sa(3), dir);

acseqa = zeros(sa);
acseqb = zeros(sa);

disp('Reading autocorrelation sequence');
[acseqa acseqb] = readFilter(acname, sa(2), sa(1), sa(3));

% Now perform the actual filtering in the frequency domain
% where it is a simple multiplication

fftgrid = fft3bcc(pfa, pfb);
fftfilt = fft3bcc(dfilta, dfiltb);
fftacseq = fft3bcc(acseqa, acseqb);

fftresult = fftgrid .* fftfilt ./ fftacseq;

[tempofa tempofb] = ifft3bcc(fftresult);

% Final result is obtained by removing the padding

ofa = real(tempofa(padsize+1:end-padsize, padsize+1:end-padsize, padsize+1:end-padsize ));
ofb = real(tempofb(padsize+1:end-padsize, padsize+1:end-padsize, padsize+1:end-padsize ));
return