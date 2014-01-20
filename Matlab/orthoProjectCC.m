function out = orthoProjectCC(fin, filterType, dir)

  %The different 1D derivative filters 
  dfilters = {...
      [1/2, 0, -1/2], ... 
      [1/24, 5/12, 0, -5/12, -1/24], ...
      [1/720, 7/90, 49/144, 0, -49/144, -7/90, -1/720],...
      [1/720, 7/90, 49/144, 0, -49/144, -7/90, -1/720],...
      [1/40320, 41/6720, 289/2880, 809/2880, 0, -809/2880, -289/2880, -41/6720, -1/40320]...
  };

% The corresponding sequences that the 1D filters are tensor multiplied
% with
dseqs = { [1/6, 2/3, 1/6], ...
    [1/120, 13/60, 11/20, 13/60, 1/120], ...
    [1/5040, 1/42, 397/1680, 151/315, 397/1680, 1/42, 1/5040], ...
    [1/5040, 1/42, 397/1680, 151/315, 397/1680, 1/42, 1/5040], ...
    [1/362880, 251/181440, 913/22680, 44117/181440, 15619/36288, 44117/181440, 913/22680, 251/181440, 1/362880]...    
};

 % AC sequences for linear and cubic splines, (duals)
 acseqs = {[1/6, 2/3, 1/6], ...
     [1/5040, 1/42, 397/1680, 151/315, 397/1680, 1/42, 1/5040] };
     
 % Sampled seqeuences for the cubic and quintic splines 
 sseqs = {[1/6, 2/3, 1/6], ...
 [1/120, 13/60, 11/20, 13/60, 1/120]};
 
 padsize = 10;

 dfIndex = 0;
 ssIndex = 0;
 acIndex = 0;
 
switch lower(filterType)
   case 'll' 
       dfIndex = 1;
       acIndex = 1;
   case 'cl'
       dfIndex = 2;
       acIndex = 1;
       ssIndex = 1;
    case 'cc'
       dfIndex = 3;
       acIndex = 2;
       ssIndex = 1;
    case 'ql'
       dfIndex = 4;
       acIndex = 1;
       ssIndex = 2;
    case 'qc'
       dfIndex = 5;
       acIndex = 2;
       ssIndex = 2;
    otherwise
        error('Unknown filter type.')
end

% Take care of prefiltering first then pad

if (ssIndex ~= 0)
    prefilter = makeFilterCC(tpFilter(sseqs{ssIndex}, sseqs{ssIndex}, sseqs{ssIndex}), size(fin,2), size(fin,1), size(fin,3));
    fin = real(ifftn(fftn(fin) ./ fftn(prefilter)));
end
fin = padarray(fin, [padsize padsize padsize], 'replicate');

% The derivative filter in the x direction
dfilter =tpFilter(dfilters{dfIndex}, dseqs{dfIndex}, dseqs{dfIndex});

% permute for the appropriate derivative if necessary
switch dir
    case 1
        dfilter = permute(dfilter, [2 1 3]);
    case 2
        dfilter = permute(dfilter, [1 3 2]);
    otherwise
end

dfilterGrid = makeFilterCC(dfilter, size(fin,2), size(fin,1), size(fin,3));
acfilterGrid = makeFilterCC(tpFilter(acseqs{acIndex}, acseqs{acIndex}, acseqs{acIndex}), size(fin,2), size(fin,1), size(fin,3));

% Now perform the filtering in the frequency domain

out = ifftn(fftn(fin) .* fftn(dfilterGrid) ./ fftn(acfilterGrid));
out = real(out(padsize+1:end-padsize, padsize+1:end-padsize, padsize+1:end-padsize ));

%----------------------
end