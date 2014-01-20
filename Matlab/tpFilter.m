function tp = tpFilter(f1, f2, f3)
% Creates a tensor product extension 3D filter from a 1D filters (row vectors)
n = numel(f1);
tp = repmat(f2'*f1, [1 1 n]) .* padarray(reshape(f3, [1 1 n]), [floor(n/2) floor(n/2)], 'replicate');
end