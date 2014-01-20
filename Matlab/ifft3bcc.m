function [fa fb] = ifft3bcc(tbcc)
  %A function to compute an Inverse Discrete Fourier Transform 
  % for a bcc grid using the FFT.
  % tbcc is a transform that consists of M cols, N rows and 2*S slices and 
  % 
  % the two interleaved Cartesian grids of a bcc lattice fa and fb are
  % returned
  
  M = size(tbcc,2);
  N = size(tbcc,1);
  S = 0.5*size(tbcc,3);
  
  % determine the offsets
  
  u_t = -1 * ceil(-M/2);
  v_t = -1 * ceil(-N/2);
  w_t = -1 * ceil(-S/2);
  
  [m n s] = meshgrid(0:1:M-1, 0:1:N-1, 0:1:S-1);

  temp = exp(-2*pi*i*( (u_t/M).*m + (v_t/N).*n + (w_t/S).*s) );
  
  % First compute the unshifted grid.
  
  fa = 0.5 *  temp .* splitfft(tbcc);
  
  % Now compute the shifted grid
  
  [u v w] = meshgrid(0:1:M-1,0:1:N-1,0:1:2*S-1);
  
  fb = 0.5 * exp(-pi*i*(u_t/M + v_t/N + w_t/S)) .* temp .* ...
      splitfft(tbcc .* exp(pi*i*(u./M + v./N + w./S)) );
  

  % A nested function to do a computation common to both unshifted 
  % and shifted cases
    function inv = splitfft(tr)
        inv1 = ifftn(tr(:,:,1:S));
        inv2 = exp(2*pi*i*s) .* ifftn(tr(:,:,S+1:end));
        inv = (inv1 + inv2);
    end

end