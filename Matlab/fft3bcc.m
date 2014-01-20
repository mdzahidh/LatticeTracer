function out = fft3bcc(fa, fb)
  % A function to compute a Discrete Fourier Transform 
  % on a bcc grid using the FFT.
  % fa and fb both consist of M cols and N rows and S slices and represent
  % the two interleaved Cartesian grids of a bcc lattice
  
  M = size(fa,2);
  N = size(fa,1);
  S = size(fa,3);
  
  % determine the offsets
  
  u_t = -1 * ceil(-M/2);
  v_t = -1 * ceil(-N/2);
  w_t = -1 * ceil(-S/2);
  
  [m n s] = meshgrid(0:1:M-1, 0:1:N-1, 0:1:S-1);
  
  % The transform can be split up into two parts
  
  out1 = zeros(N, M, 2*S); % for the unshifted grid
  out2 = zeros(N, M, 2*S); % for the shifted grid
  
  temp = exp(2*pi*i*(m.*(u_t/M) + n.*(v_t/N) +s.*(w_t/S)) );
  
  %
  % first do the unshifted case.
  %
  
  % Slices 0:S-1 are given by a regular FFT
  out1(:,:,1:S) = fftn(fa .* temp);
  
  % Slices S:2S-1 are given by the same FFT
  out1(:,:,S+1:end) = out1(:,:,1:S); 
 
  %
  % Now do the shifted case
  %
  [u v w] = meshgrid(0:1:M-1, 0:1:N-1, 0:1:S-1);

  tempfft = fftn(fb.*temp);
  
  % Slices 0:S-1 are given by a regular FFT
 
  out2(:,:,1:S) = exp(pi * i * (u_t/M + v_t/N + w_t/S)) .* ...
      exp(-pi * i * (u./M + v./N + w./S) ) .* tempfft;
  
  % Slices S:2S-1 are modulated
  out2(:,:,S+1:end) = exp(pi * i * (u_t/M + v_t/N + w_t/S)) .* ...
      exp(-pi * i * (u./M + v./N + (w+S)./S) ) .*...
      tempfft;
      

  % Now the result is simply the two FFTs added
  
  out = out1 + out2;
return