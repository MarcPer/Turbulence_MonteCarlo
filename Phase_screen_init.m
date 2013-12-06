% SUB-SCRIPT of Phase_screen_main.m
%   Initialize some variables

phz = zeros(N, N, n);       % Phase screen arrays
Iout_real = zeros(N);       % Output intensity for a given realization
Iout = zeros([N,N,leng]);   % Output intensities for all gammas
Uout_real = zeros(N);       % Output field in a single realization
sg = repmat(sg, [1 1 n]);   % Super-gaussian filter
Pslit = zeros(leng,1);      % Power integrated over slit (single counts)
