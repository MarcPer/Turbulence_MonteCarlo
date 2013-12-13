function phz = generateScreen(params, idxStrength)

narginchk(1,2);

if (nargin == 1)
    idxStrength = 1;
end

%% CHECK parameter names!!
r0 = params.friedCoherenceRadiusMatrix;
N = params.transverseGridSize;
delta = params.gridSpacingVector;
nPlanes = params.numberOfPhasePlanes;
%%
phz = zeros(N,N,nPlanes);

for iScr = 1 : nPlanes
    if isinf(r0(idxStrength, iScr))
        phz(:,:,iScr) = 0;
    else
        [phz_lo, phz_hi] = ft_sh_phase_screen( ...
            r0(idxStrength, iScr), N, delta(iScr), L0, l0);
        phz(:,:,iScr) = phz_lo + phz_hi;
    end
end