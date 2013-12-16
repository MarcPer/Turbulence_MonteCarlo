function phz = generateScreen(params)

narginchk(1,2);

%% CHECK parameter names!!
r0 = params.friedCoherenceRadiusMatrix;
N = params.transverseGridSize;
delta = params.gridSpacingVector;
nPlanes = params.numberOfPhasePlanes;
L0 = params.outerScale;
l0 = params.innerScale;
gammaIndex = params.gammaIndex;
%%
phz = zeros(N,N,nPlanes);

for iScr = 1 : nPlanes
    if isinf(r0(gammaIndex, iScr))
        phz(:,:,iScr) = 0;
    else
        [phz_lo, phz_hi] = ft_sh_phase_screen( ...
            r0(gammaIndex, iScr), N, delta(iScr), L0, l0);
        phz(:,:,iScr) = phz_lo + phz_hi;
    end
end