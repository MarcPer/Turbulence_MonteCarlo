function phz = generateScreen(params)

narginchk(1,2);

r0 = params.friedCoherenceRadiusMatrix;
[NxEff, NyEff] = params.getPhaseScreenGridSize;
delta = params.gridSpacingVector;
nPlanes = params.numberOfPhasePlanes;
L0 = params.outerScale;
l0 = params.innerScale;
gammaIndex = params.gammaIndex;

N = NxEff;      % Temporary, before non-square grid is implemented
NyEff = NxEff;  % Temporary, before non-square grid is implemented

phz = zeros(NyEff,NxEff,nPlanes);

for iScr = 1 : nPlanes
    if isinf(r0(gammaIndex, iScr))
        phz(:,:,iScr) = 0;
    else
        [phz_lo, phz_hi] = ft_sh_phase_screen( ...
            r0(gammaIndex, iScr), N, delta(iScr), L0, l0);
        phz(:,:,iScr) = phz_lo + phz_hi;
    end
end