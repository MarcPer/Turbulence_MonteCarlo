function phz = generateScreen(params)

r0 = params.friedCoherenceRadiusMatrix;
[Nx, Ny] = params.getPhaseScreenGridSize;
delta = params.gridSpacingVector;
nPlanes = params.numberOfPhasePlanes;
L0 = params.outerScale;
l0 = params.innerScale;
gammaIndex = params.gammaCurrentIndex;

N = Nx;      % Temporary, before non-square grid is implemented
Ny = Nx;     % Temporary, before non-square grid is implemented

phz = zeros(Ny,Nx,nPlanes);

for iScr = 1 : nPlanes
    if isinf(r0(gammaIndex, iScr))
        phz(:,:,iScr) = 0;
    else
        [phz_lo, phz_hi] = ft_sh_phase_screen( ...
            r0(gammaIndex, iScr), N, delta(iScr), L0, l0);
        phz(:,:,iScr) = phz_lo + phz_hi;
    end
end