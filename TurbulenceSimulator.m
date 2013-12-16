classdef TurbulenceSimulator<handle
   properties (Access = private)
       simulationParameters;
       numberOfTransverseSeparations;
       phaseScreenProfiles;
   end
   
   methods
       function ts = TurbulenceSimulator(simParams)
           ts.simulationParameters = simParams;
           ts.numberOfTransverseSeparations = length(...
               simParams.transverseSeparationInR0Units);
       end
   end
   methods(Access = public)
       function outputField = propagate(obj, phScreen)
           Uin = obj.simulationParameters.getInputField;
           sg = obj.simulationParameters.getSuperGaussianFilter;
           wvl = obj.simulationParameters.wavelength;
           delta1 = obj.simulationParameters.gridSpacingSourcePlane;
           deltan = obj.simulationParameters.gridSpacingObservationPlane;
           
           [~,~, outputField] = ang_spec_multi_prop(Uin, wvl, ...
            delta1, deltan, z, sg.*exp(1i*phScreen));
       end
       function intSep = getFieldForEachTransverseSeparation(obj)
           nSep = obj.numberOfTransverseSeparations;
           [Nx, Ny] = obj.simulationParameters.getEffectiveGridSize;
           intSep = zeros(Ny, Nx, nSep);
           
       end
       function IavgRe = getAverageOverRealizations(obj)
           nRe = obj.simulationParameters.numberOfRealizations;
           nSep = obj.numberOfTransverseSeparations;
           [Nx, Ny] = obj.simulationParameters.getEffectiveGridSize;
           IavgRe = zeros(Ny, Nx, nSep);
           
           for iRe = 1 : nRe
               obj.phaseScreenProfiles = generateScreen(obj.simulationParameters);
               outputField =  ...
                   obj.getFieldForEachTransverseSeparation();
               intensitySingleRealization = outputField.^2;
               IavgRe = IavgRe + intensitySingleRealization;
           end
           IavgRe = IavgRe/nRe;
           
       end
       function intGamma = getIntensityForEachGamma(obj)
           nGamma = length(obj.simulationParameters.gammaStrength);
           intGamma = cell(nGamma);
           
           for iGamma = 1 : nGamma
               obj.simulationParameters.gammaIndex = iGamma;
               intGamma{iGamma} = obj.getAverageOverRealizations();
           end
       end
       % Returns cell{idxGamma} = int(Ny,Nx,separationIndex)
   end
end