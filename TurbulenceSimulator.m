classdef TurbulenceSimulator<handle
    properties (Access = private)
        simulationParameters;
        numberOfTransverseSeparations;
        phaseScreenProfiles;
        abortButtonHandle;
    end
    
    properties (SetAccess = private, GetAccess = public)
        isAborted; 
    end
    
    methods
        function ts = TurbulenceSimulator(simParams)
            ts.isAborted = 0;
            ts.simulationParameters = simParams;
            ts.numberOfTransverseSeparations = length(...
                simParams.transverseSeparationInR0Units);
        end
    end
    methods(Access = public)
        function outputField = propagate(obj, phScreen)
            Uin = obj.simulationParameters.getInputField;
            sg = obj.simulationParameters.getSuperGaussianFilter;
            sg = repmat(sg, [1, 1, obj.simulationParameters.numberOfPhasePlanes]);
            wvl = obj.simulationParameters.wavelength;
            z = obj.simulationParameters.planePositions;
            delta1 = obj.simulationParameters.gridSpacingSourcePlane;
            deltan = obj.simulationParameters.gridSpacingObservationPlane;
            
            [~,~, outputField] = ang_spec_multi_prop(Uin, wvl, ...
                delta1, deltan, z, sg.*exp(1i*phScreen));
        end
        function fieldSep = getFieldForEachTransverseSeparation(obj)
            nSep = obj.numberOfTransverseSeparations;
            phz = obj.phaseScreenProfiles;
            [Nx, Ny] = obj.simulationParameters.getTransverseGridSize;
            
            fieldSep = zeros(Ny, Nx, nSep);
            
            for iSep = 1 : nSep
                obj.isAborted = UserInput.isAborted(obj.abortButtonHandle);
                if obj.isAborted
                    break;
                end
                [deltaX, ~] = obj.simulationParameters.getTransverSeparationInPixels(iSep);
                phScreen = Util.displace(phz,deltaX,0);
                phScreen = Util.crop(phScreen, Nx, Ny);
                fieldSep(:,:,iSep) = obj.propagate(phScreen);
            end
        end
        function IavgRe = getAverageOverRealizations(obj)
            obj.abortButtonHandle = UserInput.createWaitBar;
            try
                nRe = obj.simulationParameters.numberOfRealizations;
                nSep = obj.numberOfTransverseSeparations;
                [Nx, Ny] = obj.simulationParameters.getTransverseGridSize;
                IavgRe = zeros(Ny, Nx, nSep);
                
                for iRe = 1 : nRe
                    obj.isAborted = UserInput.isAborted(obj.abortButtonHandle);
                    if obj.isAborted
                        break;
                    end
                    obj.phaseScreenProfiles = generateScreen(obj.simulationParameters);
                    outputField =  ...
                        obj.getFieldForEachTransverseSeparation();
                    intensitySingleRealization = abs(outputField).^2;
                    IavgRe = IavgRe + intensitySingleRealization;
                    UserInput.updateWaitBar(obj.abortButtonHandle, iRe, nRe);
                end
                delete(obj.abortButtonHandle);
                IavgRe = IavgRe/nRe;
            catch exception
                delete(obj.abortButtonHandle);
                rethrow(exception);
            end
        end
        function [intGamma, plotInfo] = getIntensityForEachGamma(obj)
            nGamma = length(obj.simulationParameters.gammaStrength);
            intGamma = cell(nGamma);
            
            for iGamma = 1 : nGamma
                if obj.isAborted
                    break;
                end
                UserInput.printOutProgress('Turbulence strength', ...
                    iGamma, nGamma);
                obj.simulationParameters.gammaIndex = iGamma;
                intGamma{iGamma} = obj.getAverageOverRealizations();
            end
            
            tit = 'Intensity vs Turbulence Strength';
            labelX = '\gamma';
            labelY = 'intensity';
            plotInfo = struct('title', tit, ...
                'labelX', labelX, 'labelY', labelY);
        end
        % Returns cell{idxGamma} = int(Ny,Nx,separationIndex)
    end
    
end