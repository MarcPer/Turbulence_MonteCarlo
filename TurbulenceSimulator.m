classdef TurbulenceSimulator<handle
    % Takes phase screens and input field and propagates the latter.
    % Includes a method to average over many realizations and perform
    % various operations on the screens (e.g. inversion and displacement)
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
            phz = obj.inversionAndFourthOrderOperationsOnScreen();
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
                nSep = obj.numberOfTransverseSeparations;
                nRe = obj.getNumberOfRealizations;
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
                obj.abortButtonHandle = [];
                IavgRe = IavgRe/nRe;
            catch exception
                if ~isempty(obj.abortButtonHandle)
                    delete(obj.abortButtonHandle);
                    obj.abortButtonHandle = [];
                end
                rethrow(exception);
            end
        end
        function intGamma = getIntensityForEachGamma(obj,varargin)
            inParams = Util.transformInputParametersIntoStructure(varargin);
            isNormalized = false;
            if isfield(inParams,'Normalized')
                isNormalized = inParams.Normalized; 
            end
            
            nGamma = length(obj.simulationParameters.gammaStrength);
            intGamma = cell(nGamma, 1);
            
            for iGamma = 1 : nGamma
                if obj.isAborted
                    break;
                end
                UserInput.printOutProgress('Turbulence strength', ...
                    iGamma, nGamma);
                obj.simulationParameters.gammaCurrentIndex = iGamma;
                intGamma{iGamma} = obj.getAverageOverRealizations();
                if isNormalized
                    intGamma{iGamma} = Util.normalize(intGamma{iGamma});
                end
            end
        end
        % Returns cell{idxGamma} = int(Ny,Nx,separationIndex)
    end
    methods(Access = private)
        function phScreen = inversionAndFourthOrderOperationsOnScreen(obj) 
            phScreen = obj.phaseScreenProfiles;
            if ~(obj.simulationParameters.isFourthOrder)
                return;
            end
            
            if obj.simulationParameters.isInverted
                phScreen = phScreen + Util.rot90All(phScreen,2);
            else
                phScreen = 2*phScreen;
            end
        end
        function nRe = getNumberOfRealizations(obj)
            nReInput = obj.simulationParameters.numberOfRealizations;
            idxGamma = obj.simulationParameters.gammaCurrentIndex;
            r0 = obj.simulationParameters.totalFriedCoherenceRadiusByStrength;
            if isinf(r0(idxGamma))
                nRe = 1;
            else
                nRe = nReInput;
            end
        end
    end
end