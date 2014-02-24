classdef TurbulenceSimulator<handle
    % Takes phase screens and input field and propagates the latter.
    % Includes a method to average over many realizations and perform
    % various operations on the screens (e.g. inversion and displacement)
    properties (Access = private)
        numberOfTransverseSeparations;
        phaseScreenProfiles;
        abortButtonHandle;
    end
    
    properties (SetAccess = private, GetAccess = public)
        simulationParameters;
        isAborted;
        isNormalized;
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
            numOrders = numel(obj.simulationParameters.hermiteGaussOrders);
            [Nx, Ny] = obj.simulationParameters.getTransverseGridSize;
            outputField = zeros(Ny, Nx, numOrders);
            sg = obj.simulationParameters.getSuperGaussianFilter;
            sg = repmat(sg, [1, 1, obj.simulationParameters.numberOfPhasePlanes]);
            wvl = obj.simulationParameters.wavelength;
            if obj.simulationParameters.isFourthOrder
                wvl = wvl/2;
            end
            z = obj.simulationParameters.planePositions;
            delta1 = obj.simulationParameters.gridSpacingSourcePlane;
            deltan = obj.simulationParameters.gridSpacingObservationPlane;
            
            for i = 1 : numOrders
                Uin = obj.simulationParameters.getInputField(i);
                [~,~, outputField(:,:,i)] = ang_spec_multi_prop(Uin, wvl, ...
                    delta1, deltan, z, sg.*exp(1i*phScreen));
            end
        end

        function fieldSep = getFieldForEachTransverseSeparation(obj)
            nSep = obj.numberOfTransverseSeparations;
            phz = obj.inversionAndFourthOrderOperationsOnScreen();
            [Nx, Ny] = obj.simulationParameters.getTransverseGridSize;
            
            fieldSep = zeros(Ny, Nx, nSep);
            
            for iSep = 1 : nSep
                obj.isAborted = UserInput.isAborted(obj.abortButtonHandle);
                if obj.isAborted
                    break;
                end
                [deltaX, ~] = obj.simulationParameters.getTransverseSeparationInPixels(iSep);
                phScreen = Util.displace(phz,deltaX,0);
                phScreen = Util.crop(phScreen, Nx, Ny);
                fieldSep(:,:,iSep) = obj.propagate(phScreen);
            end
        end

        function outMode = getFieldForEachMode(obj)
            phz = obj.inversionAndFourthOrderOperationsOnScreen();
            [Nx, Ny] = obj.simulationParameters.getTransverseGridSize;
            
            obj.isAborted = UserInput.isAborted(obj.abortButtonHandle);
            phScreen = Util.crop(phz, Nx, Ny);
            outMode = obj.propagate(phScreen);
        end

        function intGamma = getIrradianceForEachGamma(obj,varargin)
            % Returns cell{idxGamma} = int(Ny,Nx,separationIndex)
            inParams = Util.transformInputParametersIntoStructure(varargin);
            obj.isNormalized = false;
            if isfield(inParams,'Normalized')
                obj.isNormalized = inParams.Normalized;
            end
            
            nGamma = length(obj.simulationParameters.gammaStrength);
            intGamma = obj.fillIrradianceMetaData();
            intGamma.values = cell(nGamma, 1);
            
            for iGamma = 1 : nGamma
                if obj.isAborted
                    break;
                end
                UserInput.printOutProgress('Turbulence strength', ...
                    iGamma, nGamma);
                obj.simulationParameters.gammaCurrentIndex = iGamma;
                intGamma.values{iGamma} = obj.getAverageIrradianceOverRealizations();
                if obj.isNormalized
                    intGamma.values{iGamma} = Util.normalize(intGamma.values{iGamma});
                end
            end
        end

        function pwrAndSI = getPowerAndSIOnCircularAperture(obj,apertureRadius,varargin)
            % Returns pwr(separationIndex,gammaIndex)
            inParams = Util.transformInputParametersIntoStructure(varargin);
            obj.isNormalized = false;
            if isfield(inParams,'Normalized')
                obj.isNormalized = inParams.Normalized;
            end
            
            obj.simulationParameters.circularApertureRadius = apertureRadius;
            
            nGamma = length(obj.simulationParameters.gammaStrength);
            nSep = obj.numberOfTransverseSeparations;

            pwrGamma = obj.fillCircularApertureMetaData();
            pwrGamma.values = zeros(nSep,nGamma);
            pwrGamma.params = obj.getSimulationParameters();

            scintIdx = obj.fillScintillationIndexOnCircularApertureMetaData();
            scintIdx.values = zeros(nSep,nGamma);
            scintIdx.params = obj.getSimulationParameters();
            
            for iGamma = 1 : nGamma
                if obj.isAborted
                    break;
                end
                UserInput.printOutProgress('Turbulence strength', ...
                    iGamma, nGamma);
                obj.simulationParameters.gammaCurrentIndex = iGamma;
                [pwrGamma.values(:,iGamma), scintIdx.values(:,iGamma)] = ...
                    obj.getPowerAnsSIOnCircularApertureAveragedOverRealizations();
            end

            pwrAndSI = {pwrGamma, scintIdx};
        end

        function modeOverlap = getModeMatching(obj)
            obj.numberOfTransverseSeparations = 1;
            [Nx, Ny] = obj.simulationParameters.getTransverseGridSize();
            numOrders = numel(obj.simulationParameters.hermiteGaussOrders);
            if numOrders == 1
                fprintf('WARNING: Mode-match simulation is about to be performed with a single mode.');
            end
            nGamma = length(obj.simulationParameters.gammaStrength);
            obj.setFreeSpaceConditions();
            if  length(obj.simulationParameters.gammaStrength)<2
                error('turbuSimulator:modeMatching', 'This simulation requires a non-zero turbulence strength as input parameter.');
            end

            modeOverlap = cell(nGamma-1,1);
            refModes = zeros(Ny, Nx, numOrders);

            for i = 1 : numOrders
                refModes(:,:,i) = obj.simulationParameters.getOutputField(i);
            end

            for iGamma = 2 : nGamma
                if obj.isAborted
                    break;
                end
                 UserInput.printOutProgress('Turbulence strength', ...
                    iGamma-1, nGamma-1);
                obj.simulationParameters.gammaCurrentIndex = iGamma;
                modeOverlap{iGamma-1} = obj.fillModeOverlapMetaData;
                modeOverlap{iGamma-1}.params = obj.getSimulationParameters;

                modeOverlap{iGamma-1}.values = obj.getModeMatchingAveragedOverRealizations(refModes);
            end
        end
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

        function [pwr, scintIdx] = getPowerAnsSIOnCircularApertureAveragedOverRealizations(obj)
            obj.abortButtonHandle = UserInput.createWaitBar;
            apertureRadius = obj.simulationParameters.circularApertureRadius;
            try
                nSep = obj.numberOfTransverseSeparations;
                nRe = obj.getNumberOfRealizations;
                pwr = zeros(nSep,1);
                pwr2 = zeros(nSep,1);
                
                for iRe = 1 : nRe
                    obj.isAborted = UserInput.isAborted(obj.abortButtonHandle);
                    if obj.isAborted
                        break;
                    end
                    obj.phaseScreenProfiles = generateScreen(obj.simulationParameters);
                    outputField =  ...
                        obj.getFieldForEachTransverseSeparation();
                    intensitySingleRealization = abs(outputField).^2;
                    if obj.isNormalized
                        intensitySingleRealization = Util.normalize(intensitySingleRealization);
                    end
                    pwrSingleRealization = obj.getPowerOverCircularAperture(intensitySingleRealization, ...
                        apertureRadius);
                    pwr = pwr + pwrSingleRealization(:);
                    pwr2 = pwr2 + pwrSingleRealization(:).^2;
                    UserInput.updateWaitBar(obj.abortButtonHandle, iRe, nRe);
                end
                delete(obj.abortButtonHandle);
                obj.abortButtonHandle = [];
                pwr = pwr/nRe;
                scintIdx = pwr2/nRe * pwr.^-2 - 1;
            catch exception
                if ~isempty(obj.abortButtonHandle)
                    delete(obj.abortButtonHandle);
                    obj.abortButtonHandle = [];
                end
                rethrow(exception);
            end
        end

        function mm = getModeMatchingAveragedOverRealizations(obj,refModes)
            obj.abortButtonHandle = UserInput.createWaitBar;

            numOrders = numel(obj.simulationParameters.hermiteGaussOrders);
            nRe = obj.getNumberOfRealizations;
            mm = zeros(numOrders);

            try
                for iRe = 1 : nRe
                    obj.isAborted = UserInput.isAborted(obj.abortButtonHandle);
                    if obj.isAborted
                        break;
                    end
                    obj.phaseScreenProfiles = generateScreen(obj.simulationParameters);
                    outputModes = obj.getFieldForEachMode();
                    outputModes = Util.normalize(outputModes);

                    mm = mm + abs(Util.modeInnerProduct(refModes, outputModes)).^2;
                    UserInput.updateWaitBar(obj.abortButtonHandle, iRe, nRe);
                end
                delete(obj.abortButtonHandle);
                obj.abortButtonHandle = [];
                mm = mm / nRe;
            catch exception
                if ~isempty(obj.abortButtonHandle)
                    delete(obj.abortButtonHandle);
                    obj.abortButtonHandle = [];
                end
                rethrow(exception);
            end
        end

        function irrStruct = fillIrradianceMetaData(obj)
            irrStruct = struct;
            simParams = obj.simulationParameters;
            irrStruct.columnParams = simParams.gammaStrength;
            irrStruct.rowParams = simParams.transverseSeparationInR0Units;
            
            tit = 'Irradiance Profile';
            labelColumn = '\gamma';
            labelRow = 'Separation (in units of r0)';
            labelZ = 'Irradiance';
            irrStruct.info = struct('title', tit, ...
                'labelColumn', labelColumn, 'labelRow', labelRow, ...
                'labelZ', labelZ);
        end
        function IavgRe = getAverageIrradianceOverRealizations(obj)
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
                    IavgRe = IavgRe + abs(outputField).^2;
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
        function pwrStruct = fillCircularApertureMetaData(obj)
            pwrStruct = struct;
            simParams = obj.simulationParameters;
            pwrStruct.columnParams = simParams.gammaStrength;
            pwrStruct.rowParams = simParams.transverseSeparationInR0Units;
            
            tit = 'Power Over Circular Aperture';
            labelColumn = '\gamma';
            labelRow = 'Separation (in units of r0)';
            labelZ = 'Power';
            labelLegend = obj.buildLegendCell();
            pwrStruct.info = struct('title', tit, ...
                'labelColumn', labelColumn, 'labelRow', labelRow, ...
                'labelZ', labelZ, 'labelLegend', labelLegend);
        end
        function siStruct = fillScintillationIndexOnCircularApertureMetaData(obj)
            siStruct = struct;
            siStruct.columnParams = obj.simulationParameters.gammaStrength;
            siStruct.rowParams = obj.simulationParameters.transverseSeparationInR0Units;
            
            tit = 'Scintillation Index on Circular Aperture';
            labelColumn = '\gamma';
            labelRow = 'Separation (in units of r0)';
            labelZ = 'SI';
            labelLegend = obj.buildLegendCell();
            siStruct.info = struct('title', tit, ...
                'labelColumn', labelColumn, 'labelRow', labelRow, ...
                'labelZ', labelZ, 'labelLegend', labelLegend);
        end
        function mmStruct = fillModeOverlapMetaData(obj)
            mmStruct = struct;
            mmStruct.columnParams = obj.simulationParameters.hermiteGaussOrders;
            mmStruct.rowParams = obj.simulationParameters.hermiteGaussOrders;
            
            iGamma = obj.simulationParameters.gammaCurrentIndex;

            tit = sprintf('HG mode-match, gamma = %3.3g', ...
                obj.simulationParameters.gammaStrength(iGamma));
            labelColumn = 'Transmitted Mode';
            labelRow = 'Reference Mode';
            labelZ = 'Mode-Matching';
            tickX = Util.getHGOrderLabel(obj.simulationParameters.hermiteGaussOrders);
            tickY = Util.getHGOrderLabel(obj.simulationParameters.hermiteGaussOrders);

            mmStruct.info = struct('title', tit, ...
                'labelColumn', labelColumn, 'labelRow', labelRow, ...
                'labelZ', labelZ, 'tickX', {tickX}, 'tickY', {tickY});
        end
        function pwr = getPowerOverCircularAperture(obj, irradiance, apertureRadius)
            circ = obj.simulationParameters.getCircularApertureArray(apertureRadius);
            circ = repmat(circ, [1 1 size(irradiance, 3)]);
            pwr = sum(sum(circ .* irradiance, 2), 1);
        end
        function leg = buildLegendCell(obj)
            sep = obj.simulationParameters.transverseSeparationInR0Units;
            leg = cell(1, length(sep));
            for i = 1 : length(leg)
                leg{i} = sprintf('%2.2g r0',sep(i));
            end
        end
        function params = getSimulationParameters(obj)
            simParams = fieldnames(obj.simulationParameters);
            for p = 1 : numel(simParams)
                params.(simParams{p}) = obj.simulationParameters.(simParams{p});
            end
        end
        function setFreeSpaceConditions(obj)
            if ~ismember(0, obj.simulationParameters.gammaStrength)
                obj.simulationParameters.gammaStrength = ...
                [0, obj.simulationParameters.gammaStrength];
            end
            obj.simulationParameters.gammaCurrentIndex = 1;
        end
    end
    methods(Static)

    end
end