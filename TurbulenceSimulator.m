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
        function outputField = propagate(obj, phScreen, wvl)
            numOrders = numel(obj.simulationParameters.hermiteGaussOrders);
            [Nx, Ny] = obj.simulationParameters.getTransverseGridSize;
            outputField = zeros(Ny, Nx, numOrders);
            sg = obj.simulationParameters.getSuperGaussianFilter;
            sg = repmat(sg, [1, 1, obj.simulationParameters.numberOfPhasePlanes]);
            z = obj.simulationParameters.planePositions;
            delta1 = obj.simulationParameters.gridSpacingSourcePlane;
            deltan = obj.simulationParameters.gridSpacingObservationPlane;
            
            for i = 1 : numOrders
                Uin = obj.simulationParameters.getInputField(i);
                [~,~, outputField(:,:,i)] = ang_spec_multi_prop(Uin, wvl, ...
                    delta1, deltan, z, sg.*exp(1i*phScreen));
            end
        end

        function outputField = propagateInputField(obj, Uin, phScreen, wvl)
            sg = obj.simulationParameters.getSuperGaussianFilter;
            sg = repmat(sg, [1, 1, obj.simulationParameters.numberOfPhasePlanes]);
            z = obj.simulationParameters.planePositions;
            delta1 = obj.simulationParameters.gridSpacingSourcePlane;
            deltan = obj.simulationParameters.gridSpacingObservationPlane;
               
            [~,~, outputField] = ang_spec_multi_prop(Uin, wvl, ...
                delta1, deltan, z, sg.*exp(1i*phScreen));
        end

        function fieldSep = getFieldForEachTransverseSeparation(obj)
            nSep = obj.numberOfTransverseSeparations;
            [Nx, Ny] = obj.simulationParameters.getTransverseGridSize;
            wvl = obj.simulationParameters.wavelength;
            if obj.simulationParameters.isFourthOrder
                wvl = wvl/2;
            end 

            fieldSep = zeros(Ny, Nx, nSep);
            
            for iSep = 1 : nSep
                obj.isAborted = UserInput.isAborted(obj.abortButtonHandle);
                if obj.isAborted
                    break;
                end
                phScreen = obj.inversionAndDisplacementOperationsOnScreen(iSep);
                
                fieldSep(:,:,iSep) = obj.propagate(phScreen, wvl);
            end
        end

        function outMode = getFieldForEachMode(obj)
            phz = obj.inversionAndDisplacementOperationsOnScreen();
            [Nx, Ny] = obj.simulationParameters.getTransverseGridSize;
            wvl = obj.simulationParameters.wavelength;
            if obj.simulationParameters.isFourthOrder
                wvl = wvl/2;
            end 

            obj.isAborted = UserInput.isAborted(obj.abortButtonHandle);
            phScreen = Util.crop(phz, Nx, Ny);
            outMode = obj.propagate(phScreen, wvl);
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

            intGamma.info.date = [datestr(date, 'yyyy-mm-dd'), '_', datestr(clock, 'HHMMSS')];
        end

        function [psiParity, psiChi2] = getPsiDifferenceParityAndChiSquared(obj)
            nGamma = length(obj.simulationParameters.gammaStrength);
            if ~ismember(0, obj.simulationParameters.transverseSeparationInR0Units)
                obj.simulationParameters.transverseSeparationInR0Units = [0,  obj.simulationParameters.transverseSeparationInR0Units];
                obj.numberOfTransverseSeparations = obj.numberOfTransverseSeparations + 1;
            end

            nSep = obj.numberOfTransverseSeparations;

            psiParity = obj.fillPsiParityMetaData();
            psiParity.values = zeros(nSep,nGamma);
            psiParity.params = obj.getSimulationParameters();

            psiChi2 = obj.fillPsiChi2MetaData();
            psiChi2.values = zeros(nSep,nGamma);
            psiChi2.params = obj.getSimulationParameters();

             for iGamma = 1 : nGamma
                if obj.isAborted
                    break;
                end
                UserInput.printOutProgress('Turbulence strength', ...
                    iGamma, nGamma);
                obj.simulationParameters.gammaCurrentIndex = iGamma;
                [psiParity.values(:,iGamma), psiChi2.values(:,iGamma)] = obj.getAveragePsiParityAndChi2();
            end

            psiParity.info.date = [datestr(date, 'yyyy-mm-dd'), '_', datestr(clock, 'HHMMSS')];
            psiChi2.info.date = [datestr(date, 'yyyy-mm-dd'), '_', datestr(clock, 'HHMMSS')];

        end

        function [avgPrt, avgChi2] = getAveragePsiParityAndChi2(obj)
            obj.abortButtonHandle = UserInput.createWaitBar;
            try
                nRe = obj.getNumberOfRealizations;
                NROI = round(1/4*obj.simulationParameters.regionOfInterestAtObservationPlane / obj.simulationParameters.gridSpacingObservationPlane);
                wvl = obj.simulationParameters.wavelength;
                [Nx, Ny] = obj.simulationParameters.getTransverseGridSize;
                
                nSep = obj.numberOfTransverseSeparations;
                prt = zeros(nSep, nRe);
                chi2 = zeros(nSep, nRe);

                vacuumPhase = zeros(Ny, Nx, nSep);

                for iSep = 1 : nSep
                    Uvac = obj.simulationParameters.getInputPointSource(iSep);
                    vacuumPhase(:,:,iSep) = angle(obj.propagateInputField(Uvac, 0, wvl));
                end
                
                for iRe = 1 : nRe
                    obj.isAborted = UserInput.isAborted(obj.abortButtonHandle);
                    if obj.isAborted
                        break;
                    end

                    r0 = obj.simulationParameters.totalFriedCoherenceRadiusByStrength;
                    if isinf(r0(obj.simulationParameters.gammaCurrentIndex))
                        nSep = 1;
                    else
                        nSep = obj.numberOfTransverseSeparations;
                    end

                    phScreen = generateScreen(obj.simulationParameters);
                    Uout0 = obj.propagateInputField(obj.simulationParameters.getInputPointSource(1), phScreen, wvl);
                    Uout0 = Uout0 .* exp(-1i*vacuumPhase(:,:,1));

                    for iSep = 1 : nSep
                        
                        if iSep == 1
                            Uout = Uout0;
                        else
                            Uin = obj.simulationParameters.getInputPointSource(iSep);
                            Uout = obj.propagateInputField(Uin, phScreen, wvl);
                            Uout = Uout .* exp(-1i*vacuumPhase(:,:,iSep));
                        end

                        if obj.simulationParameters.isFourthOrder && obj.simulationParameters.isInverted
                            Uout = Uout0 .* Util.rot90All(Uout,2);
                        elseif obj.simulationParameters.isFourthOrder && ~obj.simulationParameters.isInverted
                            Uout = Uout0 .^ Uout;
                        end 

                        Ucrop = Util.crop(Uout, NROI, NROI);

                        psiOdd = (Ucrop - Util.rot90All(Ucrop,2))/2;
                        psiEven = (Ucrop + Util.rot90All(Ucrop,2))/2;
                        prt(iSep,iRe) = sum( abs(psiOdd(:)).^2)/sum( abs(psiEven(:)).^2);

                        % psiVar = log(Ucrop);
                        chi2(iSep,iRe) = std(Ucrop(:));

                        if isinf(r0(obj.simulationParameters.gammaCurrentIndex))
                            prt(:,iRe) = prt(iSep,iRe);
                            chi2(:,iRe) = chi2(iSep,iRe);
                        end
                    end
                    UserInput.updateWaitBar(obj.abortButtonHandle, iRe, nRe);
                end
                avgPrt = mean(prt, 2);
                avgChi2 = mean(chi2, 2);

                delete(obj.abortButtonHandle);
                obj.abortButtonHandle = [];
                
            catch exception
                if ~isempty(obj.abortButtonHandle)
                    delete(obj.abortButtonHandle);
                    obj.abortButtonHandle = [];
                end
                rethrow(exception);
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

            pwrGamma.info.date = [datestr(date, 'yyyy-mm-dd'), '_', datestr(clock, 'HHMMSS')];
            scintIdx.info.date = [datestr(date, 'yyyy-mm-dd'), '_', datestr(clock, 'HHMMSS')];
            pwrAndSI = {pwrGamma, scintIdx};
        end

        function modeOverlap = getModeMatching(obj)
            obj.numberOfTransverseSeparations = 1;
            [Nx, Ny] = obj.simulationParameters.getTransverseGridSize();
            numOrders = numel(obj.simulationParameters.hermiteGaussOrders);
            if numOrders == 1
                fprintf('WARNING: Mode-match simulation is about to be performed with just a single mode.');
                beep;
            end
            nGamma = length(obj.simulationParameters.gammaStrength);
            obj.setFreeSpaceConditions();
            if  length(obj.simulationParameters.gammaStrength)<2
                error('turbuSimulator:modeMatching', 'This simulation requires a non-zero turbulence strength as input parameter.');
            end

            modeOverlap = cell(nGamma,1);
            refModes = zeros(Ny, Nx, numOrders);

            for i = 1 : numOrders
                refModes(:,:,i) = obj.simulationParameters.getOutputField(i);
            end

            for iGamma = 1 : nGamma
                if obj.isAborted
                    break;
                end
                 UserInput.printOutProgress('Turbulence strength', ...
                    iGamma, nGamma);
                obj.simulationParameters.gammaCurrentIndex = iGamma;
                modeOverlap{iGamma} = obj.fillModeOverlapMetaData;
                modeOverlap{iGamma}.params = obj.getSimulationParameters;

                modeOverlap{iGamma}.values = obj.getModeMatchingAveragedOverRealizations(refModes);
                modeOverlap{iGamma}.info.date = [datestr(date, 'yyyy-mm-dd'), '_', datestr(clock, 'HHMMSS')];
            end
        end

        function modeParity = getModeParity(obj)
            obj.numberOfTransverseSeparations = 1;
            numOrders = numel(obj.simulationParameters.hermiteGaussOrders);
            if numOrders == 1
                fprintf('WARNING: Parity simulation is about to be performed with just a single mode.');
                beep;
            end
            nGamma = length(obj.simulationParameters.gammaStrength);
            
            modeParity = cell(nGamma,1);

            for iGamma = 1 : nGamma
                if obj.isAborted
                    break;
                end
                 UserInput.printOutProgress('Turbulence strength', ...
                    iGamma, nGamma);
                obj.simulationParameters.gammaCurrentIndex = iGamma;
                modeParity{iGamma} = obj.fillModeParityMetaData;
                modeParity{iGamma}.params = obj.getSimulationParameters;

                modeParity{iGamma}.values = obj.getParityAveragedOverRealizations;
                modeParity{iGamma}.info.date = [datestr(date, 'yyyy-mm-dd'), '_', datestr(clock, 'HHMMSS')];
            end
        end

        function params = getSimulationParameters(obj)
            simParams = fieldnames(obj.simulationParameters);
            for p = 1 : numel(simParams)
                params.(simParams{p}) = obj.simulationParameters.(simParams{p});
            end
        end
    end

    methods(Access = private)
        function phScreen = inversionAndDisplacementOperationsOnScreen(obj, varargin)
            idxGamma = obj.simulationParameters.gammaCurrentIndex;
            r0 = obj.simulationParameters.totalFriedCoherenceRadiusByStrength;

            phScreen = obj.phaseScreenProfiles;
            [Nx, Ny] = obj.simulationParameters.getTransverseGridSize;
            [NxEff, ~] = obj.simulationParameters.getPhaseScreenGridSize;
            phScreen = Util.crop(phScreen,NxEff, Ny);

            if isinf(r0(idxGamma)) || ~(obj.simulationParameters.isFourthOrder)
                phScreen = Util.crop(phScreen, Nx, Ny);
                return;
            end

            separation = 0;
            if ~isempty(varargin)
                if isfloat(varargin{1})
                    separation = round(varargin{1});
                end
            end

            phScreen2 = phScreen;
            if obj.simulationParameters.isInverted
                phScreen2 = Util.rot90All(phScreen2,2);
            end

            if separation
                [deltaX, ~] = obj.simulationParameters.getTransverseSeparationInPixels(separation);
                phScreen = Util.displace(phScreen, round(deltaX/2),0);
                phScreen2 = Util.displace(phScreen2, round(deltaX/2),0);
            end
                
           phScreen = Util.crop(phScreen + phScreen2, Nx, Ny);
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
                scintIdx = pwr2/nRe .* pwr.^-2 - 1;
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

        function prt = getParityAveragedOverRealizations(obj)
            obj.abortButtonHandle = UserInput.createWaitBar;

            numOrders = numel(obj.simulationParameters.hermiteGaussOrders);
            nRe = obj.getNumberOfRealizations;
            prt = zeros(numOrders,2);

            try
                for iRe = 1 : nRe
                    obj.isAborted = UserInput.isAborted(obj.abortButtonHandle);
                    if obj.isAborted
                        break;
                    end
                    obj.phaseScreenProfiles = generateScreen(obj.simulationParameters);
                    outputModes = obj.getFieldForEachMode();
                    outputModes = Util.normalize(outputModes);

                    prt = prt + Util.getParityComponents(outputModes);
                    UserInput.updateWaitBar(obj.abortButtonHandle, iRe, nRe);
                end
                delete(obj.abortButtonHandle);
                obj.abortButtonHandle = [];
                prt = prt / nRe;
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
            pwrStruct.columnParams = simParams.structureConstantSquared;
            pwrStruct.rowParams = simParams.transverseSeparationInR0Units;
            
            tit = 'Power Over Circular Aperture';
            labelColumn = 'C_n^2';
            labelRow = 'Separation (in units of r0)';
            labelZ = 'Power';
            labelLegend = obj.buildLegendCell();
            pwrStruct.info = struct('title', tit, ...
                'labelColumn', labelColumn, 'labelRow', labelRow, ...
                'labelZ', labelZ, 'labelLegend', {labelLegend});
        end
        function siStruct = fillScintillationIndexOnCircularApertureMetaData(obj)
            siStruct = struct;
            siStruct.columnParams = obj.simulationParameters.structureConstantSquared;
            siStruct.rowParams = obj.simulationParameters.transverseSeparationInR0Units;
            
            tit = 'Scintillation Index on Circular Aperture';
            labelColumn = 'C_n^2';
            labelRow = 'Separation (in units of r0)';
            labelZ = 'SI';
            labelLegend = obj.buildLegendCell();
            siStruct.info = struct('title', tit, ...
                'labelColumn', labelColumn, 'labelRow', labelRow, ...
                'labelZ', labelZ, 'labelLegend', {labelLegend});
        end

        function psiParity = fillPsiParityMetaData(obj)
            psiParity = struct;
            psiParity.columnParams = obj.simulationParameters.structureConstantSquared;
            psiParity.rowParams = obj.simulationParameters.transverseSeparationInR0Units;
            
            tit = 'Parity ratio of the propagated log irradiance';
            labelColumn = 'C_n^2';
            labelRow = 'Separation (in units of r0)';
            labelZ = 'Parity ratio';
            labelLegend = obj.buildLegendCell();
            psiParity.info = struct('title', tit, ...
                'labelColumn', labelColumn, 'labelRow', labelRow, ...
                'labelZ', labelZ, 'labelLegend', {labelLegend});
        end

        function psiParity = fillPsiChi2MetaData(obj)
            psiParity = struct;
            psiParity.columnParams = obj.simulationParameters.gammaStrength;
            psiParity.rowParams = obj.simulationParameters.transverseSeparationInR0Units;
            
            tit = 'Chi-squared of the propagated log irradiance';
            labelColumn = '\gamma';
            labelRow = 'Separation (in units of r0)';
            labelZ = 'Chi2';
            labelLegend = obj.buildLegendCell();
            psiParity.info = struct('title', tit, ...
                'labelColumn', labelColumn, 'labelRow', labelRow, ...
                'labelZ', labelZ, 'labelLegend', {labelLegend});
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

        function prtStruct = fillModeParityMetaData(obj)
            prtStruct = struct;
            prtStruct.columnParams = [0,1];
            prtStruct.rowParams = obj.simulationParameters.hermiteGaussOrders;
            
            iGamma = obj.simulationParameters.gammaCurrentIndex;

            tit = sprintf('Parity Component, gamma = %3.3g', ...
                obj.simulationParameters.gammaStrength(iGamma));
            labelColumn = 'Transmitted Mode Parity';
            labelRow = 'Collected Parity';
            labelZ = 'Parity';
            tickX = {'Even', 'Odd'};
            tickY = {'Even', 'Odd'};

            prtStruct.info = struct('title', tit, ...
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