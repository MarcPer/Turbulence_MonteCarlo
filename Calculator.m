classdef Calculator
    %Calculator Computes all relevant quantities
    %   All class methods take either field or intensity profiles as input
    %   arguments and compute some other quantity.
    %
    %   Public methods also return information that can be used by the
    %   Plotter class, such as plot labels and title.
    %
    %   The following plot informations are returned:
    %       title
    %       labelX
    %       labelY
    
    properties
    end
    
    methods(Static)
        function pwrSlit = computePowerThroughSlit(intProfile, simParams)
        % computePowerThroughSlit returns power for each gamma and
        % separation.
        %
        % SYNTAX:
        % [pwrSlit, plotInfo] = computePowerThroughSlit(intProfile,
        % simParams)
        %
        % OUTPUT:
        % pwrSlit: Matrix for power through slit (i,j) -> [Separation, Gamma]
        % plotInfo: Structure with plot title and labels
        %
        % INPUT:
        % intProfile: Cell (indexed by Gamma), each element consisting of
        % an array indexed as [profile row, profile column, separation]
        % simParams: Instance of SimulationParameters
            slitWidthPx = round(simParams.slitWidth / ...
                simParams.gridSpacingObservationPlane);
            if (slitWidthPx == 0)
                fprintf('Slit width set to smaller than pixel size.');
                fprintf('Setting to 1 pixel = %3.3g m', simParams.gridSpacingObservationPlane);
                slitWidthPx = 1;
            end
            
            pwrSlit = struct;
            pwrSlit.columnParams = simParams.gammaStrength;
            pwrSlit.rowParams = simParams.transverseSeparationInR0Units;
            
            if simParams.isFourthOrder
                pwrSlit.values = Calculator.coincidenceSlitIntegratePx(intProfile, slitWidthPx);
            else
                pwrSlit.values = Calculator.intensitySlitIntegratePx(intProfile, slitWidthPx);
            end
            
            tit = 'Intensity vs Turbulence Strength';
            labelColumn = '\gamma';
            labelRow = 'Separation (in r0)';
            labelZ = 'Power through slit';
            dt = [datestr(date, 'yyyy-mm-dd'), '_', datestr(clock, 'HHMMSS')];
            pwrSlit.info = struct('title', tit, ...
                'labelColumn', labelColumn, 'labelRow', labelRow, ...
                'labelZ', labelZ, 'date', dt);
            
        end
        function pwrSlit = computePowerThroughCircularAperture(apertureDiameter, intProfile, simParams)
        % computePowerThroughSlit returns power for each gamma and
        % separation.
        %
        % SYNTAX:
        % [pwrSlit, plotInfo] = computePowerThroughSlit(apertureDiameter,
        % intProfile, simParams)
        %
        % OUTPUT:
        % pwrSlit: Matrix for power through slit (i,j) -> [Separation, Gamma]
        % plotInfo: Structure with plot title and labels
        %
        % INPUT:
        % apertureDiameter: Scalar diameter of the circular aperture
        % intProfile: Cell (indexed by Gamma), each element consisting of
        % an array indexed as [profile row, profile column, separation]
        % simParams: Instance of SimulationParameters
            apertureDiamPx = round(apertureDiameter / ...
                simParams.gridSpacingObservationPlane);
            if (apertureDiamPx == 0)
                fprintf('Aperture diameter set to smaller than pixel size.');
                fprintf('Setting to 1 pixel = %3.3g m', simParams.gridSpacingObservationPlane);
                apertureDiamPx = 1;
            end
            
            pwrSlit = struct;
            pwrSlit.columnParams = simParams.gammaStrength;
            pwrSlit.rowParams = simParams.transverseSeparationInR0Units;
            
            if simParams.isFourthOrder
                pwrSlit.values = Calculator.coincidenceSlitIntegratePx(intProfile, apertureDiamPx);
            else
                pwrSlit.values = Calculator.intensitySlitIntegratePx(intProfile, apertureDiamPx);
            end
            
            tit = 'Power over circular aperture vs Turbulence Strength';
            labelColumn = '\gamma';
            labelRow = 'Separation (in r0)';
            labelZ = 'Power through circular aperture';
            dt = [datestr(date, 'yyyy-mm-dd'), '_', datestr(clock, 'HHMMSS')];
            pwrSlit.info = struct('title', tit, ...
                'labelColumn', labelColumn, 'labelRow', labelRow, ...
                'labelZ', labelZ, 'date', dt);
            
        end

        function er = computeErrorRateVsRelativeLengths(data)
            if ~iscell(data)
                error('calculator:computeERr0', 'Input argument should be a cell.');
            end
            nGamma = numel(data{1}.params.gammaStrength);
            er = struct;
            er.info = struct('title', 'Error-Rate', 'labelColumn', '2\sigma_1/r_0', ...
                    'labelZ', 'ER');
            if nGamma < 2
                return;
            end
            er.columnParams = sqrt(2)*data{1}.params.waistAtSourcePlane ...
                ./ data{1}.params.totalFriedCoherenceRadiusByStrength;
            er.values = zeros(nGamma,1);

            for iGamma = 1 : nGamma - 1
                [nOrd, ~] = size(data{iGamma}.values);
                er.values(iGamma+1) = sum( sum( data{iGamma}.values .* ~eye(nOrd) ,2) ,1) ...
                / sum( diag(data{iGamma}.values) );
            end
        end

        function er = computeErrorRateVsCn2(data)
            if ~iscell(data)
                error('calculator:computeERCn2', 'Input argument should be a cell.');
            end
            nGamma = numel(data{1}.params.gammaStrength);
            er = struct;
            er.info = struct('title', 'Error-Rate', 'labelColumn', 'C_n^2', ...
                    'labelZ', 'ER');
            if nGamma < 2
                return;
            end
            er.columnParams = data{1}.params.structureConstantSquared;
            er.values = zeros(nGamma,1);

            for iGamma = 1 : nGamma - 1
                [nOrd, ~] = size(data{iGamma}.values);
                er.values(iGamma+1) = sum( sum( data{iGamma}.values .* ~eye(nOrd) ,2) ,1) ...
                / sum( diag(data{iGamma}.values) );
            end
        end

        function dp = computeDetectionProbabilityVsRelativeLengths(data)
            if ~iscell(data)
                error('calculator:computeDPr0', 'Input argument should be a cell.');
            end
            nGamma = numel(data{1}.params.gammaStrength);
            nOrd = numel(data{1}.info.tickX);

            dp = struct;
            dp.info = struct('title', 'Error-Rate', 'labelColumn', '2\sigma_1/r_0', ...
                    'labelZ', 'ER');

            leg = data{1}.info.tickX;
            dp.info.labelLegend = leg;

            if nGamma < 2
                return;
            end
            dp.columnParams = sqrt(2)*data{1}.params.waistAtSourcePlane ...
                ./ data{1}.params.totalFriedCoherenceRadiusByStrength;
            dp.values = ones(nGamma,nOrd);

            for iGamma = 1 : nGamma - 1
                [nOrd, ~] = size(data{iGamma}.values);
                dp.values(iGamma+1,:) = sum(data{iGamma}.values ,1);
            end
        end

    end
    
    methods(Static, Access = private)
        function coincSlitPwr = coincidenceSlitIntegratePx(intProfile, slitWidth)
            %CoincSlitIntegratePx Integrates array over effective slit that is the
            %   convolution of a slit of width 'a' (in pixels) with itself.
            %
            %   SYNTAX:
            %   y = CoincSlitIntegrate(A,a);
            %
            %   y = CoincSlitIntegrate(A,a) returns the result of an integration of
            %   matrix 'A' over the effective slit that is given by the convolution of
            %   a rectangular slit of width 'a' with itself. The resulting amplitude
            %   apperture takes the form of a triangle of width '2a' and height 1 in
            %   the horizontal direction.
            %
            %   The slit is automatically horizontally displaced so that its horizontal
            %   maximum coincides with the maximum element of 'A'.
            %
            %   This function is recurrent in coincidence signals integrated over a
            %   slit of width 'a'.
            
            if ~iscell(intProfile)
                intProfile = {intProfile};
            end
            
            numCells = length(intProfile);
            [hght, wdth, numSeps] = Calculator.getSize(intProfile);
            [~, colMax] = Calculator.getIndexForMaxima(intProfile);
            
            coincSlitPwr = zeros(max(numSeps),numCells);
            
            for iCell = 1 : numCells
                cSlit = triang(2*slitWidth -1)';
                lenSlit = length(cSlit);
                
                % Adjust the slit function to be of the same size as matrix 'A'
                if (lenSlit > wdth(iCell))
                    lenDiff = lenSlit - wdth(iCell);
                    cSlit([1:floor(lenDiff/2), lenSlit-floor(lenDiff/2)+1:end]) = [];
                else
                    cSlit = [cSlit, zeros(1, wdth(iCell)-lenSlit)];
                end
                
                % Center slit function horizontally on maximum element of 'A'
                [~, xSlit] = max(cSlit);
                cSlit = repmat(cSlit, [hght(iCell), 1, numSeps(iCell)]);
                
                cSlit = Util.displace(cSlit, colMax(iCell,:) - xSlit, 0);
                                
                coincSlitPwr(:,iCell) = sum( sum( ...
                    cSlit .* intProfile{iCell}, 2), 1);
                
            end
            % Normalization
            c0 = coincSlitPwr(:,1);
            c0 = repmat(c0, [1 numCells]);
            coincSlitPwr = coincSlitPwr ./ c0;
        end

        function singSlitPwr = intensitySlitIntegratePx(intProfile, slitWidth)
            %intensitySlitIntegratePx Integrates array over a slit of width 'a' 
            % (in pixels) with itself.
            %
            %   SYNTAX:
            %   y = intensitySlitIntegratePx(A,a);
            %
            %   y = intensitySlitIntegratePx(A,a) returns the result of an integration of
            %   matrix 'A' over a slit of width 'a'.
            %
            %   The slit is automatically horizontally displaced so that its horizontal
            %   maximum coincides with the maximum element of 'A'.
            
            if ~iscell(intProfile)
                intProfile = {intProfile};
            end
            numCells = length(intProfile);
            [hght, wdth, numSeps] = Calculator.getSize(intProfile);
            [~, colMax] = Calculator.getIndexForMaxima(intProfile);
            singSlitPwr = zeros(max(numSeps),numCells);
            
            for iCell = 1 : numCells
                % Adjust the slit function to be of the same size as matrix 'A'
                if (slitWidth >= wdth(iCell))
                    slit = ones(1,wdth(iCell));
                else
                    slit = [ones(1,slitWidth), zeros(1,wdth(iCell)-slitWidth)];
                end
                
                % Center slit function horizontally on maximum element of 'A'
                xSlit = ceil(slitWidth/2);
                
                slit = repmat(slit, [hght(iCell), 1, numSeps(iCell)]);
                slit = Util.displace(slit, colMax(iCell,:)-xSlit, 0);
                                
                singSlitPwr(:,iCell) = sum( sum( ...
                    slit .* intProfile{iCell}, 2), 1);
                
            end
            % Normalization
            s0 = singSlitPwr(:,1);
            s0 = repmat(s0, [1 numCells]);
            singSlitPwr = singSlitPwr ./ s0;
        end 

        function [hght, wdth, numSep] = getSize(A)
            if ~iscell(A)
                A = {A};
            end
            numCells = length(A); 
            
            hght = zeros(numCells, 1);
            wdth = zeros(numCells, 1);
            numSep = zeros(numCells, 1);
            
            for iCell = 1 : numCells
                [hght(iCell), wdth(iCell), numSep(iCell)] = ...
                    size(A{iCell});
            end
        end
        function [rowMax, colMax] = getIndexForMaxima(A)
            if ~iscell(A)
                A = {A};
            end
            numCells = length(A);
            numSeps = length(A{1}(1,1,:));
            
            rowMax = zeros(numCells,numSeps);
            colMax = zeros(numCells,numSeps);
            
            for iCell = 1 : numCells
                for iSep = 1 : numSeps
                    [~, rowMax(iCell,iSep)] = max( max(A{iCell}(:,:,iSep),[],2), [],1);
                    [~, colMax(iCell,iSep)] = max( max(A{iCell}(:,:,iSep),[],1), [],2);
                end
            end
        end



    end
end

