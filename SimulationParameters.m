classdef SimulationParameters<handle
    %   Distance unit: m
    
    properties(Access = public)
    	% USER INPUT
		% Geometry
        isFourthOrder;
        propagationDistance;
        numberOfPhasePlanes;
		turbulenceRegionStartPosition;
		turbulenceRegionEndPosition;
        % Turbulence Statistics
        gammaStrength;
        gammaIndex;
        outerScale;
        innerScale;
        % Beam
        wavelength;
        waistAtObservationPlane;    % e-2 intensity radius
        % Simulation
        numberOfRealizations;
        transverseGridSize;
        gridSpacingSourcePlane;
        gridSpacingObservationPlane;
		transverseSeparationInR0Units;
    end
    
    properties(GetAccess = public, SetAccess = private)
		% DERIVED PARAMETERS
		waveNumber;
		planePositions;
		gridSpacingVector;
		friedCoherenceRadiusMatrix; %{i,j} -> Turb. strength, Prop. plane
        totalFriedCoherenceRadiusByStrength;
        regionOfInterestAtSourcePlane;
        regionOfInterestAtObservationPlane;
    end
    
    properties(Access = private)
       complexBeamParameter; 
    end
       
    methods(Access = public)
        function simParams = SimulationParameters(varargin)
           simParams.gammaIndex = 1;
           simParams.checkIfFourthOrder(varargin);
           fData = simParams.openParametersFile();
           simParams.readParameters(fData);
           simParams.setDefaultValueForBlankParameters();
           simParams.computeDerivedQuantities();
        end
        function sg = getSuperGaussianFilter(obj)
            [x1, y1] = obj.getMeshGridAtSourcePlane();
            [NxEff, NyEff] = obj.getEffectiveGridSize;
            sg = exp(-(x1/(0.47*NxEff*delta1)).^16 ... 
                -(y1/(0.47*NyEff*delta1)).^16);
        end
        function Uin = getInputField(obj)
            [x1,y1] = obj.getMeshGridAtSourcePlane();
            qz = obj.complexBeamParameter;
            k = obj.waveNumber;
            Uin = 1/qz*exp(1i*k/(2*qz)*(x1.^2+y1.^2)); 
        end
        function [NxEff, NyEff] = getEffectiveGridSize(obj)
            [Nx, Ny] = obj.getTransverseGridSize;
            r0sw = obj.totalFriedCoherenceRadiusByStrength;
            r0sw(isinf(r0sw)) = 0;
            r0sw = r0sw(obj.gammaIndex);
            maxSep = max(obj.transverseSeparationInR0Units);
            extraGridLength = maxSep * r0sw/ min(obj.gridSpacingVector);
            
            NxEff = Nx + extraGridLength;
            NyEff = Ny; % Transverse separation only in x direction for now.
            
            % Get smaller power of 2 numbers that exceeds NxEff and NyEff
            NxEff = 2.^(ceil(log2(NxEff)));
            NyEff = 2.^(ceil(log2(NyEff)));            
        end
        function [Nx, Ny] = getTransverseGridSize(obj)
            Nx = obj.transverseGridSize;
            Ny = obj.transverseGridSize;
        end
        function [x1, y1] = getMeshGridAtSourcePlane(obj)
            [NxEff, NyEff] = obj.getEffectiveGridSize;
            delta1 = obj.gridSpacingSourcePlane;
            [x1,y1] = meshgrid((-NxEff/2 : NxEff/2-1) * delta1, ...
                (-NyEff/2 : NyEff/2-1) * delta1);
        end
		function fail = constraintAnalysis(obj)
		modelSensitivity = 4; % see pag. 173
        
        L = obj.propagationDistance;
        wvl = obj.wavelength;
        wn = obj.waistAtObservationPlane;
        k = obj.waveNumber;
        D1 = obj.regionOfInterestAtSourcePlane;
        D2 = obj.regionOfInterestAtObservationPlane;
        r0sw = min(obj.totalFriedCoherenceRadiusByStrength);
        
		rad = L*(1 + (k*wn^2/(2*L))^2);	% Beam radius of curvature
		
		% Effective ROI
		D1p = D1 + modelSensitivity*wvl*L/r0sw;
		D2p = D2 + modelSensitivity*wvl*L/r0sw;
		
		figHndl = figure;
		params = struct('figureHandle', figHndl, ...
			'radiusOfCurvature', rad, ...
			'effectiveSourceROI', D1p, 'effectiveObsROI', D2p);
		
        obj.plotConstraint2(params);
        obj.plotConstraint1(params);
        obj.plotConstraint3(params);
        obj.plotCurrentValue();
        
        fprintf('Veryfing sampling requirements...\n');
		fail1 = obj.checkConstraint1(params);
		fail2 = obj.checkConstraint2(params);
		fail3 = obj.checkConstraint3(params);
		fail4 = obj.checkConstraint4();
		
		fail = (fail1 || fail2 || fail3 || fail4);
        end

    end
      
    methods(Access = private)
        function checkIfFourthOrder(obj,isFourthOrder)
            if (numel(isFourthOrder) == 0)
                obj.isFourthOrder = 0;
            elseif (numel(isFourthOrder) == 1)
                obj.isFourthOrder = isFourthOrder;
            else
                if (strcmpi('fourthorder',isFourthOrder{1}))
                    obj.isFourthOrder = isFourthOrder{2};
                end
            end
        end
        function readParameters(obj,fData)
            hashTable = regexp(fData, '(\w+)\s*:\s*([\w\.\-\,\s]+)$', 'tokens', 'lineanchors');
            for i = 1 : length(hashTable)
                if (length(hashTable{i}) == 2)
                    splitStr = regexprep(hashTable{i}{2},'\s+', '');
                    splitStr = regexp(splitStr, '[\,]+', 'split');
                    obj.(hashTable{i}{1}) = sort(str2double(splitStr));
                end
            end
        end
        function setDefaultValueForBlankParameters(obj)
            z0 = obj.turbulenceRegionStartPosition;
            z1 = obj.turbulenceRegionEndPosition;
            transvSep = obj.transverseSeparationInR0Units;
            
           if (isnan(z0) || isempty(z0))
              obj.turbulenceRegionStartPosition = 0;
           end
           if (isnan(z1) || isempty(z1))
              obj.turbulenceRegionEndPosition = obj.propagationDistance;
           end
           if (isempty(transvSep))
               obj.transverseSeparationInR0Units = 0;
           end
        end
        function computeDerivedQuantities(obj)
           wvl = obj.wavelength;
           L = obj.propagationDistance;
           npl = obj.numberOfPhasePlanes;
           delta1 = obj.gridSpacingSourcePlane;
           deltan = obj.gridSpacingObservationPlane;
           wn = obj.waistAtObservationPlane;
           
           obj.waveNumber =  2*pi/wvl;
           k = obj.waveNumber;
           obj.planePositions = linspace(0, L, npl);
           z = obj.planePositions;
           
           obj.complexBeamParameter = -L -1i*k*wn^2/2;
           
           obj.gridSpacingVector = (1-z/L)*delta1 + z/L*deltan;
           w1 = wn*sqrt(1 + 2*L/ (k*wn^2));
           obj.regionOfInterestAtSourcePlane = 4*w1;
           obj.regionOfInterestAtObservationPlane = 4*wn;
           
           obj.computeFriedCoherenceRadiusMatrix();
        end
        function computeFriedCoherenceRadiusMatrix(obj)
            L = obj.propagationDistance;
            z = obj.planePositions;
            npl = obj.numberOfPhasePlanes;
            zmin = obj.turbulenceRegionStartPosition;
            zmax = obj.turbulenceRegionEndPosition;
            g = obj.gammaStrength';
            if (obj.isFourthOrder)
                k = obj.waveNumber/2;
            else
                k = obj.waveNumber;
            end
            
            zmask = (z > zmin & z < zmax);
            ztmin_idx = find(zmask, 1, 'first');
            ztmax_idx = find(zmask, 1, 'last');
            s = sum((1-z(ztmin_idx:ztmax_idx)/L).^(5/3));
            
            r0 = (0.423 * k^2/(7.75 * L^(5/3) * s) * g.^2).^(-3/5);
            r0 = repmat(r0, 1, npl);
            r0 = r0 .* repmat(zmask, length(g), 1);
            r0(isinf(1./r0)) = inf;
            r0(isnan(r0)) = inf;
            
            obj.friedCoherenceRadiusMatrix = r0;
            
            aux_mat = repmat( (z .* zmask)/L,length(g),1);
            obj.totalFriedCoherenceRadiusByStrength = ...
                sum(r0.^(-5/3) .* aux_mat.^(5/3), 2).^(-3/5);
        end
        function plotConstraint1(obj, params)
            L = obj.propagationDistance;
            wvl = obj.wavelength;
            D1p = params.effectiveSourceROI;
            D2p = params.effectiveObsROI;
            
            d1 = linspace(0, 1.1*wvl*L/D2p, 100);
            dn = linspace(0, 1.1*wvl*L/D1p, 100);
            [d1, dn] = meshgrid(d1, dn);
            
            deltan_max = -D2p/D1p*d1 + wvl*L/D1p;
            plot(d1(1,:), deltan_max(1,:), 'k--', 'Linewidth', 2);
            axis([0 d1(end) 0 dn(end)]);
            set(gca, 'Color', 'none', 'Layer', 'top');
        end
        function plotConstraint2(obj, params)
            L = obj.propagationDistance;
            wvl = obj.wavelength;
            D1p = params.effectiveSourceROI;
            D2p = params.effectiveObsROI;
            
            d1 = linspace(0, 1.1*wvl*L/D2p, 100);
            dn = linspace(0, 1.1*wvl*L/D1p, 100);
            [d1, dn] = meshgrid(d1, dn);
            
            N2 = log2((wvl * L + D1p*dn + D2p*d1) ./ (2 * d1 .* dn));
            contourf(d1, dn, N2);
            %clabel(C,hh, 'FontSize', 15, 'Rotation', 0, ...
            %    'FontWeight', 'bold');
            xlabel('\delta_1 [m]');
            ylabel('\delta_n [m]');
            colorbar;
            hold all;
        end
        function plotConstraint3(obj, params)
            L = obj.propagationDistance;
            wvl = obj.wavelength;
            D1 = obj.regionOfInterestAtSourcePlane;
            
            rad = params.radiusOfCurvature;
            D1p = params.effectiveSourceROI;
            D2p = params.effectiveObsROI;
            
            d1 = linspace(0, 1.1*wvl*L/D2p, 100);
            dn = linspace(0, 1.1*wvl*L/D1p, 100);
            [d1, ~] = meshgrid(d1, dn);
            
            dnmin3 = (1+L/rad)*d1 - wvl*L/D1;
            dnmax3 = (1+L/rad)*d1 + wvl*L/D1;
            plot(d1(1,:), dnmax3(1,:), 'k-.');
            set(gca, 'Color', 'none', 'Layer', 'top');
            plot(d1(1,:), dnmin3(1,:), 'k-.');
            set(gca, 'Color', 'none', 'Layer', 'top');    
        end
        function plotCurrentValue(obj)
            delta1 = obj.gridSpacingSourcePlane;
            deltan = obj.gridSpacingObservationPlane;
            plot(delta1, deltan, 'w*', 'MarkerSize', 20);
            set(gca, 'Color', 'none', 'Layer', 'top'); 
        end
        function fail1 = checkConstraint1(obj,params)
            wvl = obj.wavelength;
            L = obj.propagationDistance;
            deltan = obj.gridSpacingObservationPlane;
            delta1 = obj.gridSpacingSourcePlane;
            
            D1p = params.effectiveSourceROI;
            D2p = params.effectiveObsROI;
            
            fail1 = (deltan >= -D2p/D1p*delta1 + wvl*L/D1p);
            fprintf('Constraint 1: ');
            if ~fail1
                fprintf('Satisfied\n');
            else
                fprintf('Not satisfied [deltan = %3.2e should be smaller than %3.2e]\n', ...
                    deltan, -D2p/D1p*delta1 + wvl*L/D1p);
            end
        end
        function fail2 = checkConstraint2(obj,params)
            wvl = obj.wavelength;
            L = obj.propagationDistance;
            [Nx, ~] = obj.getEffectiveGridSize;
            N = min(Nx);
            deltan = obj.gridSpacingObservationPlane;
            delta1 = obj.gridSpacingSourcePlane;
            
            D1p = params.effectiveSourceROI;
            D2p = params.effectiveObsROI;
            
            Nmin = (wvl*L + D1p*deltan + D2p*delta1)./ (2*delta1 .* deltan);

            fail2 = (N <= Nmin);
            fprintf('Constraint 2: ');
            if ~fail2
                fprintf('Satisfied\n');
            else
                fprintf('Not satisfied [N = %u should be greater than %u]\n', ...
                    N, round(Nmin));
            end
        end
        function fail3 = checkConstraint3(obj,params)
            wvl = obj.wavelength;
            L = obj.propagationDistance;
            deltan = obj.gridSpacingObservationPlane;
            delta1 = obj.gridSpacingSourcePlane;
            D1 = obj.regionOfInterestAtSourcePlane;
            
            rad = params.radiusOfCurvature;
            
            success3 = (deltan > (1+L/rad)*delta1 - wvl*L/D1) & ...
                (deltan < (1+L/rad)*delta1 + wvl*L/D1);

            fail3 = ~success3;
            fprintf('Constraint 3: ');
            if ~fail3
                fprintf('Satisfied\n');
            else
                fprintf(['Not satisfied [deltan = %3.2e should be between ', ...
                    '%3.2e and %3.2e]\n'], deltan, ...
                    (1+L/rad)*delta1 - wvl*L/D1, (1+L/rad)*delta1 + wvl*L/D1);
            end
        end
        function fail4 = checkConstraint4(obj)
            wvl = obj.wavelength;
            [Nx, ~] = obj.getEffectiveGridSize;
            N = min(Nx);
            deltan = obj.gridSpacingObservationPlane;
            delta1 = obj.gridSpacingSourcePlane;
            delta = obj.gridSpacingVector;
            
            zmax = min([delta1 deltan])^2 * N / wvl;
            
            fail4 = (max(delta) >= zmax);
            fprintf('Constraint 4: ');
            if ~fail4
                fprintf('Satisfied\n');
            else
                fprintf(['Not satisfied [Max(delta) = %3.2e ', ...
                    'should be smaller than %3.2e]\n'], ...
                    max(delta), zmax);
            end
        end
    end
        
    methods(Static)
        function fData = openParametersFile()
            fileName = 'inputParameters.dat';
            fileName = uigetfile(fileName);
            fid = fopen(fileName);
            try
                fData = fread(fid, inf, '*char');
                fData = fData';
            catch exception
                fprintf('Error reading inputParameters.dat.');
                fclose(fid);
                rethrow(exception);
            end
            fclose(fid);
        end
        
    end
end

