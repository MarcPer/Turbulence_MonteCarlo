classdef SimulationParameters<handle
    %   Distance unit: m
    
    properties(Access = public)
    	% USER INPUT
        % Simulation Information
        simulationType;
        isFourthOrder;
        isInverted;
        shutdownOrHibernate;
		% Geometry
        propagationDistance;
        numberOfPhasePlanes;
		turbulenceRegionStartPosition;
		turbulenceRegionEndPosition;
        slitWidth;
        circularApertureRadius;
        % Turbulence Statistics
        gammaRaw;
        gammaCurrentIndex;
        outerScale;
        innerScale;
        % Beam
        wavelength;
        waistRadius;
        waistPosition;
        hermiteGaussOrders;
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
        waistAtSourcePlane;         % e-2 intensity radius
        waistAtObservationPlane;    % e-2 intensity radius
        gammaStrength;
		planePositions;
		gridSpacingVector;
		friedCoherenceRadiusMatrix;   % {i,j} -> Turb. strength, Prop. plane
        totalFriedCoherenceRadiusByStrength;
        structureConstantSquared;
        regionOfInterestAtSourcePlane;
        regionOfInterestAtObservationPlane;
    end
    
    properties(Access = private)
       complexBeamParameter; 
    end
       
    methods(Access = public)
        function simParams = SimulationParameters(pStruct)
            fieldNames = fieldnames(pStruct);
            for i = 1 : numel(fieldNames)
                simParams.(fieldNames{i}) = pStruct.(fieldNames{i});
            end

            simParams.simulationType = upper(simParams.simulationType);
            simParams.checkValidParameters;
            simParams.setDefaultValueForBlankParameters;
            simParams.computeDerivedQuantities;
        end
        
        function isAbort = checkConstraints(obj)
            isFail = obj.constraintAnalysis;
            isAbort = UserInput.abortWhenConstraintFail(isFail, obj.shutdownOrHibernate); 
        end

        function sg = getSuperGaussianFilter(obj)
            delta1 = obj.gridSpacingSourcePlane;
            [x1, y1] = obj.getMeshGridAtSourcePlane();
            [Nx, Ny] = obj.getTransverseGridSize;
            sg = exp(-(x1/(0.47*Nx*delta1)).^16 ... 
                -(y1/(0.47*Ny*delta1)).^16);
        end

        function Uin = getInputField(obj, varargin)
            [x1,y1] = obj.getMeshGridAtSourcePlane();
            qz = obj.complexBeamParameter;
            k = obj.waveNumber;
            if obj.isFourthOrder
                k = 2*k;
            end

            if isempty(varargin)
                Uin = HG(0,qz,k,x1) .* HG(0,qz,k,y1);
                return;
            end

            n = varargin{1};
            if ~isfloat(n)
                Uin = HG(0,qz,k,x1) .* HG(0,qz,k,y1);
                return;
            end

            hgOrderXY = obj.hermiteGaussOrders;
            [nx, ny] = Util.getHermiteGaussOrders(hgOrderXY(n));
            Uin = HG(nx,qz,k,x1) .* HG(ny,qz,k,y1);
        end

        function Uout = getOutputField(obj, varargin)
            [xL,yL] = obj.getMeshGridAtObservationPlane();
            L = obj.propagationDistance;
            qz = obj.complexBeamParameter + L;
            k = obj.waveNumber;
            if obj.isFourthOrder
                k = 2*k;
            end

            if isempty(varargin)
                Uout = HG(0,qz,k,xL) .* HG(0,qz,k,yL);
                return;
            end

            n = varargin{1};
            if ~isfloat(n)
                Uout = HG(0,qz,k,xL) .* HG(0,qz,k,yL);
                return;
            end

            hgOrderXY = obj.hermiteGaussOrders;
            [nx, ny] = Util.getHermiteGaussOrders(hgOrderXY(n));
            Uout = HG(nx,qz,k,xL) .* HG(ny,qz,k,yL);
        end

        function Uin = getInputPointSource(obj, separationIndex)
            [x1,y1] = obj.getMeshGridAtSourcePlane();
            L = obj.propagationDistance;
            k = obj.waveNumber;
            iGamma = obj.gammaCurrentIndex;
            Dwindow = obj.regionOfInterestAtObservationPlane/2;
            r0 = obj.totalFriedCoherenceRadiusByStrength;
            r0(isinf(r0)) = 0;
            xc = obj.transverseSeparationInR0Units(separationIndex) * r0(iGamma);

            arg = k*Dwindow/(2*pi*L);
            Uin = 2*pi*L/k * exp(-1i*k/(2*L) * (x1.^2 + y1.^2)) .* exp(1i*k/(2*L) * xc^2) .* exp(-1i*k/L * xc * x1) .* ...
                arg^2 .* sinc(arg*x1) .* sinc(arg*y1);
        end

        function vacuumPhase = getOutputVacuumPhaseProfile(obj, separationIndex)
            [xn,yn] = obj.getMeshGridAtObservationPlane;
            L = obj.propagationDistance;
            wvl = obj.wavelength;
            iGamma = obj.gammaCurrentIndex;
            r0 = obj.totalFriedCoherenceRadiusByStrength(iGamma);
            r0(isinf(r0)) = 0;
            xc = obj.transverseSeparationInR0Units(separationIndex) * r0;

            vacuumPhase = exp(-1i*pi/(wvl*L) * ((xn-xc).^2 + yn.^2));
        end

        function [NxEff, NyEff] = getPhaseScreenGridSize(obj)
            [Nx, Ny] = obj.getTransverseGridSize;
            % r0sw = obj.totalFriedCoherenceRadiusByStrength;
            % r0sw(isinf(r0sw)) = 0;
            % r0sw = r0sw(obj.gammaCurrentIndex);
            % maxSep = max(obj.transverseSeparationInR0Units);
            % extraGridLength = maxSep/2 * r0sw/ min(obj.gridSpacingVector);
            
            % NxEff = Nx + extraGridLength;
            % NyEff = Ny; % Transverse separation only in x direction for now.
            
            % % Get smaller power of 2 numbers that exceed NxEff and NyEff
            % NxEff = 2.^(ceil(log2(NxEff)));
            % NyEff = 2.^(ceil(log2(NyEff)));            
            NxEff = Nx;
            NyEff = Ny;
        end

        function [NxMax, NyMax] = getMaximumScreenGridSize(obj)
            [Nx, Ny] = obj.getTransverseGridSize;
            r0sw = obj.totalFriedCoherenceRadiusByStrength;
            r0sw(isinf(r0sw)) = 0;
            maxSep = max(obj.transverseSeparationInR0Units);
            extraGridLength = maxSep/2 * r0sw/ min(obj.gridSpacingVector);
            
            NxMax = Nx + max(extraGridLength);
            NyMax = Ny; % Transverse separation only in x direction for now.
            
            % Get smaller power of 2 numbers that exceed NxMax and NyMax
            NxMax = 2.^(ceil(log2(NxMax)));
            NyMax = 2.^(ceil(log2(NyMax)));            
        end        

        function [D1p, D2p] = getEffectiveROI(obj)
            modelSensitivity = 4; % see pag. 173
            
            L = obj.propagationDistance;
            wvl = obj.wavelength;
            r0sw = obj.totalFriedCoherenceRadiusByStrength;

            [nx, ny] = Util.getMaximumHermiteGaussOrders(obj.hermiteGaussOrders);

            D1 = obj.regionOfInterestAtSourcePlane;
            D2 = obj.regionOfInterestAtObservationPlane;
            
            if all(isnan(r0sw))
                D1p = D1;
                D2p = D2;
            else
                D1p = max(D1 + modelSensitivity*wvl*L ./ r0sw);
                D2p = max(D2 + modelSensitivity*wvl*L ./ r0sw);
            end

            D1p = sqrt(max(nx,ny)+1)*D1p;
            D2p = sqrt(max(nx,ny)+1)*D2p;
        end

        function [deltaX, deltaY] = getTransverseSeparationInPixels(obj, sepIndex)
            r0sw = obj.totalFriedCoherenceRadiusByStrength;
            r0sw(isinf(r0sw)) = 0;
            r0sw = r0sw(obj.gammaCurrentIndex);
            sep = obj.transverseSeparationInR0Units(sepIndex);
            deltaX = round(sep * r0sw ./ obj.gridSpacingVector); 
            deltaY = round(sep * r0sw ./ obj.gridSpacingVector);
        end

        function [Nx, Ny] = getTransverseGridSize(obj)
            Nx = obj.transverseGridSize;
            Ny = obj.transverseGridSize;
        end

        function [x1, y1] = getMeshGridAtSourcePlane(obj)
            [Nx, Ny] = obj.getTransverseGridSize;
            delta1 = obj.gridSpacingSourcePlane;
            [x1,y1] = meshgrid((-Nx/2 : Nx/2-1) * delta1, ...
                (-Ny/2 : Ny/2-1) * delta1);
        end

        function [xn, yn] = getMeshGridAtObservationPlane(obj)
            [Nx, Ny] = obj.getTransverseGridSize;
            deltan = obj.gridSpacingObservationPlane;
            [xn,yn] = meshgrid((-Nx/2 : Nx/2-1) * deltan, ...
                (-Ny/2 : Ny/2-1) * deltan);
        end

        function circ = getCircularApertureArray(obj, apertureRadius)
            [xn, yn] = obj.getMeshGridAtObservationPlane;
            circ = ( xn.^2 + yn.^2 <= apertureRadius^2);
        end

		function fail = constraintAnalysis(obj)
        [D1p, D2p] = obj.getEffectiveROI;
        
        rad = 1 / real(1/obj.complexBeamParameter); % Radius of curvature at source plane
		
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
		
		fail = (fail1 || fail2 || fail3);
        end
        
        function apertureRadius = setPointDetectorAtObservationPlane(obj,isPointDtc)
            apertureRadius = NaN;
            if isPointDtc
                r0 = obj.totalFriedCoherenceRadiusByStrength;
                r0 = ~isinf(r0) .* r0;
                if min(r0) < obj.gridSpacingObservationPlane
                    fprintf('WARNING: Point detector condition not satisfied for strongest turbulence.\n');
                    fprintf('Setting detector area to 1 pixel.\n');
                    beep();
                    apertureRadius = obj.gridSpacingObservationPlane;
                    obj.regionOfInterestAtObservationPlane = 4*apertureRadius;
                    return;
                end
                apertureRadius = min(r0);
                obj.regionOfInterestAtObservationPlane = 4*min(r0);
            end
        end
    end
      
    methods(Access = private)
        function checkValidParameters(obj)
            obj.checkNonNegativeParams;
            obj.checkOutOfRangeParams;
            obj.checkIntegerParams;
        end

        function checkNonNegativeParams(obj)
            % not implemented
        end

        function checkOutOfRangeParams(obj)
            % not implemented
        end

        function checkIntegerParams(obj)
            % not implemented
        end

        function setDefaultValueForBlankParameters(obj)
            obj.setValueIfEmpty('turbulenceRegionStartPosition', 0);
            obj.setValueIfEmpty('turbulenceRegionEndPosition', obj.propagationDistance);
            obj.setValueIfEmpty('innerScale', 0);
            obj.setValueIfEmpty('outerScale', Inf);
            obj.setValueIfEmpty('transverseSeparationInR0Units', 0);
            obj.setValueIfEmpty('hermiteGaussOrders', 0);
        end

        function setValueIfEmpty(obj, field, value)
            obj.(field) (isnan(obj.(field)) | isempty(obj.(field))) = value;
        end

        function computeDerivedQuantities(obj)
           wvl = obj.wavelength;
           L = obj.propagationDistance;
           delta1 = obj.gridSpacingSourcePlane;
           deltan = obj.gridSpacingObservationPlane;

           obj.numberOfPhasePlanes = obj.getMinimumNumberOfPhasePlanes;

           [obj.waistAtSourcePlane, obj.waistAtObservationPlane] = obj.getBeamWidthsAtSourceAndObservationPlanes;

           obj.waveNumber =  2*pi/wvl;

           obj.planePositions = linspace(0, L, obj.numberOfPhasePlanes);
           z = obj.planePositions;
           
           obj.complexBeamParameter = obj.computeComplexParameterAtSourcePlane;
           
           obj.gridSpacingVector = (1-z/L)*delta1 + z/L*deltan;
           
           obj.computeGammaStrength();
           obj.computeFriedCoherenceRadiusMatrix();
           obj.computeStructureConstantSquared();

           r0sw = obj.totalFriedCoherenceRadiusByStrength;
           r0sw = ~isinf(r0sw) .* r0sw;
           obj.regionOfInterestAtSourcePlane = 4*obj.waistAtSourcePlane + max(obj.transverseSeparationInR0Units) * max(r0sw);
           obj.regionOfInterestAtObservationPlane = 4*obj.waistAtObservationPlane;
        end

        function npl = getMinimumNumberOfPhasePlanes(obj)
            wvl = obj.wavelength;
            [Nx, Ny] = obj.getPhaseScreenGridSize;
            N = min(Nx, Ny);
            deltan = obj.gridSpacingObservationPlane;
            delta1 = obj.gridSpacingSourcePlane;

            zmax = min([delta1 deltan])^2 * N / wvl;
            npl = ceil(obj.propagationDistance/zmax)+1;
            npl = max(npl, 4);
        end

        function [w1, wn] = getBeamWidthsAtSourceAndObservationPlanes(obj)
            k = 2*pi/obj.wavelength;
            if obj.isFourthOrder
                k = 2*k;
            end
            z0 = obj.waistPosition;
            w0 = obj.waistRadius;
            
            L = obj.propagationDistance;
            w1 = w0 * sqrt(1 + (2*z0/(k*w0^2)).^2);
            wn = w0 * sqrt(1 + (2*(L-z0)/(k*w0^2)).^2);
        end

        function q = computeComplexParameterAtSourcePlane(obj)
            z0 = obj.waistPosition;
            w0 = obj.waistRadius;

            if obj.isFourthOrder
                k = 2*obj.waveNumber;
            end

            q = -z0 -1i*k*w0^2/2;
        end

        function computeGammaStrength(obj)
            obj.gammaStrength = obj.gammaRaw * 10^(-36/5);
        end

        function computeFriedCoherenceRadiusMatrix(obj)
            L = obj.propagationDistance;
            z = obj.planePositions;
            npl = obj.numberOfPhasePlanes;
            zmin = obj.turbulenceRegionStartPosition;
            zmax = obj.turbulenceRegionEndPosition;
            g = obj.gammaStrength';
            k = obj.waveNumber;
                        
            zmask = (z >= zmin & z <= zmax);
            ztmin_idx = find(zmask, 1, 'first');
            ztmax_idx = find(zmask, 1, 'last');
            s = sum((1-z(ztmin_idx:ztmax_idx)/L).^(5/3));
            
            r0 = (0.423 * k^2/(2.601 * L^(5/3) * s) * g.^(5/3)).^(-3/5);
            r0 = repmat(r0, 1, npl);
            r0(isinf(1./r0)) = 1e-5;
            r0(~repmat(zmask, length(g), 1)) = inf;
            r0(isnan(r0)) = inf;
            
            obj.friedCoherenceRadiusMatrix = r0;
            
            aux_mat = repmat( (z .* zmask)/L,length(g),1);
            obj.totalFriedCoherenceRadiusByStrength = ...
                sum(r0.^(-5/3) .* aux_mat.^(5/3), 2).^(-3/5);
        end
        function computeStructureConstantSquared(obj)
            k = obj.waveNumber;
            L = obj.propagationDistance;
            z = obj.planePositions;
            r0 = obj.totalFriedCoherenceRadiusByStrength;
            zmin = obj.turbulenceRegionStartPosition;
            zmax = obj.turbulenceRegionEndPosition;
            dz = abs(z(2)-z(1));

            zmask = (z >= zmin & z <= zmax);
            ztmin_idx = find(zmask, 1, 'first');
            ztmax_idx = find(zmask, 1, 'last');
            s = sum((z(ztmin_idx:ztmax_idx)/L).^(5/3));

            obj.structureConstantSquared = r0.^(-5/3)/(0.423 * k^2 * dz * s);
        end
        function plotConstraint1(obj, params)
            L = obj.propagationDistance;
            wvl = obj.wavelength;
            D1p = params.effectiveSourceROI;
            D2p = params.effectiveObsROI;
            
            d1 = linspace(0, 1.1*wvl*L/D2p, 100);
            dn = linspace(0, 1.1*wvl*L/D1p, 100);
            
            deltan_max = -D2p/D1p*d1 + wvl*L/D1p;
            plot(d1, deltan_max, 'k--', 'Linewidth', 2);
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
            
            deltaMax = abs(-D2p/D1p*delta1 + wvl*L/D1p);
            
            fail1 = (deltan >= deltaMax);
            fprintf('Constraint 1: ');
            if ~fail1
                fprintf('Satisfied\n');
            else
                fprintf('Not satisfied [deltan = %3.2e should be smaller than %3.2e]\n', ...
                    deltan, deltaMax);
            end
        end
        function fail2 = checkConstraint2(obj,params)
            wvl = obj.wavelength;
            L = obj.propagationDistance;
            [Nx, ~] = obj.getPhaseScreenGridSize;
            N = min(Nx);
            deltan = obj.gridSpacingObservationPlane;
            delta1 = obj.gridSpacingSourcePlane;
            
            D1p = params.effectiveSourceROI;
            D2p = params.effectiveObsROI;
            
            Nmin = (wvl*L + D1p*deltan + D2p*delta1)/ (2*delta1 * deltan);

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
        
    end
end

