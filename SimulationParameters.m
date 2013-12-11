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
           simParams.checkIfFourthOrder(varargin);
           fData = simParams.openParametersFile();
           simParams.readParameters(fData);
           simParams.setDefaultValueForBlankParameters();
           simParams.computeDerivedQuantities();
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
    end
        
    methods(Static)
        function fData = openParametersFile()
            fid = fopen('inputParameters.dat');
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
        function sg = superGaussianFilter(obj)
            [x1, y1] = obj.getMeshGridAtSourcePlane();
            N = obj.transverseGridSize;
            sg = exp(-(x1/(0.47*N*delta1)).^16 ... 
                -(y1/(0.47*N*delta1)).^16);
        end
        function Uin = getInputField(obj)
            [x1,y1] = obj.getMeshGridAtSourcePlane();
            qz = obj.complexBeamParameter;
            k = obj.waveNumber;
            Uin = 1/qz*exp(1i*k/(2*qz)*(x1.^2+y1.^2)); 
        end
        function [x1, y1] = getMeshGridAtSourcePlane(obj)
            N = obj.transverseGridSize;
            delta1 = obj.gridSpacingSourcePlane;
            [x1,y1] = meshgrid((-N/2 : N/2-1) * delta1);
        end
    end
end

