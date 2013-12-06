classdef SimulationParameters<handle
    %   Distance unit: m
    
    properties(GetAccess = public, SetAccess = public)
        % Geometry
        propagationDistance = 3000;
        numberOfPhasePlanes = 20;
        % Turbulence
        gammaStrength = [0; 4.29; 7.99; 12.1; 16.6; 21.3; 27.9] *1e-6;
        outerScale = Inf;
        innerScale = 0;
        % Beam
        wavelength = 325e-9;
        beamWidthObsPlane = 60e-6;
        % Simulation
        numberOfRealizations = 50;
        transverseGridSize = 512;
        gridSpacingSourcePlane = 5.5e-5;
        gridSpacingObservationPlane = 7e-6;
    end
    
       
    methods(Access = public)
        function MM = MatrizMueller(fileID)
            MM.fileID = fileID;
            MM.getHeader();
            MM.data = struct;
            MM.getTableHeader();
            MM.retrieveData();
        end
        
        function createMatrices(MM, numMat)
            if nargin == 1
                numMat = 1;
            end
            
            numAngles = length(MM.data.angle);
            
            MM.matrix = zeros(4,4,numMat,numAngles);
            m = zeros(5,1);
            tableHeader = fieldnames(MM.data);
            
            for ang = 1 : numAngles
                if ismember('error', tableHeader)
                    errorVect = [MM.data.(tableHeader{5})(ang); ...
                        MM.data.(tableHeader{7})(ang); ...
                        MM.data.(tableHeader{9})(ang) ; ...
                        MM.data.(tableHeader{11})(ang); ...
                        MM.data.(tableHeader{13})(ang)];
                else
                    errorVect = zeros(5,1);
                end

                for i = 1 : numMat
                    m(1) = -MM.data.(tableHeader{4})(ang);
                    m(2) = MM.data.(tableHeader{6})(ang);
                    m(3) = MM.data.(tableHeader{8})(ang);
                    m(4) = MM.data.(tableHeader{10})(ang);
                    m(5) = MM.data.(tableHeader{12})(ang);

                    m = m + errorVect .* randn(5,1);

                    MM.matrix(:,:,i,ang) = ...
                        [1 m(1) 0 0; m(1) m(2) 0 0; ...
                        0 0 m(3) m(4); 0 0 -m(4) m(5)];
                end
            end
        end
        
        function numUnphysical = filterUnphysical(MM)
            % [~, ~, numRand, ~] = size(MM.matrix);
            numMat = length(MM.matrix(1,1,:));
            MM.eigenValues = zeros(4, numMat);
            negPos = [];
            
            for i = 1 : numMat
                hmat = Hmatrix(MM.matrix(:,:,i));
                eigenVal = real(eig(hmat));
                if (any(eigenVal < 0))
                    negPos = [negPos; i];
                end
                eigenVal = eigenVal / sum(eigenVal);
                MM.eigenValues(:,i) = eigenVal;
            end
            
            MM.matrix(:,:,negPos) = [];
            MM.eigenValues(:,negPos) = [];
            MM.data.angle(negPos) = [];
            numUnphysical = length(negPos);
            
        end
        
        function generateDmEm(MM)
            numMat = length(MM.matrix(1,1,:));
            MM.DmEm = zeros(numMat, 2);
            
            for i = 1 : numMat
                [MM.DmEm(i,1), MM.DmEm(i,2)] = calcDmEm(MM.eigenValues(:,i));
            end
            
        end
        
    end
      
        
        
     methods(Access = private)   
        function getHeader(MM)
            header = fgetl(MM.fileID);
            re = regexp(header, '([ \w]+) - ([\d\.]+)','tokens');
            if (isempty(re))
                MM.material = 'NotFound';
                MM.wavelength = 'NotFound';
            else
                MM.material = re{1}{1};
                MM.wavelength = str2double(re{1}{2});
            end
            
        end
        
        function getTableHeader(MM)
            % Scan for table header
            [isAngle, isTheta] = deal(0);
            while(~isAngle && ~isTheta)
                currentLine = fgetl(MM.fileID);
                isAngle = ~isempty(regexpi(currentLine, 'angle.+error', 'ONCE'));
                isTheta = ~isempty(regexpi(currentLine, 'theta.+error', 'ONCE'));
            end
            
            % Use column names for keys in the MM object structure
            reHeader = regexp(currentLine, '[\w\d-()\/]+', 'match');
            keyCount = 1;
            
            for i = 1 : length(reHeader)
                reHeader{i} = regexprep(reHeader{i}, '-', 'minus');
                reHeader{i} = regexprep(reHeader{i}, '\/', 'div');
                reHeader{i} = regexprep(reHeader{i}, '[()]', '_');
                reHeader{i} = regexprep(reHeader{i}, 'theta', 'angle');
                if (i>1)
                    fieldNames = fieldnames(MM.data);
                    if ismember(reHeader{i},fieldNames)
                        keyCount = keyCount + 1;
                        reHeader{i} = [reHeader{i}, num2str(keyCount)];
                    end
                end
                MM.data.(reHeader{i}) = [];
            end
        end
        
        function retrieveData(MM)
            timeOut = 0;
            reHeader = fieldnames(MM.data);
            pattern = '-?[\d\.]+';
            
            while (timeOut < 1000)
                timeOut = timeOut + 1;
                string = fgetl(MM.fileID);
                if (isfloat(string))
                    break
                end
                if (~isempty(regexp(string,'[a-zA-Z_]','ONCE')))
                    break
                end
                if (ischar(string))
                    re = regexprep(string, '\-+\s', '0\t');
                    re = regexprep(re, '\-+$', '0');
                    re = regexp(re,pattern,'match');
                    if (~isempty(re))
                        for i = 1 : length(re)
                            MM.data.(reHeader{i}) = ...
                                [MM.data.(reHeader{i}); str2double(re{i})];
                        end
                    end
                end
            end
        end
        
        
    end
    
end

