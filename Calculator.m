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
        function [pwrSlit, plotInfo] = computePowerThroughSlit(intProfile, simParams)
            slitWidthPx = round(simParams.slitWidth / ...
                simParams.gridSpacingObservationPlane);
            
            if simParams.isFourthOrder
                pwrSlit = coincidenceSlitIntegratePx(intProfile, slitWidthPx);
            else
                pwrSlit = intensitySlitIntegratePx(intProfile, slitWidthPx);
            end
            
            tit = 'Intensity vs Turbulence Strength';
            labelX = '\gamma';
            labelY = 'intensity';
            plotInfo = struct('title', tit, ...
                'labelX', labelX, 'labelY', labelY);
            
        end
    end
    
    methods(Access = private)
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
            hght = zeros(numCells, 1);
            wdth = zeros(numCells, 1);
            numSeps = zeros(numCells, 1);
            
            maxColumn = zeros(numCells, numSeps);
            coincSlitPwr = zeros(numCells, numSeps);
            
            for iCell = 1 : numCells
                [hght(iCell), wdth(iCell), numSeps(iCell)] = ...
                    size(intProfile{iCell});
                [~, maxColumn(iCell)] = max( max(intProfile{iCell}, [], 1) );
                
                cSlit = triang(2*slitWidth -1)';
                lenSlit = length(cSlit);
                
                % Adjust the slit function to be of the same size as matrix 'A'
                if (lenSlit > wdth)
                    lenDiff = lenSlit - wdth;
                    cSlit([1:floor(lenDiff/2), lenSlit-floor(lenDiff/2)+1:end]) = [];
                else
                    lenDiff = wdth - lenSlit;
                    cSlit = [zeros(1, floor(lenDiff/2)), cSlit, ...
                        zeros(1, ceil(lenDiff/2))];
                end
                
                % Center slit function horizontally on maximum element of 'A'
                [~, xSlit] = max(cSlit);
                cSlit = Displace(cSlit, maxColumn(iCell) - xSlit, 0);
                
                cSlit = repmat(cSlit, hght, 1);
                coincSlitPwr(iCell) = sum( sum( ...
                    cSlit .* intProfile{iCell} ));
                
            end
            
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
            hght = zeros(numCells, 1);
            wdth = zeros(numCells, 1);
            maxColumn = zeros(numCells, 1);
            singSlitPwr = zeros(numCells, 1);
            
            for iCell = 1 : numCells
                [hght(iCell), wdth(iCell)] = size(intProfile{iCell});
                [~, maxColumn(iCell)] = max( max(intProfile{iCell}, 1) );
                                
                % Adjust the slit function to be of the same size as matrix 'A'
                if (slitWidth >= wdth)
                    slit = ones(1,wdth);
                else
                    slit = [zeros(1, maxColumn(iCell)-floor(slitWidth/2)), ...
                        ones(1, slitWidth), ...
                        zeros(1, wdth-maxColumn(iCell)-ceil(slitWidth/2))];
                end
                
                % Center slit function horizontally on maximum element of 'A'
                [~, xSlit] = max(slit);
                slit = Displace(slit, maxColumn(iCell)-xSlit, 0);
                
                slit = repmat(slit, hght, 1);
                singSlitPwr(iCell) = sum( sum( ...
                    slit .* intProfile{iCell} ));
                
            end
            
        end 
    end
end

