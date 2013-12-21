classdef Plotter<handle
    %Plotter Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods (Static)
        function plotIntensityProfilesForEachGamma(gammaValues, intProfiles)
            figure;
            set(gcf, 'Units', 'normalized', 'Position', [0.625 0.28 0.365 0.61]);
            lenG = length(gammaValues);
            [numRows, numCols] = Util.findOptimumSubplotGrid(lenG);
            
            for iG = 1 : lenG
                subplot( numRows, numCols, iG);
                imagesc(intProfiles{iG}(:, :, 1));
                axis off;
                title(['\gamma = ', num2str(gammaValues(iG)*1e6), ' \mum']);
            end
        end
        
    end
end

