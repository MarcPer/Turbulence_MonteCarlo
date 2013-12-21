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
        function plot2D(data)
            if isstruct(data)
                x = data.data.columnParams;
                y = data.data.values;
            else
                y = data;
                x = 1 : length(data);
            end
            figure;
            plot(x,y, 'LineWidth', 2);
            grid on;
            if isstruct(data)
                if isfield(data.info,'title')
                    title(data.info.title);
                end
                if isfield(data.info,'labelColumn')
                    xlabel(data.info.labelColumn);
                end
                if isfield(data.info,'labelZ')
                    ylabel(data.info.labelZ);
                end
            end
        end
    end
end

