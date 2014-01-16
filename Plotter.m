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
            Plotter.drawPlotInformation(data);
        end
        function drawPlotInformation(data)
            if ~isstruct(data)
                return
            end
            
            if isfield(data.info,'title')
                str = data.info.title;
                title(str);
            end
            if isfield(data.info,'labelColumn')
                str = data.info.labelColumn;
                xlabel(str);
            end
            if isfield(data.info,'labelZ')
                str = data.info.labelZ;
                ylabel(str);
            end
            if isfield(data.info,'labelLegend')
                legend(data.info.labelLegend)
            end
        end
    end
end

