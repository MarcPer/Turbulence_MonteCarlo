classdef Plotter<handle
    %Plotter Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods (Static)
        function plotIntensityProfilesForEachGamma(data)
            figure;
            set(gcf, 'Units', 'normalized', 'Position', [0.625 0.28 0.365 0.61]);
            row = data.rowParams;       % Separation
            col = data.columnParams;    % Gamma
            
            for iRow = 1 : length(row)
                for iCol = 1 : length(col)
                    i = 3*(iRow-1)+iCol;
                    subplot( length(row), length(col), i);
                    imagesc(data.values{iCol}(:, :, iRow));
                    axis off;
                    str = sprintf('gamma = %2.2g um, sep = %2.2g r0', col(iCol), row(iRow));
                    title(str);
                end
            end
        end
        function plot2D(data)
            if isstruct(data)
                x = data.columnParams;
                y = data.values;
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

