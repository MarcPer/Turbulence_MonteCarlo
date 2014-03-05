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
        function plot2D(data, varargin)
            if isstruct(data)
                x = data.columnParams;
                y = data.values;
            else
                y = data;
                x = 1 : length(data);
            end
            if length(y) < 2
                return;
            end
            if size(varargin)==0
                figure;
                plot(x,y, 'LineWidth', 2);
                grid on;
                Plotter.drawPlotInformation(data);
                return;
            end
            if isfloat(varargin{1})
                hold all;
                plot(varargin{1}, x,y, 'LineWidth', 2);
                hold off;
            end
        end

        function plotSemiLogX(data)
            if isstruct(data)
                x = data.columnParams;
                y = data.values;
            else
                y = data;
                x = 1 : length(data);
            end
            if length(y) < 2
                return;
            end
            figure;
            semilogx(x,y, 'LineWidth', 2);
            grid on;
            Plotter.drawPlotInformation(data);
        end

        function plotBars(data)
            figure;
            nGamma = numel(data);
            [nrow, ncol] = Util.findOptimumSubplotGrid(nGamma);

            for iGamma = 1 : nGamma

                subplot(nrow, ncol, iGamma);
                bar3(data{iGamma}.values);
                Plotter.drawBarPlotInformation(data{iGamma});
            end
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

        function drawBarPlotInformation(data)
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
            if isfield(data.info,'labelRow')
                str = data.info.labelRow;
                ylabel(str);
            end
            if isfield(data.info, 'tickX')
                tck = data.info.tickX;
                set(gca, 'XTickLabel', tck);
            end
            isfield(data.info, 'tickY')
                tck = data.info.tickY;
                set(gca, 'YTickLabel', tck);
            end

    end
end

