classdef Plotter<handle
    %Plotter Reads data structure, plots results and logs meta-information
    %   The Plotter class creates a log so that information about a
    %   particular plot can be found. This information associates the data
    %   shown with its source, so the used parameters can be retrieved.
    
    properties
        plotInfo;
    end
    
    methods (Access = public)
        function plotter = Plotter()
            plotter.plotInfo = struct();
            plotter.plotInfo.dataDates = {};
            plotter.plotInfo.dataTitles = {};
        end

        function plotIntensityProfilesForEachGamma(obj, data)
            figure;
            set(gcf, 'Units', 'normalized', 'Position', [0.625 0.28 0.365 0.61]);
            row = data.rowParams;       % Separation
            col = data.columnParams;    % Gamma
            
            obj.setPlotInfo(data);
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

        function plot2D(obj, data, varargin)
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

            obj.setPlotInfo(data);
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

        function plotSemiLogX(obj, data)
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

            obj.setPlotInfo(data);
            figure;
            semilogx(x,y, 'LineWidth', 2);
            grid minor;
            Plotter.drawPlotInformation(data);
        end

        function plotSemiLogY(obj, data, varargin)
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

            obj.setPlotInfo(data);
            if size(varargin)==0
                figure;
                semilogy(x,y, 'LineWidth', 2);
                grid minor;
                Plotter.drawPlotInformation(data);
                return;
            end
            if isfloat(varargin{1})
                hold all;
                semilogy(varargin{1}, x,y, 'LineWidth', 2);
                hold off;
            end
        end

        function plotLogLog(obj, data, varargin)
            if isstruct(data)
                x = data.columnParams;
                y = data.values';
            else
                y = data;
                x = 1 : length(data);
            end

            if length(y) < 2
                return;
            end

            obj.setPlotInfo(data);
            if size(varargin)==0
                figure;
                h = loglog(x,y, 's-', 'LineWidth', 2);
                disp(get(h, 'Color'));
                for i = 1 : numel(h)
                    set(h(i), 'MarkerFaceColor', get(h(i), 'Color'));
                end
                grid on;
                grid minor;
                Plotter.drawPlotInformation(data);
                return;
            end
            if isfloat(varargin{1})
                hold all;
                loglog(varargin{1}, x,y, 'LineWidth', 2);
                hold off;
            end
        end

        function plotBars(obj, data)
            figure;
            nGamma = numel(data);
            [nrow, ncol] = Util.findOptimumSubplotGrid(nGamma);

            obj.setPlotInfo(data);
            for iGamma = 1 : nGamma
                subplot(nrow, ncol, iGamma);
                bar3(data{iGamma}.values);
                Plotter.drawBarPlotInformation(data{iGamma});
            end
        end

        function setPlotInfo(obj, data)
            if iscell(data)
                data = data{end};
            end

            if ~isstruct(data)
                return
            end

            if isfield(data.info, 'date')
                obj.plotInfo.dataDates = [obj.plotInfo.dataDates; {data.info.date}];
            end

            if isfield(data.info, 'title')
                obj.plotInfo.dataTitles = [obj.plotInfo.dataTitles; {data.info.title}];
            end
        end

        function clearPlotInfo(obj)
            obj.plotInfo = struct();
        end

        function exportPlot(obj)
            [fileName, pathName] = uiputfile;
            if ~fileName
                return;
            end

            saveName = fullfile(pathName, fileName);
            logFile = fullfile(pathName, 'plotLog.txt');
            fid = fopen(logFile, 'a');

            try
                saveas(gcf, saveName);
                obj.logPlotInfo(fid, fileName);
            catch exception
                fclose(fid);
                rethrow(exception);
            end
            
            fclose(fid);
        end
    end

    methods (Access = private)
        function logPlotInfo(obj, fid, fileName)
            dt = datestr(date, 'yyyy-mm-dd');
            tm = datestr(clock, 'HHMMSS');
            fprintf(fid, 'Plot exported to %s -- (%s, %s)\n', fileName, dt, tm);

            for i = 1 : numel(obj.plotInfo.dataDates)
                fprintf(fid, '\t%s: %s\n', obj.plotInfo.dataDates{i}, obj.plotInfo.dataTitles{i});
            end

            fprintf(fid, '----------------------\n\n');
        end

    end

    methods (Static)
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
            if isfield(data.info, 'tickY')
                tck = data.info.tickY;
                set(gca, 'YTickLabel', tck);
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


    end
end

