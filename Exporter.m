classdef Exporter
    %Exporter Contains methods to export data to disk
    
    properties
    end
    
    methods(Static)
        function exportToDisk(results, simParams, inputFile)
            fileName = Exporter.getFullFilename(results, simParams.simulationType);
            fid = Exporter.getFileId(fileName);

            save(fileName, 'results');
            try
                for i = 1 : length(results)
                    Exporter.writeTitle(fid, Exporter.getDataSet(results,i));
                    Exporter.writeSimulationType(fid, simParams);
                    Exporter.saveHeader(fid, Exporter.getDataSet(results,i));
                    Exporter.saveDataArray(fid, Exporter.getDataSet(results,i));
                    Exporter.saveFigure(fileName, Exporter.getDataSet(results,i));
                end
                Exporter.saveTimeStamp(results, fid);
                Exporter.saveParametersFromFile(fid, inputFile, simParams);
                Exporter.saveSetPrivateProperties(fid, simParams);
            catch exception
                fclose(fid);
                rethrow(exception);
            end
            fclose(fid);
        end 
    end
    
    methods(Static, Access = private)
        function dataSet = getDataSet(data, idx)
            if ~iscell(data)
                dataSet = data;
                return;
            end
            dataSet = data{idx};
        end

        function writeTitle(fid, data)
            fprintf(fid,'%s\n', data.info.title);
        end

        function writeSimulationType(fid, simParams)
            if simParams.isFourthOrder
                isFourthOrderString = 'Fourth-order';
            else
                isFourthOrderString = 'Second-order';
            end
            if simParams.isInverted
                isInvertedString = 'with inversion';
            else
                isInvertedString = 'without inversion';
            end
            fprintf(fid,'Simulation type: %s, %s.\n', ...
                isFourthOrderString, isInvertedString);
        end

        function saveHeader(fid, data)
            fprintf(fid, 'Row: %s\n', data.info.labelRow); 
            fprintf(fid, 'Column: %s\n', data.info.labelColumn);
            fprintf(fid, 'Value: %s\n\n', data.info.labelZ);
        end
        
        function saveDataArray(fid, data)
            if iscell(data.values)
               return; 
            end
            fprintf(fid, 'r\\c');
            fprintf(fid, '\t%4.3g', data.columnParams);
            fprintf(fid, '\n');
            dataToWrite = [data.rowParams', data.values];
            for ln = 1 : size(dataToWrite, 1)
                fprintf(fid, '%3.3g\t', dataToWrite(ln,:));
                fprintf(fid, '\n');
            end
            fprintf(fid, '\n--------------------------------- \n');
        end
        function saveTimeStamp(data, fid)
            if iscell(data)
                data = data{end};
            end

            dateAndTime = regexp(data.info.date, '_', 'split');
            dt = dateAndTime{1};
            tm = dateAndTime{2};

            fprintf(fid,'\n\nSimulation ended at %s of %s.\n', tm, dt);
        end
        function saveParametersFromFile(fidWriteTo, inputFileName,simParams)
            fprintf(fidWriteTo, '\nParameters --\n');
            fidReadFrom = fopen(inputFileName);
            try
                currentLine = fgetl(fidReadFrom);
                while (currentLine ~= -1)
                    Exporter.writeToFileIfInPublicProperties(fidWriteTo,currentLine,simParams);
                    currentLine = fgetl(fidReadFrom);
                end
            catch exception
                fclose(fidReadFrom);
                rethrow(exception);
            end
            fclose(fidReadFrom);
        end
        function writeToFileIfInPublicProperties(fidWriteTo, currentLine, simParams)
            re = regexp(currentLine, '(\w)+\s*:','tokens');
            if isempty(re)
                fprintf(fidWriteTo, '%s\n', currentLine);
                return;
            end
            
            pblProp = Util.getSetPublicProperties(simParams);
            if ~ismember(re{1}{1}, pblProp) 
                fprintf(fidWriteTo, '%s\n', currentLine);
                return;
            end
            
            fprintf(fidWriteTo, '%s: ', re{1}{1});
            if (strcmpi(re{1}{1},'hermiteGaussOrders'))
                if length(simParams.(re{1}{1})) == 1
                    fprintf(fidWriteTo, '%02d', simParams.(re{1}{1}));
                else
                    fprintf(fidWriteTo, '%02d, ', simParams.(re{1}{1})(1:end-1));
                    fprintf(fidWriteTo, '%02d', simParams.(re{1}{1})(end));
                end
                fprintf(fidWriteTo, '\n');
                return;
            end

            if length(simParams.(re{1}{1})) == 1
                fprintf(fidWriteTo, '%4.3g', simParams.(re{1}{1}));
            else
                fprintf(fidWriteTo, '%4.3g, ', simParams.(re{1}{1})(1:end-1));
                fprintf(fidWriteTo, '%4.3g', simParams.(re{1}{1})(end));
            end
            fprintf(fidWriteTo, '\n');
        end
        function saveSetPrivateProperties(fid, simParams)
            pvtProp = Util.getSetPrivateProperties(simParams);
            if isempty(pvtProp)
                return;
            end
            fprintf(fid, '\n# Derived quantities\n');
            for p = 1 : length(pvtProp)
                 fprintf(fid, '%s: ', pvtProp{p});
                if length(simParams.(pvtProp{p})) == 1
                    fprintf(fid, '%4.4g', simParams.(pvtProp{p}));
                else
                    fprintf(fid, '%4.4g, ', simParams.(pvtProp{p})(1:end-1));
                    fprintf(fid, '%4.4g', simParams.(pvtProp{p})(end));
                end
                fprintf(fid, '\n');
            end
        end
        function saveFigure(fileName,data)
            if ~iscell(data.values)
                return;
            end
            saveas(gcf, [fileName, '.tiff']);
        end

        function fileName = getFullFilename(data, simulationType)
            [~, outFolder] = ConfigHelper.getIOFolders;
            dt =  datestr(date, 'yyyy-mm-dd');
            expFolder = fullfile(outFolder,dt);

            if ~exist(expFolder, 'dir')
                mkdir(expFolder);
            end

            tm = datestr(clock, 'hhMMss');
            fileName = fullfile(expFolder, [dt, ' ', simulationType, ' ', tm]);
        end

        function fid = getFileId(fileName)
            fid = fopen([fileName, '.dat'], 'w');
            if (fid == -1)
                error('exporter:fileOpen', 'Could not open %s.', fileName);
            end
        end
    end
end

