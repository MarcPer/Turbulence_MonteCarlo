classdef Exporter
    %Exporter Contains methods to export data to disk
    
    properties
    end
    
    methods(Static)
        function exportToDisk(ioPaths,data,userInput,simParams)
            fileName = Exporter.getFullFilename(ioPaths, 'dat');
            fid = Exporter.getFileId(fileName);
            try
                Exporter.writeTitle(fid,data);
                Exporter.writeSimulationType(fid,userInput);
                Exporter.saveHeader(fid,data);
                Exporter.saveDataArray(fid,data);
                Exporter.saveTimeStamp(fid);
                Exporter.saveParametersFromFile(fid,ioPaths.inputParametersFileName,simParams);
                Exporter.saveSetPrivateProperties(fid, simParams);
                Exporter.saveFigure(ioPaths,data);
            catch exception
                fclose(fid);
                rethrow(exception);
            end
            fclose(fid);
        end 
    end
    
    methods(Static, Access = private)
        function writeTitle(fid, data)
            fprintf(fid,'%s\n', data.info.title);
        end
        function writeSimulationType(fid, userInput)
            if userInput.isFourthOrder
                isFourthOrderString = 'Fourth-order';
            else
                isFourthOrderString = 'Second-order';
            end
            if userInput.isInverted
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
        end
        function saveTimeStamp(fid)
            dt = datestr(date, 'yyyy-mm-dd');
            tm = datestr(clock, 'HH:MM:SS');
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
                    fprintf(fid, '%4.3g', simParams.(pvtProp{p}));
                else
                    fprintf(fid, '%4.3g, ', simParams.(pvtProp{p})(1:end-1));
                    fprintf(fid, '%4.3g', simParams.(pvtProp{p})(end));
                end
                fprintf(fid, '\n');
            end
        end
        function saveFigure(ioPaths,data)
            if ~iscell(data.values)
                return;
            end
            filename = Exporter.getFullFilename(ioPaths,'tiff');
            saveas(gcf,filename);
        end
        function fileName = getFullFilename(ioPaths, extension)
            filePath = ioPaths.getExportPath;
            if ~exist(filePath, 'dir')
                mkdir(filePath)
            end           
            fileName = ioPaths.getExportFileName(extension);
            fileName = fullfile(filePath, fileName);
        end
        function fid = getFileId(fileName)
            fid = fopen(fileName, 'w');
            if (fid == -1)
                error('exporter:fileOpen', 'Could not open %s.', fileName);
            end
        end
    end
end

