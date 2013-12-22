classdef Exporter
    %Exporter Contains methods to export data to disk
    
    properties
    end
    
    methods(Static)
        function exportToDisk(ioPaths,data)
            filePath = ioPaths.getExportPath;
            if ~exist(filePath, 'dir')
                mkdir(filePath)
            end           
            fileName = ioPaths.getExportFileName;
            fileName = fullfile(filePath, fileName);
            fid = fopen(fileName, 'w');
            if (fid == -1)
                error('exporter:fileOpen', 'Could not open %s.', fileName);
            end
            try
                Exporter.saveHeader(fid,data);
                Exporter.saveData(fid,data);
                Exporter.saveTimeStamp(fid);
                Exporter.saveParameters(fid,ioPaths.inputParametersFileName);
            catch exception
                fclose(fid);
                rethrow(exception);
            end
            fclose(fid);
        end 
    end
    
    methods(Static, Access = private)
        function saveHeader(fid, data)
            fprintf(fid, 'Row: %s\n', data.info.labelRow); 
            fprintf(fid, 'Column: %s\n', data.info.labelColumn);
            fprintf(fid, 'Value: %s\n\n', data.info.labelZ);
        end
        function saveData(fid, data)
            fprintf(fid, 'r\\c');
            fprintf(fid, '\t%3.3g', data.data.columnParams);
            fprintf(fid, '\n');
            dataToWrite = [data.data.rowParams, data.data.values];
            for ln = 1 : size(dataToWrite, 1)
                fprintf(fid, '%3.3g\t', dataToWrite(ln,:));
                fprintf('\n');
            end
        end
        function saveTimeStamp(fid)
            dt = datestr(date, 'yyyy-mm-dd');
            tm = datestr(clock, 'HH:MM:SS');
            fprintf(fid,'\n\nSimulação finalizada às %s do dia %s.\n', tm, dt);
        end
        function saveParameters(fidWriteTo, inputFileName)
            fprintf(fidWriteTo, '\nParâmetros --\n');
            fidReadFrom = fopen(inputFileName);
            try
                currentLine = fgetl(fidReadFrom);
                while (currentLine ~= -1)
                    fprintf(fidWriteTo, '%s\n', currentLine);
                    currentLine = fgetl(fidReadFrom);
                end
            catch exception
                fclose(fidReadFrom);
                rethrow(exception);
            end
            fclose(fidReadFrom);
        end
    end
end

