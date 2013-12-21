classdef Exporter
    %Exporter Contains methods to export data to disk
    
    properties
    end
    
    methods(Static)
        function exportToDisk(ioPaths,data)
            filePath = ioPaths.getExportPath;
            fileName = ioPaths.getExportFileName;
            fileName = fullfile(filePath, fileName);
            fid = fopen(fileName, 'w');
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
        function saveHeader(fid, data)
            fprintf(fid, 'Row: %s\n', data.info.labelRow); 
            fprintf(fid, 'Column: %s\n', data.info.labelColumn);
            fprintf(fid, 'Value: %s\n\n', data.info.labelZ);
        end
        function saveData(fid, data)
            fprintf(fid, '\t%s', data.data.columnParams);
            fprintf(fid, '%3.3g\n', [data.data.rowParams; data.data.values']);
        end
        function saveParameters(fid, inputFileName)
            fprintf(fid,  
        end
    end
    
end

