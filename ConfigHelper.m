classdef ConfigHelper
	methods(Static)
		function inFolder = getInputParametersFolder()
			fileID = fopen('turbSimulator.conf', 'r');
			if (fileID == -1)
				inFolder = ConfigHelper.writeInputParametersFolder();
				return;
			end

			try
				confParams = textscan(fileID, '%s', 'Delimiter', '\n');
			catch exception
				fclose(fileID);
			end
			fclose(fileID);

			if (isempty(confParams) || length(confParams{1})<2 )
				inFolder = ConfigHelper.writeInputParametersFolder();
				return;
			end

			inFolder = confParams{1}{2};
		end

		function inFolder = writeInputParametersFolder()
			inFolder = uigetdir();
			fileID = fopen(fullfile(inFolder, 'turbSimulator.conf'), 'w');
			try
				fprintf(fileID, 'Input parameters directory:\n%s', inFolder);
			catch exception
				fclose(fileID);
				rethrow(exception);
			end
			fclose(fileID);
		end

	end

end
