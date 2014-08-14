classdef ConfigHelper
	methods(Static)
		function inFile = getInputParametersFile()
			inFolder = ConfigHelper.getInputFolder();
			[inFile, inFolder] = uigetfile(fullfile(inFolder,'*.json'), 'Choose input parameters file');
			if (inFile == 0)
				error('turbSimulator:noInputFile', 'No input file selected')
			end
			inFile = fullfile(inFolder, inFile);
			ConfigHelper.writeInputParametersFolder(inFolder);
		end

		function inFolder = getInputFolder()
			inFolder = '';
			fileID = fopen('turbSimulator.conf', 'r');

			% Configuration file does not exist
			if (fileID == -1)
				return;
			end

			% If it does exist
			try
				confParams = textscan(fileID, '%s', 'Delimiter', '\n');
			catch exception
				fclose(fileID);
			end
			fclose(fileID);

			% but folder information is not there
			if (isempty(confParams) || length(confParams{1})<2 )
				return;
			end

			% Otherwise get the information
			inFolder = confParams{1}{2};
		end

		function writeInputParametersFolder(inFolder)
			fileID = fopen('turbSimulator.conf', 'w');
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
