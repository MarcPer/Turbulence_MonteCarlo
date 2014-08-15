clear; close all;
set(0,'ShowHiddenHandles','on');
delete(get(0,'Children'));
addpath('jsonlab');

inputFile = ConfigHelper.getInputParametersFile;
ConfigHelper.setOutputFolder;
simParams = SimulationParameters(loadjson(inputFile));
isAbort = simParams.checkConstraints;
if isAbort
	return;
end

turbSimulator = TurbulenceSimulator(simParams);
results = turbSimulator.simulate;