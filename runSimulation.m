clear; close all;
set(0,'ShowHiddenHandles','on');
delete(get(0,'Children'));
addpath('jsonlab');

inputFile = ConfigHelper.getInputParametersFile;
ConfigHelper.setOutputFolder;
simParams = SimulationParameters(loadjson(inputFile));
shutdownOrHibernate = UserInput.confirmShutdownOrHibernate(simParams.shutdownOrHibernate);

isAbort = simParams.checkConstraints;
if isAbort
	return;
end

turbSimulator = TurbulenceSimulator(simParams);
try
	results = turbSimulator.simulate;
catch exception
    UserInput.shutdownComputer(shutdownOrHibernate);
    set(0,'ShowHiddenHandles','on')
    delete(get(0,'Children'))
    rethrow(exception);
end