%% Clear and close all
close all; clear all;
%
% This script simulates the effect of turbulence, modelling it as a
% series of phase screens. Technical details are taken from
%
% Numerical Simulation of Optical Wave Propagation
%            (with examples in MATLAB)
%               Jason D. Schmidt
%
% <LENGTH UNITS: m>

ioPaths = IOPaths;  % Imports functions and defines methods for export paths
usrIn = UserInput;  % Stores user input information
pltr = Plotter(ioPaths);    % Stores data source information for later reference

% Ask user to choose simulation type
usrIn.getSimulationType;

% Shutdown computer at the end of the script?
%   (0 = NO, 1 = SHUTDOWN, 2 = HIBERNATE)
usrIn.enumShutdown = 0;

%% SETUP AND SIMULATION PARAMETERS
% Setup geometry
simParams = SimulationParameters(ioPaths,'FourthOrder', usrIn.isFourthOrder, 'Inverted', usrIn.isInverted);
simParams.setPointDetectorAtObservationPlane(false);
constraintReturnCode = simParams.constraintAnalysis;

isAbort = usrIn.abortWhenConstraintFail(constraintReturnCode, usrIn.enumShutdown);
if isAbort
    return
end
close all; clear isAbort;

turbSimulator = TurbulenceSimulator(simParams);
clear simParams;

%% RUN SIMULATION
try
    % apertureRadius = turbSimulator.simulationParameters.gridSpacingObservationPlane;
    % pwrAndSI = turbSimulator.getPowerAndSIOnCircularAperture(apertureRadius,'Normalized', true);
    pwrGamma = turbSimulator.getIrradianceForEachGamma('Normalized', true);
    slitPwr = Calculator.computePowerThroughSlit(pwrGamma.values, turbSimulator.simulationParameters);
    % modeMatching = turbSimulator.getModeMatching;
    % parity = turbSimulator.getModeParity;
catch exception
    % SHUTDOWN COMPUTER?
    usrIn.shutdownComputer;
    set(0,'ShowHiddenHandles','on')
    delete(get(0,'Children'))
    rethrow(exception);
end

if turbSimulator.isAborted
    fprintf('Simulation aborted.\n');
    return;
end

%% PLOT RESULTS

% Beam transverse profiles at observation plane ('measured gammas only')
close all;

try
    % pltr.plotIntensityProfilesForEachGamma(pwrGamma);
    % pltr.plot2D(pwrAndSI{1});
    % pltr.plot2D(pwrAndSI{2});
    % pltr.plotBars(parity);
    % pltr.plot2D(Calculator.computeErrorRateVsRelativeLengths(parity));
    pltr.plot2D(slitPwr);
        
    % EXPORT RESULTS (only if simulation was completed)
    Exporter.exportToDisk(ioPaths, slitPwr, usrIn, turbSimulator.simulationParameters, 'slitPwr');
catch exception
    % SHUTDOWN COMPUTER?
    usrIn.shutdownComputer;
    rethrow(exception);
end
usrIn.shutdownComputer;