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

ioPaths = IOPaths;   % Imports functions and defines methods for export paths
usrIn = UserInput;

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
    apertureRadius = turbSimulator.simulationParameters.gridSpacingObservationPlane;
    pwrAndSI = turbSimulator.getPowerAndSIOnCircularAperture(apertureRadius,'Normalized', true);
    % pwrGamma = turbSimulator.getIrradianceForEachGamma('Normalized', true);
    % modeMatching = turbSimulator.getModeMatching;
%     parity = turbSimulator.getModeParity;
catch exception
    % SHUTDOWN COMPUTER?
    usrIn.shutdownComputer;
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
%     Plotter.plotIntensityProfilesForEachGamma(pwrGamma);
    Plotter.plot2D(pwrAndSI{1});
    Plotter.plot2D(pwrAndSI{2});
%     Plotter.plotBars(parity);
%     Plotter.plot2D(Calculator.computeErrorRateVsRelativeLengths(parity));
        
    % EXPORT RESULTS (only if simulation was completed)
    Exporter.exportToDisk(ioPaths, pwrAndSI, usrIn, turbSimulator.simulationParameters, 'pwrAndSI');
catch exception
    % SHUTDOWN COMPUTER?
    usrIn.shutdownComputer;
    rethrow(exception);
end
usrIn.shutdownComputer;