%% Clear and close all
close all; clear all;
%
% This script simulates the effect of turbulence by modelling it as a
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
constraintReturnCode = simParams.constraintAnalysis;

isAbort = usrIn.abortWhenConstraintFail(constraintReturnCode, usrIn.enumShutdown);
if isAbort
    return
end
close all; clear isAbort;

turbSimulator = TurbulenceSimulator(simParams);

%% RUN SIMULATION
intProfileGamma = turbSimulator.getIntensityForEachGamma('Normalized', true);

if turbSimulator.isAborted
    fprintf('Simulation aborted.');
    return;
end

%% COMPUTE RELEVANT VALUES
slitWidth = 50e-6;
pwrSlit = Calculator.computePowerThroughSlit(intProfileGamma, simParams);

%% PLOT RESULTS

% Beam transverse profiles at observation plane ('measured gammas only')
close all;
Plotter.plotIntensityProfilesForEachGamma(simParams.gammaStrength, intProfileGamma);
Plotter.plot2D(pwrSlit);

% %% EXPORT RESULTS (only if simulation was completed)
Exporter.exportToDisk(ioPaths, pwrSlit, usrIn);
 
% %% SHUTDOWN COMPUTER?
% wtime = 90;         % Time limit to abort shutdown
% 
% if usrIn.enumShutdown
%     if usrIn.enumShutdown == 1
%         shut_mode = 'shutdown';
%         shut_cmd = 'shutdown /s';
%     else
%         shut_mode = 'hibernation';
%         shut_cmd = 'shutdown /h';
%     end
%     
%     fprintf('Computer %s enabled.\n', shut_mode);
%     tic;
%     reply_shut = 'k';
%     while (~sum(strcmpi(reply_shut, {'y', 'n', ''})) && toc < wtime)
%         if usrIn.enumShutdown == 1
%         reply_shut = waitinput('Proceed with shutdown? Y/N [Y]:    ', ...
%             wtime, 's');
%         else
%         reply_shut = waitinput('Proceed with hibernation? Y/N [Y]:    ', ...
%             wtime, 's');
%         end
%         fprintf('\n');
%     end
%     if (isempty(reply_shut) || ~sum(strcmpi(reply_shut, {'y', 'n'})))
%         reply_shut = 'Y';
%     end
%     if strcmpi(reply_shut,'y')
%         if usrIn.enumShutdown == 1
%             fprintf('\nShutting down.\n');
%         else
%             fprintf('\nHibernating.\n');
%         end
%         pause(3);
%         usrIn.enumShutdown = 0;
%         system(shut_cmd);
%     end
%     if usrIn.enumShutdown == 1
%         fprintf('Shutdown aborted.\n');
%     elseif usrIn.enumShutdown == 2
%         fprintf('Hibernation aborted.\n');
%     end
% end
