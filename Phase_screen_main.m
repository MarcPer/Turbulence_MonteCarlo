%% Clear and close all
close all; clear all;

%% INITIALIZATION CELL
% This script simulates the effect of turbulence by modelling it as a
% series of phase screens. Technical details are taken from
%
% Numerical Simulation of Optical Wave Propagation
%            (with examples in MATLAB)
%               Jason D. Schmidt
%
% <LENGTH UNITS: m>

% Set running folder
user_folder = userpath;
user_folder(end) = [];
codeFolder = cd(user_folder);
cd ../..
user_folder = fullfile(pwd, 'Dropbox');
cd(codeFolder);

% Import additional functions
fun_folder = fullfile(user_folder, 'MATLAB', ...
    'Numerical simulation of optical wave propagation');
addpath(fun_folder);
fun_folder = fullfile(user_folder, 'MATLAB', 'LuckyImaging');
addpath(fun_folder);
fun_folder = fullfile(user_folder, 'MATLAB', 'Downloaded', 'dlmcell');
addpath(fun_folder);
fun_folder = fullfile(user_folder, 'MATLAB', 'Downloaded', 'waitinput');
addpath(fun_folder);
fun_folder = fullfile(user_folder, 'MATLAB');
addpath(fun_folder);
clear fun_folder;

% Ask user to choose simulation type
Phase_screen_choose_simulation;

% Shutdown computer at the end of the script?
%   (0 = NO, 1 = SHUTDOWN, 2 = HIBERNATE)
shut = 2;
Phase_screen_shutdown_input;

%% SETUP AND SIMULATION PARAMETERS
% Setup geometry
Phase_screen_params;
params_flag = 1;

% Check sampling constraints
Phase_screen_constraint_analysis

abort = 0;
reply = 'k';
if constr_err
    while ~sum(strcmpi(reply, {'y', 'n', ''})) 
        reply = input('Abort simulation? Y/N [Y]: ', 's');
    end
    if isempty(reply)
        reply = 'Y';
    end
    if sum(strcmpi(reply,{'y', 'Y'}))
        abort = 1;
    end
end

if ~abort
    fprintf('Press any key to continue...\n');
    pause;
end

if abort
    if shut == 1
        fprintf('Shutdown aborted.\n');
    end
    if shut == 2
        fprintf('Hibernation aborted.\n');
    end
    shut = 0;
end

close(hndl);

% Initialize variables
Phase_screen_init;

params_flag = 0;        % Unflag Phase_screen_params.m script
%% RUN SIMULATION
% Loop over turbulence strengths
for idxg = 1 : leng
    % Check for cancel button press
    if abort
        break
    end
        
    fprintf('Turbulence strength: %u of %u\n', idxg, leng);
    h = waitbar(0, '0%', 'Name', 'Simulating turbulence...', ...
         'CreateCancelBtn', 'setappdata(gcbf,''canceling'',1)');
    setappdata(h,'canceling',0);
    % Loop over realizations
    for idxreal = 1 : nreals
                        
        % Loop over phase screens
        for idxscr = ztmin_idx : ztmax_idx
            if isinf(r0(idxg, idxscr))
                phz(:,:,idxscr) = 0;
            else
                [phz_lo, phz_hi] = ft_sh_phase_screen( ...
                    r0(idxg, idxscr), N, delta(idxscr), L0, l0);
                phz(:,:,idxscr) = phz_lo + phz_hi;
                if coinc        % Is it a coincidence simulation?
                    if invers   % Is there momentum inversion?
                        phz(:,:,idxscr) = phz(:,:,idxscr) + ...
                            rot90(phz(:,:,idxscr),2);
                    else
                        phz(:,:,idxscr) = 2*phz(:,:,idxscr);
                    end
                end
            end
        end
    
        % Simulate turbulent propagation
        [xn yn Uout_real] = ang_spec_multi_prop(Uin, wvl, ...
            delta1, deltan, z, sg.*exp(1i*phz));
        
        % Derive intensity pattern
        Iout_real = abs(Uout_real).^2;
        Iout_real = Iout_real/sum(Iout_real(:));
        
        % Incoherent sum of patterns from different realizations
        Iout(:,:,idxg) = Iout(:,:,idxg) + Iout_real;
        
        % Check for cancel button press
        if getappdata(h,'canceling')
            abort = 1;
            break
        end
        
        waitbar(idxreal/nreals, h, sprintf('%3.1f%%',idxreal/nreals*100));
       
    end
    % Output intensity for given gamma
    Iout(:,:,idxg) = Iout(:,:,idxg)/nreals;
    delete(h);
end

%% PLOT RESULTS (if computation was not aborted)

if ~abort
% Beam transverse profiles at observation plane ('measured gammas only')
figure(1);
set(gcf, 'Units', 'normalized', 'Position', [0.625 0.28 0.365 0.61]);

for g0i = 1 : leng0
    subplot( ceil(sqrt(leng0)), ceil(sqrt(leng0)), g0i);
    if ~isempty(find(Iout(:, :, g0indx(g0i)),1))
        imagesc(Iout(:, :, g0indx(g0i)));
        axis off;
        title(['\gamma = ', num2str(g0(g0i)*1e6), ' \mum']);
    end
end

% Power going through slit
if coinc
    for idxg = 1 : leng
    A = Iout(:,:,idxg);
    Pslit(idxg) = CoincSlitIntegrate(A,a);
    end
else
    for idxg = 1 : leng
        A = Iout(:,:,idxg);
        Pslit(idxg) = SlitIntegrate(A,a);
    end
end

Pslit = Pslit/Pslit(1);
figure(2);
set(gcf, 'Units', 'normalized', 'Position', [2/3 0.05 1/3 1/4]);
plot(g, Pslit, 'LineWidth', 1.5);
end

%% EXPORT RESULTS (only if simulation was completed)

if ~abort
    % File name and folder
    folder_save = fullfile(user_folder, 'Data', 'Computation', ...
        'Phase_screen', date);
    clock_str = clock;
    date_str = [num2str(clock_str(1)), '-', num2str(clock_str(2)), ...
        '-', num2str(clock_str(3))];
    time_str = [num2str(clock_str(4)), ':', num2str(clock_str(5)), ...
        ':', num2str(round(clock_str(6)))];

    if ~exist(folder_save, 'dir')
        mkdir(folder_save)
        curr_file = 1;
    else
        prev_files = dir(fullfile(folder_save, '*.dat'));
        if isempty(prev_files)
            curr_file = 1;
        else
            num_str = regexp(prev_files(end).name, '#(\d+).', 'tokens');
            curr_file = str2num(num_str{1}{1});
            curr_file = curr_file + 1;
        end
    end

    exp_file = [date_str, ' PSdata #', sprintf('%03u',curr_file),'.dat'];
    fid = fopen(fullfile(folder_save,exp_file), 'w');

    % File info and model parameters
    header = ['File saved at ', time_str, ' on the ', date_str];
    fprintf(fid, '%s\n', header);
    
    fprintf(fid, 'Simulation for: %s\n', choice);

    data_params = char('Wavelength', 'Beam waist', 'Slit width', ...
        'Inner scale', 'Outer scale', 'Grid spacing (source plane)', ...
        'Grid spacing (obs. plane)', 'Grid size', ...
        'Num. of realizations', 'Prop. distance', 'Num. of screens', ...
        'Screen intervals', 'Distance between planes', ...
        'Num. of prop. planes');

    fprintf(fid, '%s\n\n', '<< LENGTH UNITS: m >>');
    format_spec = [data_params(1,:), '\t%3.2e\n', ...
                data_params(2,:), '\t%3.1e\n', ...
                data_params(3,:), '\t%3.1e\n', ...
                data_params(4,:), '\t%3.1e\n', ...
                data_params(5,:), '\t%3.1e\n', ...
                data_params(6,:), '\t%3.1e\n', ...
                data_params(7,:), '\t%3.1e\n', ...
                data_params(8,:), '\t%u\n', ...
                data_params(9,:), '\t%u\n', ...
                data_params(10,:), '\t%3.1e\n', ...
                data_params(11,:), '\t%u\n', ...
                data_params(12,:), '\t(%3.1e - %3.1e)\n', ...
                data_params(13,:), '\t%3.1e\n', ...
                data_params(14,:), '\t%u\n'];
    
    fprintf(fid, format_spec, wvl, wn, a0, ...
        l0, L0, delta1, deltan, N, nreals, Dz, ...
        nscr, zscr_min, zscr_max, mean(z(2:n)-z(1:n-1)),nz);
    fprintf(fid, '\nPower over slit\n');
    fprintf(fid, '%u\t', Pslit);
    
    % Close file
    fclose(fid);
end

%% SHUTDOWN COMPUTER?
wtime = 90;         % Time limit to abort shutdown

if shut
    if shut == 1
        shut_mode = 'shutdown';
        shut_cmd = 'shutdown /s';
    else
        shut_mode = 'hibernation';
        shut_cmd = 'shutdown /h';
    end
    
    fprintf('Computer %s enabled.\n', shut_mode);
    tic;
    reply_shut = 'k';
    while (~sum(strcmpi(reply_shut, {'y', 'n', ''})) && toc < wtime)
        if shut == 1
        reply_shut = waitinput('Proceed with shutdown? Y/N [Y]:    ', ...
            wtime, 's');
        else
        reply_shut = waitinput('Proceed with hibernation? Y/N [Y]:    ', ...
            wtime, 's');
        end
        fprintf('\n');
    end
    if (isempty(reply_shut) || ~sum(strcmpi(reply_shut, {'y', 'n'})))
        reply_shut = 'Y';
    end
    if strcmpi(reply_shut,'y')
        if shut == 1
            fprintf('\nShutting down.\n');
        else
            fprintf('\nHibernating.\n');
        end
        pause(3);
        shut = 0;
        system(shut_cmd);
    end
    if shut == 1
        fprintf('Shutdown aborted.\n');
    elseif shut == 2
        fprintf('Hibernation aborted.\n');
    end
end
