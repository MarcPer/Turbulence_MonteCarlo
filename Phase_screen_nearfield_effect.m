%% Clear and close all
close all; clear all;

%% INITIALIZATION CELL
% This script is an adaptation of Phase_screen_main that focuses on the
% effect of separating overlapping phase-screens. It is an important
% consideration when delta-like near-field correlation is not assumed in
% the biphoton, as was done before.
%
% Details on most of the script can be found in the above-referred code.
% <LENGTH UNITS: m>

% Set running folder
user_folder = userpath;
user_folder(end) = [];
cd(user_folder);
cd ../..
user_folder = fullfile(pwd, 'Dropbox');
cd(fullfile(user_folder, 'MATLAB', ...
    'Turb_propagation', 'Phase_screen'));

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

% ------------------------ Code added to this particular script
g0 = 27.9e-6;           % Turbulence strength set to strongest measured
s = sum((1-zt(ztmin_idx:ztmax_idx)/Dz).^(5/3));
r0temp = (0.423 * ktemp^2/(7.75 * Dz^(5/3) * s) * g0.^2).^(-3/5);
% Separation as a multiple of the coherence radius r0temp
rsep = [0; 0.3; 0.6; 1; 1.4] * r0temp ;
lensep = length(rsep);
Iout = zeros([N,N,lensep]);   % Output intensities for all gammas
Pslit = zeros(lensep,1);      % Power integrated over slit (single counts)


params_flag = 0;        % Unflag Phase_screen_params.m script
%% RUN SIMULATION
% Loop over overlapping phase screen separations
for idxsep = 1 : lensep
    % Check for cancel button press
    if abort
        break
    end
        
    fprintf('Phase screen separation: %u of %u\n', idxsep, lensep);
    h = waitbar(0, '0%', 'Name', 'Simulating turbulence...', ...
         'CreateCancelBtn', 'setappdata(gcbf,''canceling'',1)');
    setappdata(h,'canceling',0);
    % Loop over realizations
    for idxreal = 1 : nreals
                        
        % Loop over phase screens
        for idxscr = ztmin_idx : ztmax_idx
            if isinf(r0temp)
                phz(:,:,idxscr) = 0;
            else
                [phz_lo, phz_hi] = ft_sh_phase_screen( ...
                    r0temp, N, delta(idxscr), L0, l0);
                phz(:,:,idxscr) = phz_lo + phz_hi;
                phz(:,:,idxscr) = phz(:,:,idxscr) + ...
                Displace(rot90(phz(:,:,idxscr),2), ...
                round(rsep(idxsep)/delta(idxscr)),0);                
            end
        end
    
        % Simulate turbulent propagation
        [xn yn Uout_real] = ang_spec_multi_prop(Uin, wvl, ...
            delta1, deltan, z, sg.*exp(1i*phz));
        
        % Derive intensity pattern
        Iout_real = abs(Uout_real).^2;
        Iout_real = Iout_real/sum(Iout_real(:));
        
        % Incoherent sum of patterns from different realizations
        Iout(:,:,idxsep) = Iout(:,:,idxsep) + Iout_real;
        
        % Check for cancel button press
        if getappdata(h,'canceling')
            abort = 1;
            break
        end
        
        waitbar(idxreal/nreals, h, sprintf('%3.1f%%',idxreal/nreals*100));
       
    end
    % Output intensity for given separation rsep
    Iout(:,:,idxsep) = Iout(:,:,idxsep)/nreals;
    delete(h);
end

%% PLOT RESULTS (if computation was not aborted)

if ~abort
% Beam transverse profiles at observation plane
figure(1);
set(gcf, 'Units', 'normalized', 'Position', [0.625 0.28 0.365 0.61]);

for rsepi = 1 : lensep
    subplot( ceil(sqrt(lensep)), ceil(sqrt(lensep)), rsepi);
    if ~isempty(find(Iout(:, :, rsepi),1))
        imagesc(Iout(:, :, rsepi));
        axis off;
        title(['Separation = ', num2str(rsep(rsepi))]);
    end
end

% Power going through slit
if coinc
    for idxsep = 1 : lensep
    A = Iout(:,:,idxsep);
    Pslit(idxsep) = CoincSlitIntegrate(A,a);
    end
else
    for idxsep = 1 : lensep
        A = Iout(:,:,idxsep);
        Pslit(idxsep) = SlitIntegrate(A,a);
    end
end

Pslit = Pslit/Pslit(1);
figure(2);
set(gcf, 'Units', 'normalized', 'Position', [2/3 0.05 1/3 1/4]);
plot(rsep, Pslit, 'LineWidth', 1.5);
end

%% EXPORT RESULTS (only if simulation was completed)

if ~abort
    % File name and folder
    folder_save = fullfile(user_folder, 'Data', 'Computation', ...
        'Phase_screen', 'Near-field', date);
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

    exp_file = [date_str, ' PSdata_NF #', sprintf('%03u',curr_file),'.dat'];
    fid = fopen(fullfile(folder_save,exp_file), 'w');

    % File info and model parameters
    header = ['File saved at ', time_str, ' on the ', date_str];
    fprintf(fid, '%s\n', header);
    
    data_params = char('Wavelength', 'Beam waist', 'Slit width', ...
        'Inner scale', 'Outer scale', 'Grid spacing (source plane)', ...
        'Grid spacing (obs. plane)', 'Grid size', ...
        'Num. of realizations', 'Prop. distance', 'Num. of screens', ...
        'Screen intervals', 'Distance between planes', ...
        'Num. of prop. planes', 'Turbulence strength (g)');

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
                data_params(14,:), '\t%u\n', ...
                data_params(15,:), '\t%3.1e\n'];
    
    fprintf(fid, format_spec, wvl, wn, a0, ...
        l0, L0, delta1, deltan, N, nreals, Dz, ...
        nscr, zscr_min, zscr_max, mean(z(2:n)-z(1:n-1)),nz,g0);
    fprintf(fid, '\nScreen separation (r0 units)\t Power over slit\n');
    fprintf(fid, '%3.1e\t%3.2e\n', [rsep'/r0temp; Pslit']);
    
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
