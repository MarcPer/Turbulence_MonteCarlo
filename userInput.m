classdef UserInput
    methods(Static)
        function shutdownComputer(enumShutdown)
            wtime = 90;         % Time limit to abort shutdown
            
            if ~(enumShutdown)
                return;
            end
            
            if enumShutdown == 1
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
                reply_shut = waitinput(sprintf('Proceed with %s? Y/N [Y]:\n', shut_mode), ...
                    wtime, 's');
            end
            if (isempty(reply_shut) || ~sum(strcmpi(reply_shut, {'y', 'n'})))
                reply_shut = 'Y';
            end
            
            switch lower(reply_shut)
                case 'y'
                    if enumShutdown == 1
                        fprintf('\nShutting down.\n');
                    elseif enumShutdown == 2
                        fprintf('\nHibernating.\n');
                    end
                    pause(3);
                    system(shut_cmd);
                otherwise
                    if enumShutdown == 1
                        fprintf('Shutdown aborted.\n');
                    elseif enumShutdown == 2
                        fprintf('Hibernation aborted.\n');
                    end
             end
 
        end

        function isAbort = abortWhenConstraintFail(fail, shutCtrl)
            isAbort = 0;
            if ~fail
                fprintf('Press any key to continue...\n');
                pause;
                return
            end
            switch( questdlg('Abort?','Constraint analysis failed', 'Yes', 'No', 'Yes'));
                case 'Yes'
                    isAbort = 1;
            end
            UserInput.shutDownAborted(isAbort, shutCtrl);
        end

        function shutDownAborted(isAbort, shutCtrl)
            if isAbort
                switch shutCtrl
                    case 1
                        fprintf('Shutdown aborted...\n');
                    case 2
                        fprintf('Hibernation aborted...\n');
                end
            end
        end

        function h = createWaitBar()
            h = waitbar(0, '0%', 'Name', 'Simulating turbulence...', ...
                'CreateCancelBtn', 'setappdata(gcbf,''canceling'',1)');
            setappdata(h,'canceling',0);
        end

        function updateWaitBar(waitBarHandle,iRe, nRe)
            waitbar(iRe/nRe, waitBarHandle, sprintf('%3.1f%%',iRe/nRe*100)); 
        end

        function abort = isAborted(waitBarHandle)
            abort = 0;
            if getappdata(waitBarHandle,'canceling')
                abort = 1;
            end
        end

        function printOutProgress(statusString, idx, fIdx)
            fprintf('%s: %u of %u\n', statusString, idx, fIdx);
        end

        function shutDown = confirmShutdownOrHibernate(enumShutdown)
            shutDown = enumShutdown;

            switch enumShutdown
                case 0
                    return;
                case 1
                    shutString = 'SHUTDOWN';
                case 2
                    shutString = 'HIBERNATION';
                otherwise
                    shutString = 'HIBERNATION';
            end

            prompt = sprintf('Computer is scheduled for %s. Press CANCEL to unschedule.\n', shutString);
            choiceShut = questdlg(prompt, 'Confirm Shutdown/Hibernation', 'OK', 'Cancel', 'Cancel');

            switch choiceShut
                case 'Cancel'
                    shutDown = 0;
                    fprintf('%s cancelled.\n', shutString);
                case 'OK'
                    fprintf('%s confirmed.\n', shutString);
            end

        end

    end 
end