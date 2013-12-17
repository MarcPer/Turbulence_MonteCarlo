classdef UserInput<handle
    properties
        isShutDown;
        isFourthOrder;
        isInverted;
    end
    methods(Static)
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
    end
    
end