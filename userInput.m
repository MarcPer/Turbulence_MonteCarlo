classdef userInput<handle
    properties
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
            userInput.shutDownAborted(isAbort, shutCtrl);
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
    end
    
end