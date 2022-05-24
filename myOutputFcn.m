
function stop = myOutputFcn(x,optimValues,state,displayflag)
    stop = false;
    switch state
        case 'iter'
            % Make updates to plot or guis as needed
            xPhys = reshape(x, nely, nelx, nelz);
            %% PRINT RESULTS
            fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f ch.:%7.3f\n',optimValues.iteration,optimValues.fval, ...
                mean(xPhys(:)),optimValues.stepsize);
            %% PLOT DENSITIES
            if displayflag, figure(10); clf; display_3D(xPhys); end
            title([' It.:',sprintf('%5i',optimValues.iteration),...
                ' Obj. = ',sprintf('%11.4f',optimValues.fval),...
                ' ch.:',sprintf('%7.3f',optimValues.stepsize)]);
        case 'init'
            % Setup for plots or guis
            if displayflag
                figure(10)
            end
        case 'done'
            % Cleanup of plots, guis, or final plot
            figure(10); clf; display_3D(xPhys);
        otherwise
    end % switch
end % myOutputFcn
