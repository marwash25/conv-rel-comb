function [x, info] = FISTA(fx, gradf, g, proxg, parameter)
    
    if parameter.verbose
        fprintf('%s\n', repmat('*', 1, 68));
    end
    
    % Set the clock.
    time1       = tic;
    timestart   = toc(time1);
    
    % Initialize x, y and t.
    x_next       = parameter.x0;
    y_next       = parameter.x0;
    t_next       = 1;
    F_new        = fx(parameter.x0)+g(parameter.x0);
    
    % Main loop.
    for iter = 1:parameter.maxit
        
        x           = x_next;      
        y           = y_next;
        t           = t_next;
        F_old       = F_new;
        % Compute error and save data to be plotted later on.
        info.itertime(iter ,1)  = toc(time1) - timestart;
        info.fx(iter, 1)        = F_old;
                
        if parameter.verbose
            % Print the information.
            fprintf('Iter = %4d, f(x) = %5.3e\n', ...
                    iter, info.fx(iter, 1));
        end
        % Start the clock.
        timestart   = toc(time1);
        
        % Update next iteration
        x_next = proxg( y - (1/parameter.Lips) * gradf(y)); %should be proxg(upd,L) (defined as prox_g(u,L)=arg min .5||u-x||_2^2 + 1/L g(x))
         
        % evaluate the new value of f(x_{k+1}).
        F_new       = fx(x_next)+g(x_next); 
        
        if(parameter.AR) %if we want to do Adaptive restart
       
        % Compare the old_f(x) and new_f(x) to decide to restart or not.
        
        if( F_old < F_new ) 
            y = x; t= 1;
            x_next = proxg( y - (1/parameter.Lips) * gradf(y));
        end
        
        end
        
        t_next = .5*(1+ sqrt(4*t^2+1));
        gamma  = (t - 1)/t_next;
        y_next = x_next + gamma *(x_next - x);
      
        

        % Check stopping criterion.
        if (abs(F_new - F_old)/F_old) <= parameter.tolx 
            break;
        end

    end

    % Finalization.
    info.iter           = iter;
    info.time           = cumsum(info.itertime);
    info.totaltime      = info.time(iter);
    
end
%**************************************************************************
% END OF THE IMPLEMENTATION.
%**************************************************************************
