function [state, options, optchanged] = ga_output_function(options, state, flag)
    optchanged = false;
    persistent history_x;
    
    if strcmp(flag, 'iter')
        % Store parameter values at each iteration
        history_x = [history_x; state.Generation, state.Population];
    end

    % Plot the parameters after each generation
    if state.Generation == options.MaxGenerations
        figure;
        for i = 1:size(history_x, 2)  % Loop through each design variable
            subplot(2, 2, i);
            plot(history_x(:, 1), history_x(:, i+1), '-o');
            xlabel('Generation');
            ylabel(['x' num2str(i)]);
            title(['Parameter ' num2str(i) ' Evolution']);
        end
    end
end
