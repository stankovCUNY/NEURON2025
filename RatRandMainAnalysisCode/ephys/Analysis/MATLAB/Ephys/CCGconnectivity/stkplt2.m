function [bias, scale] = stkplt(action, addButtons)
%STKPLT Stack plots vertically in the current axes
%   [BIAS, SCALE] = STKPLT stacks or unstacks plots in the current axes
%   [BIAS, SCALE] = STKPLT('init') toggles between stacked and unstacked plots (default)
%   [BIAS, SCALE] = STKPLT('init', 1) also adds audio playback buttons
%   STKPLT('play') plays the audio for the selected plot
%
%   Outputs:
%     BIAS - offsets applied to each plot
%     SCALE - scaling factors applied to each plot
%
%   Original Author: Bert deVries, late 1990s
%   Updated for MATLAB R2024b

    % Default action
    if nargin < 1
        action = 'init';
    end

    % Find valid plots (exclude vertical lines)
    lines = findobj(gca, 'Type', 'line');
    plots = [];
    
    for i = 1:length(lines)
        xData = get(lines(i), 'XData');
        if ~isempty(xData) && any(xData ~= xData(1))
            plots = [plots, lines(i)];
        end
    end
    
    numPlots = length(plots);
    
    % Process actions
    if strcmp(action, 'play')
        handlePlayAction(plots, numPlots);
    else
        [bias, scale] = handleStackAction(plots, numPlots, action, nargin > 1 && addButtons == 1);
    end
end

function handlePlayAction(plots, numPlots)
    % Play the selected signal as audio
    playButtons = findobj(gcf, 'Type', 'uicontrol', 'String', 'Play');
    
    for i = 1:length(playButtons)
        if get(playButtons(i), 'Value')
            set(playButtons(i), 'Value', 0);
            y = get(plots(numPlots + 1 - i), 'YData');
            soundsc(y - mean(y));
        end
    end
end

function [bias, scale] = handleStackAction(plots, numPlots, ~, addButtons)
    % Handle stacking/unstacking plots
    
    % Initialize outputs
    bias = zeros(numPlots, 1);
    scale = zeros(numPlots, 1);
    
    % Create single group of plots
    groupIDs = ones(1, numPlots);
    numGroups = max(groupIDs);
    
    % Calculate spacing parameters
    vertSpace = 0.1 / numGroups;
    plotHeight = (1 - (numGroups + 1) * vertSpace) / numGroups;
    
    % Create app data if it doesn't exist
    if isempty(getappdata(gca, 'StackState'))
        setappdata(gca, 'StackState', 'off');
    end
    
    % Toggle stack state
    if strcmp(getappdata(gca, 'StackState'), 'off')
        % Stack the plots
        for i = 1:numPlots
            % Get data
            y = get(plots(i), 'YData');
            
            % Calculate range, skipping NaN values
            yMax = max(y(~isnan(y)));
            yMin = min(y(~isnan(y)));
            
            % Calculate position
            offset = groupIDs(i) * vertSpace + (groupIDs(i) - 1) * plotHeight;
            
            % Scale to fit in allocated space
            scale(i) = plotHeight / (yMax - yMin + eps); % Add eps to prevent division by zero
            bias(i) = offset - scale(i) * yMin;
            
            % Update plot
            set(plots(i), 'YData', scale(i) * y + bias(i), 'UserData', [scale(i), bias(i)]);
        end
        
        % Update axis properties
        setappdata(gca, 'StackState', 'on');
        set(gca, 'YTick', []);
        
        % Add play buttons if requested
        if addButtons
            addPlayButtons(plots, numPlots, bias);
        end
    else
        % Unstack the plots
        deletePlayButtons();
        
        for i = 1:numPlots
            % Get stored scaling information
            stackInfo = get(plots(i), 'UserData');
            
            if ~isempty(stackInfo)
                % Extract scaling parameters
                plotScale = stackInfo(1);
                plotBias = stackInfo(2);
                
                % Get data and restore original values
                y = get(plots(i), 'YData');
                set(plots(i), 'YData', (y - plotBias) / plotScale, 'UserData', []);
            end
        end
        
        % Update axis properties
        setappdata(gca, 'StackState', 'off');
        set(gca, 'YTickMode', 'auto');
    end
end

function addPlayButtons(plots, numPlots, bias)
    % Add play buttons for each plot
    
    % Get figure and axes properties
    set(gcf, 'Units', 'normalized');
    axisPos = get(gca, 'Position');
    yLim = get(gca, 'YLim');
    
    % Create play buttons
    for i = 1:numPlots
        % Calculate button position
        buttonY = axisPos(2) + (bias(i) - yLim(1)) / (yLim(2) - yLim(1)) * axisPos(4) - 0.035;
        
        % Create button with modern style
        uicontrol('Style', 'togglebutton', ...
                 'String', 'Play', ...
                 'Callback', @(~,~) stkplt('play'), ...
                 'Units', 'normalized', ...
                 'Position', [0.01, buttonY, 0.09, 0.07], ...
                 'FontWeight', 'bold');
    end
end

function deletePlayButtons()
    % Remove play buttons
    delete(findobj(gcf, 'Style', 'togglebutton', 'String', 'Play'));
end