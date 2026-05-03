function f = plotDurationStartTime(dataSet, varargin)
% plotDurationStartTime - Scatter plot for cell cycle periods vs start time
%  
% Usage:
%   f = plotDurationStartTime(dataSet);                 
%   f = plotDurationStartTime(dataSet, 'newFig', false);     % to olverlay different conditions
%   f = plotDurationStartTime(dataSet, 'newFig', false, ...
%                                    'MarkerSize', 5, 'Color', 'r');
%
% Name-Value options:
%   'newFig' (logical) : Create a new figure panel (default: true)
%   'Axes'        (axes)    : Target axes handle for plotting (default: gca)
%   'MarkerSize'  (scalar)  : Marker size for scatter (default: 3)
%   'Color'       (char/str): Color for markers/lines (default: 'k')
    
    visualization.plotFiguresPreamble;
    
    % ---- Parse inputs ----
    p = inputParser;
    addRequired(p, 'data', @(s)isstruct(s) && isfield(s,'cycle') && isfield(s,'info') && isfield(s,'FrameToMin'));
    addParameter(p, 'newFig', true, @(x)islogical(x) || isnumeric(x));
    addParameter(p, 'Axes', [], @(ax) isempty(ax) || isa(ax,'matlab.graphics.axis.Axes'));
    addParameter(p, 'MarkerSize', 3, @(x)isnumeric(x) && isscalar(x) && x>0);
    addParameter(p, 'Color', 'k', @(c)ischar(c) || isstring(c) || (isnumeric(c) && numel(c)==3));
    parse(p, dataSet, varargin{:});
    LINEWIDTH = 2;
    
    % ---- Target Figure Pannel ----
    newFig      = p.Results.newFig;
    if newFig, f = figure; else, hold on; end
    ax          = p.Results.Axes;
    if isempty(ax), ax = gca; end
    msize       = p.Results.MarkerSize;
    colorSpec   = p.Results.Color;
    
    % ---- Figure Properties for export ----
    f.Units = 'centimeters';
    f.Position = [2,2,5,5];
    
    % ---- Data preprocessing ----
    % add IGNORED column in case the table doesn't have the column (accept
    % all droplets)
    if ~ismember("IGNORED", dataSet.info.Properties.VariableNames)
        disp("IGNORED is added to dataSet.Info");
        dataSet.info.IGNORED = zeros(size(dataSet.info,1),1); 
    end
    
    % ---- Period length calculation ----
    idx = dataSet.info(dataSet.info.IGNORED == 0, :);
    ptm = nan * zeros(size(idx, 1), 99);
    for i = 1:size(idx)
        pid = idx.POS_ID(i);
        did = idx.TRACK_ID(i);
        pcycle = dataSet.cycle(bitand(dataSet.cycle.POS_ID == pid, dataSet.cycle.TRACK_ID == did), :);
        peakTime = [pcycle.START_FRAME; pcycle.END_FRAME(end)]' * dataSet.FrameToMin;
        ptm(i, 1:length(peakTime)) = peakTime;
    end
    ptm = ptm(:, any(~isnan(ptm)));
    
    % period time
    dptm = diff(ptm, 1, 2);
    
    % ---- Draw a scatter plot ----
    x = ptm(:, 1:end - 1);
    x = x(:);
    y = dptm(:);
    plot(x, y, "o", "Color", colorSpec, "MarkerSize", msize, "LineWidth", LINEWIDTH, "MarkerFaceColor", "w");
    hold on;
    lowessCont = malowess(x(bitand(~isnan(x), ~isnan(y))), y(bitand(~isnan(x), ~isnan(y))), "Span", 0.25);
    xCont = x(bitand(~isnan(x), ~isnan(y)));
    
    % ---- Draw a trend line ----
    [x, i] = sort(xCont);
    plot(x, lowessCont(i), "k-", "LineWidth", LINEWIDTH * 5, "Color", [1, 1, 1, 0.75]);
    plot(x, lowessCont(i), "k-", "LineWidth", LINEWIDTH, "Color", colorSpec);
    
    xlabel("Cycle start time (min)");
    ylabel("Cycle duration (min)");
    set(ax, "XGrid","off","YGrid","on");

end


