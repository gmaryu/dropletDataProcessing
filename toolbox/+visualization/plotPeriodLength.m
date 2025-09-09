function [f, plMat] = plotPeriodLength(dataSet, varargin)

% plotPeriodLength - Swarm plot for cell cycle periods with cycle number
%  
% Usage:
%   f = plotPeriodLength(dataSet);                 
%   f = plotPeriodLength(dataSet, 'MarkerSize', 5, 'Color', 'r');
%
% Name-Value options:
%   'Axes'        (axes)    : Target axes handle for plotting (default: gca)
%   'MarkerSize'  (scalar)  : Marker size for scatter (default: 3)
%   'Color'       (char/str): Color for markers/lines (default: 'k')
   
    visualization.plotFiguresPreamble;
    close all;
    
    % ---- Parse inputs ----
    p = inputParser;
    addRequired(p, 'data', @(s)isstruct(s) && isfield(s,'cycle') && isfield(s,'info') && isfield(s,'FrameToMin'));
    addParameter(p, 'Axes', [], @(ax) isempty(ax) || isa(ax,'matlab.graphics.axis.Axes'));
    addParameter(p, 'MarkerSize', 8, @(x)isnumeric(x) && isscalar(x) && x>0);
    addParameter(p, 'Color', 'k', @(c)ischar(c) || isstring(c) || (isnumeric(c) && numel(c)==3));
    addParameter(p, 'newFig', true, @(x)islogical(x) || isnumeric(x));
    parse(p, dataSet, varargin{:});
    
    % ---- Target Figure Pannel ----
    newFig      = p.Results.newFig;
    if newFig, f = figure; else, hold on; end
    ax          = p.Results.Axes;
    if isempty(ax), ax = gca; end
    msize       = p.Results.MarkerSize;
    colorSpec   = p.Results.Color;
    
    % ---- Data preprocessing ----
    % add IGNORED column in case the table doesn't have the column (accept
    % all droplets)
    if ~ismember("IGNORED", dataSet.info.Properties.VariableNames)
        disp("IGNORED is temporally added to dataSet.Info");
        dataSet.info.IGNORED = zeros(size(dataSet.info,1),1);
    end

    if ~ismember("IGNORED", dataSet.cycle.Properties.VariableNames)
        disp("IGNORED is temporally added to dataSet.Info");
        dataSet.cycle.IGNORED = zeros(size(dataSet.cycle,1),1);
    end
    
    % Peak time collection
    idxCont = dataSet.info(dataSet.info.IGNORED == 0, :);
    maxCycle = max(dataSet.cycle.CYCLE_ID(dataSet.cycle.IGNORED == 0));
    ptmCont = nan * zeros(size(idxCont, 1), maxCycle+1);
    for i = 1:size(idxCont)
        pid = idxCont.POS_ID(i);
        did = idxCont.TRACK_ID(i);
        pcycle = dataSet.cycle((dataSet.cycle.POS_ID == pid & dataSet.cycle.TRACK_ID == did & ...
                    ~dataSet.cycle.IGNORED), :);
        peakTime = [pcycle.START_FRAME; pcycle.END_FRAME(end)]' * dataSet.FrameToMin;
        ptmCont(i, 1:length(peakTime)) = peakTime;
    end 
    
    % Calculate period length from the peak time matrix
    plMat = diff(ptmCont, 1, 2);
    
    % ---- Draw a swarm plot ----
    x=[];
    y=[];
    for c = 1:size(plMat,2)
        x = [x c*ones(1,size(plMat,1))];
        y = [y plMat(:,c)'];
    end
    swarmchart(x,y,msize,colorSpec);
    xlabel('Cycle Number');
    ylabel('Period length (min)');
    set(ax, "XGrid","off","YGrid","on");
    
    % ---- Figure Properties for export ----
    f.Units = 'centimeters';
    f.Position = [2,2,maxCycle+1,5];

end