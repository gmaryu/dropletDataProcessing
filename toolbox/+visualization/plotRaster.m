function [props, peakMatrix] = plotRaster(dataSet, varargin)
% plotRaster - Raster plot with optional separator lines
%
% Usage:
%   [props, peakMatrix] = plotRaster(dataSet);                          % without lines
%   [props, peakMatrix] = plotRaster(dataSet, 'line_option', true);     % with separator lines
%   [props, peakMatrix] = plotRaster(dataSet, 'line_option', true, ...
%                                    'MarkerSize', 10, 'Color', 'k');
%
% Name-Value options:
%   'line_option' (logical) : Draw separator lines for accumulated droplet counts per position (default: false)
%   'Axes'        (axes)    : Target axes handle for plotting (default: gca)
%   'MarkerSize'  (scalar)  : Marker size for scatter (default: 10)
%   'lineColor'       (char/str): Color for lines (default: 'k')

    % ---- Parse inputs ----
    p = inputParser;
    addRequired(p, 'dataSet', @(s)isstruct(s) && isfield(s,'cycle') && isfield(s,'info') && isfield(s,'FrameToMin'));
    addParameter(p, 'line_option', false, @(x)islogical(x) || isnumeric(x));
    addParameter(p, 'Axes', [], @(ax) isempty(ax) || isa(ax,'matlab.graphics.axis.Axes'));
    addParameter(p, 'MarkerSize', 5, @(x)isnumeric(x) && isscalar(x) && x>0);
    addParameter(p, 'lineColor', 'k', @(c)ischar(c) || isstring(c) || (isnumeric(c) && numel(c)==3));
    addParameter(p, 'newFig', true, @(x)islogical(x) || isnumeric(x));
    parse(p, dataSet, varargin{:});

    % ---- Target Figure Pannel ----
    newFig      = p.Results.newFig;
    if newFig, f = figure; else, hold on; end
    line_option = logical(p.Results.line_option);
    ax          = p.Results.Axes;
    if isempty(ax), ax = gca; end
    msize       = p.Results.MarkerSize;
    colorSpec   = p.Results.lineColor;

    % add IGNORED column in case tables doesn't have the column name
    % (accept all doplets info)
    if ~ismember("IGNORED", dataSet.info.Properties.VariableNames)
        disp("IGNORED is temporally added to dataSet.Info");
        dataSet.info.IGNORED = zeros(size(dataSet.info,1),1);
    end
    if ~ismember("IGNORED", dataSet.cycle.Properties.VariableNames)
        disp("IGNORED is temporally added to dataSet.cycle");
        dataSet.cycle.IGNORED = zeros(size(dataSet.cycle,1),1);
    end
    

    % ---- Initialize return values ----
    peakMatrix = [];
    vec_pos = [];
    vec_trackid = [];
    vec_maxCycles = [];
    vec_medianDiameter = [];

    count = 0;
    max_cyc     = max(dataSet.cycle.CYCLE_ID(dataSet.cycle.IGNORED == 0));  % maximum cycle number in this condition
    total_peaks = dataSet.cycle;                % peak information (table expected)
    info        = dataSet.info;                 % additional information (table expected)
    unique_pos  = unique(total_peaks.POS_ID);   % list of unique position IDs
    count_vector = nan(size(unique_pos));       % cumulative droplet counts for each position (for separator lines)

    % ---- Construct peak matrix ----
    for ppos = 1:numel(unique_pos)
        tmp_cycle_pos = total_peaks(total_peaks.POS_ID == unique_pos(ppos) & total_peaks.IGNORED == 0, :);
        tmp_info_pos  = info(info.POS_ID == unique_pos(ppos), :);
        unique_droplets = unique(tmp_cycle_pos.TRACK_ID);

        for d = 1:numel(unique_droplets)
            tmp_droplet = tmp_cycle_pos(tmp_cycle_pos.TRACK_ID == unique_droplets(d), :);
            tmp_info_droplet = tmp_info_pos(tmp_info_pos.TRACK_ID == unique_droplets(d), :);

            % NaN vector with length = max_cyc+1
            tmp_peaks = nan(max_cyc+1, 1);

            % Collect unique start/end frames
            allpeaks = unique([tmp_droplet.START_FRAME; tmp_droplet.END_FRAME]);
            nput = min(numel(allpeaks), max_cyc+1);
            if nput > 0
                tmp_peaks(1:nput) = allpeaks(1:nput);
            end

            % Append droplet as a row in peakMatrix
            peakMatrix = [peakMatrix; tmp_peaks.']; %#ok<AGROW>

            % Store droplet attributes
            vec_pos            = [vec_pos; unique_pos(ppos)]; %#ok<AGROW>
            vec_trackid        = [vec_trackid; unique_droplets(d)]; %#ok<AGROW>
            vec_maxCycles      = [vec_maxCycles; tmp_info_droplet.CYCLE_NUMBER + 1]; %#ok<AGROW>
            vec_medianDiameter = [vec_medianDiameter; tmp_info_droplet.ORIGINAL_MED_DIAMETER]; %#ok<AGROW>
        end

        % Update droplet count for this position
        count = count + numel(unique_droplets);
        count_vector(ppos) = count;
    end

    % ---- Construct props output ----
    props.POS_ID         = vec_pos;
    props.TRACK_ID       = vec_trackid;
    props.maxCycles      = vec_maxCycles;
    props.medianDiameter = vec_medianDiameter;

    % ---- Plot ----
    axes(ax); hold(ax, 'on');

    nRows = size(peakMatrix, 1);
    nCols = size(peakMatrix, 2);

    % Scatter plot for each column (cycle)
    % NaNs will be ignored automatically by scatter
    for c = 1:nCols
        t = peakMatrix(:, c) * dataSet.FrameToMin;    % frame -> minutes
        y = (1:nRows).';                              % droplet index
        scatter(ax, t, y, msize, 'o', 'filled');
    end
    ylim(ax, [1, nRows]);
    xlabel(ax, 'Time (min)');

    % Draw separator lines if enabled
    if line_option
        yline(ax, count_vector, '--', 'Color', colorSpec);
    end

    % ---- Figure Properties for export ----
    f.Units = 'centimeters';
    f.Position = [2,2,9,6];

    hold(ax, 'off');
end
