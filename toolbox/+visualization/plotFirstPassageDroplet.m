function plotFirstPassageDroplet(dataSet, posid, did, FPthresh, varargin)

% ----- parse options -----
p = inputParser;
addParameter(p,'newFig',true,@islogical)
addParameter(p,'IgnoreVar',"IGNORED",@(x)ischar(x)||(isstring(x)&&isscalar(x)));
addParameter(p,'Colors',[], ...
        @(c) (isnumeric(c)   && size(c,2)==3) || ...
             (ischar(c)      && isscalar(c))  || ...
             (isstring(c)    && isscalar(c)));
addParameter(p,'MarkerSize',8,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p,'TileLayout','auto',@(v)(ischar(v)||isstring(v)) || (isnumeric(v)&&numel(v)==2));
parse(p,varargin{:});
opt = p.Results;

newFig = opt.newFig;
ignName = string(opt.IgnoreVar);
col = opt.Colors;

visualization.plotFiguresPreamble;
markerStyle = ['o','^','diamond','*','square'];
N = numel(FPthresh);
if N > length(markerStyle)
    disp('Too much threshold value')
end

%%
if newFig
    f = figure;
    f.Units = 'centimeters';
    f.Position = [0,0,5,8];
else
    gcf;
end

%%

cycles = dataSet.cycle(dataSet.cycle.IGNORED==0,:);
tmpCycles = cycles(cycles.POS_ID == posid & cycles.TRACK_ID == did,:);
lastStart = max(tmpCycles.START_FRAME)*dataSet.FrameToMin;

% xlim
d = floor(log10(lastStart)) + 1;
pos = d - 1;
xlim_r = ceil(lastStart / 10^pos) * 10^pos; % xlim right side


% ---- figure plot ----
hold on

plot(tmpCycles.START_MINUTE, tmpCycles.DURATION - min(tmpCycles.DURATION), "o-", "Color", col, "MarkerSize", 3, "LineWidth", LINEWIDTH, "MarkerFaceColor", "w");

for n = 1:N
    th = FPthresh(n);
    plot([0, xlim_r], [th, th], "Color", [0.75, 0.75, 0.75], "LineWidth", LINEWIDTH);
    fp = visualization.compressFirstPassage(cycles, th / dataSet.FrameToMin);
    plot(fp.START_MINUTE((fp.POS_ID == posid & fp.TRACK_ID == did)), fp.DURATION((fp.POS_ID == posid & fp.TRACK_ID == did)) - min(tmpCycles.DURATION), markerStyle(n), ...
        "Color", "k", "MarkerSize", 5, "LineWidth", LINEWIDTH,'MarkerFaceColor',col);
end

xlabel("Cycle start time (min)");
ylabel("Excess cycle duration (min)");
