function plotFirstPassageCondition(dataSet, varargin)

% ----- parse options -----
p = inputParser;
addParameter(p,'FPthresh',[10,30,60],@(c)isnumeric(c));
addParameter(p,'XVar',"",@(x)ischar(x)||(isstring(x)&&isscalar(x)));
addParameter(p,'YVar',"",@(x)ischar(x)||(isstring(x)&&isscalar(x)));
addParameter(p,'MarkerSize',8,@(x)isnumeric(x)&&isscalar(x)&&x>0);
parse(p,varargin{:});
opt = p.Results;

Xvar = opt.XVar;
Yvar = opt.YVar;
FPthresh = opt.FPthresh;
mkSize= opt.MarkerSize;

visualization.plotFiguresPreamble;
N = numel(FPthresh);

%% 
cycles = dataSet.cycle(dataSet.cycle.IGNORED==0,:);

% ---- figure ----
f = figure();
f.Units = 'centimeters';
f.Position = [0,0,6,6];
t = tiledlayout(N,1, ...
    "TileSpacing","none","Padding","compact");

for n = 1:N
    nexttile;
    th = FPthresh(n);
    fp = visualization.compressFirstPassage(cycles, th / dataSet.FrameToMin);
    scatter(fp.(Xvar), fp.(Yvar),mkSize,'MarkerEdgeColor','flat');
    ylabel(sprintf('FP:%d',th));
    set(gca,'YScale','log');
    if n ~= N
        xticklabels([]);
    end
end

xlabel(t,sprintf("%s",Xvar),'Interpreter','none');
ylabel(t,sprintf("%s",Yvar),'Interpreter','none');

