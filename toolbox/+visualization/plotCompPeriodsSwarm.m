function plotCompPeriodsSwarm(groupData, varargin)
% plotCompPeriodsSwarm - Swarm (and optional box) plot by Group × Category
%                       with alternating gray category bands.
%
% INPUT
%   groupData : 1×Ng struct array with fields:
%       .values  : (cell OR numeric matrix). If matrix, NaN is treated as missing.
%                  'obs-by-cat' => rows: observations, cols: categories (default)
%                  'cat-by-obs' => rows: categories,  cols: observations
%       .cats    : (optional) category IDs (numeric or string/cellstr).
%                  If omitted, 1..nCat is used after normalization.
%       .name    : (optional) group label. Auto-filled as 'G1','G2',...
%
% NAME-VALUE OPTIONS
%   'ShowBox'     (false) : overlay boxchart
%   'Gap'         (0.25)  : within-category horizontal spacing
%   'Jitter'      (0.30)  : swarm jitter width
%   'Colors'      ([]   ) : Ng×3 RGB; defaults to lines(Ng). If shorter, it repeats.
%   'CatOrder'    ([]   ) : explicit category order (numeric or string/cellstr)
%   'MatrixOrient'('obs-by-cat'|'cat-by-obs')
%
% NOTE
% - Categories may be numeric OR string. Internally we normalize and match
%   with '==' (numeric) or strcmp (string) accordingly.
% - Missing values in matrix inputs must be NaN (they are dropped).

% ---------- parse options ----------
p = inputParser;
addParameter(p,'ShowBox',false,@islogical);
addParameter(p,'Gap',0.25,@(x)isnumeric(x)&&isscalar(x));
addParameter(p,'Jitter',0.30,@(x)isnumeric(x)&&isscalar(x));
addParameter(p,'Colors',[],@(c)isnumeric(c)&&size(c,2)==3);
addParameter(p,'CatOrder',[],@(x)isnumeric(x)||isstring(x)||iscellstr(x));
addParameter(p,'MatrixOrient','obs-by-cat',@(s)ischar(s)||isstring(s));
parse(p,varargin{:});
opt = p.Results;

Ng = numel(groupData);
if Ng==0, warning('groupData is empty. Nothing to plot.'); return; end

% ---------- default colors (repeat if needed) ----------
if isempty(opt.Colors)
    opt.Colors = lines(Ng);
elseif size(opt.Colors,1) < Ng
    % repeat rows to cover all groups
    reps = ceil(Ng/size(opt.Colors,1));
    opt.Colors = repmat(opt.Colors, reps, 1);
end

% ---------- normalize inputs: matrix -> cell per category; cats autofill ----------
for g = 1:Ng
    V = groupData(g).values;
    % check field existence on the g-th element (clearer than isfield on the whole array)
    hasCats = isfield(groupData(g),'cats') && ~isempty(groupData(g).cats);

    if isnumeric(V)
        % Convert matrix into cell-of-vectors per category, dropping NaNs
        switch lower(string(opt.MatrixOrient))
            case "obs-by-cat"
                nCat = size(V,2);
                C = arrayfun(@(j) V(~isnan(V(:,j)), j), 1:nCat, 'uni', 0);
            case "cat-by-obs"
                nCat = size(V,1);
                C = arrayfun(@(i) V(i, ~isnan(V(i,:))).', 1:nCat, 'uni', 0);
            otherwise
                error('MatrixOrient must be "obs-by-cat" or "cat-by-obs".');
        end
        groupData(g).values = C;
        if ~hasCats, groupData(g).cats = 1:nCat; end
    else
        % values already cell: ensure cats exists and matches length
        nCat = numel(V);
        if ~hasCats, groupData(g).cats = 1:nCat; end
    end

    if ~isfield(groupData(g),'name') || isempty(groupData(g).name)
        groupData(g).name = sprintf('G%d', g);
    end
end

% ---------- union of categories across groups (numeric OR string) ----------
allCats = [];
isCatString = false;
for g = 1:Ng
    c = groupData(g).cats;
    if iscell(c) || isstring(c)
        c = cellstr(string(c));             % normalize to cellstr
        isCatString = true;
    end
    if isempty(allCats)
        allCats = c;
    else
        allCats = union(allCats, c, 'stable');
    end
end

% Optional explicit order: keep provided order, append any leftover categories
if ~isempty(opt.CatOrder)
    want = opt.CatOrder;
    if iscell(want) || isstring(want), want = cellstr(string(want)); isCatString = true; end
    present = intersect(want, allCats, 'stable');
    rest    = setdiff(allCats, present, 'stable');
    allCats = [present, rest];
end

% If any group used string-like cats, ensure master is cellstr
if isCatString && ~(iscell(allCats))
    allCats = cellstr(string(allCats));
end
Ncat = numel(allCats);
if Ncat==0, warning('No categories to plot.'); return; end


% ---------- plotting ----------
f = figure; 
f.Units = 'centimeters';                      % set physical units for reproducible size
f.Position = [2, 2, max(6, Ncat), 5];         % width grows with #categories (min width 6 cm)
ax = axes('Parent', f); hold(ax,'on');
centers = linspace(-opt.Gap, opt.Gap, Ng);    % within-category offsets for groups
groupLegendHandles = gobjects(Ng,1);

for g = 1:Ng
    cats = groupData(g).cats;
    if iscell(cats) || isstring(cats), cats = cellstr(string(cats)); end

    for j = 1:numel(cats)
        % locate x-index of this category in the master category list
        if isCatString
            idx = find(strcmp(allCats, cats{j}), 1);
        else
            idx = find(allCats == cats(j), 1);
        end
        if isempty(idx), continue; end

        x0 = idx + centers(g);
        yv = groupData(g).values{j};
        if isempty(yv), continue; end

        % draw swarm points
        if ~isgraphics(groupLegendHandles(g))
            s = swarmchart(ax, x0*ones(numel(yv),1), yv, 10, ...
                'filled', 'MarkerFaceColor', opt.Colors(g,:), ...
                'MarkerEdgeColor','none', ...
                'DisplayName', groupData(g).name, ...      
                'HandleVisibility','on');                  
            groupLegendHandles(g) = s;                     
        else
            s = swarmchart(ax, x0*ones(numel(yv),1), yv, 10, ...
                'filled', 'MarkerFaceColor', opt.Colors(g,:), ...
                'MarkerEdgeColor','none', ...
                'HandleVisibility','off');                 
        end
        s.XJitter = 'randn'; s.XJitterWidth = opt.Jitter;

        % optional box overlay (transparent white box with black edges)
        if opt.ShowBox
            boxchart(ax, x0*ones(size(yv)), yv, ...
                'BoxFaceColor','w','BoxFaceAlpha',0.75, ...
                'LineWidth',1,'MarkerStyle','none', ...
                'HandleVisibility','off');                
        end
    end
end

% axis limits and alternating gray bands
xlim(ax, [0.5, Ncat+0.5]);
yl = ylim(ax);
for i = 1:Ncat
    if mod(i,2)==0
        h = patch(ax, [i-0.5 i-0.5 i+0.5 i+0.5], [yl(1) yl(2) yl(2) yl(1)], ...
            'k', 'FaceAlpha', 0.08, 'LineStyle', 'none');
        uistack(h, 'bottom');   % send band behind the data
    end
end

% ticks/labels
set(ax, 'XTick', 1:Ncat);
if isCatString
    set(ax, 'XTickLabel', allCats);
end
box(ax, 'on');
grid(ax, 'off');
ax.Layer = 'top';
xlabel(ax, 'Cycles');
ylabel(ax, 'Period Length (min)');
legend(ax, groupLegendHandles(isgraphics(groupLegendHandles)), ...
       'Location','bestoutside');
end
