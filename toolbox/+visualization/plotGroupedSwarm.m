function plotGroupedSwarm(groupData, varargin)
% plotGroupedSwarm - Swarm (and optional box) plot by Group × Category.
%                    Legend shows one entry per group (color key).
%
% INPUT
%   groupData : 1×Ng struct array with fields:
%       .values : table containing at least columns:
%                 - CYCLE_ID (category id)
%                 - IGNORED  (logical or 0/1; rows with IGNORED==1 are skipped)
%                 - the measurement column specified by 'columnName'
%       .cats   : category ids present in this group (numeric or string/cellstr)
%       .name   : (optional) group label for legend; auto 'G1','G2',...
%
% NAME-VALUE OPTIONS
%   'columnName' (required) : char/string scalar; table variable to plot
%   'ShowBox'    (false)    : overlay boxchart for each (group,category)
%   'Gap'        (0.25)     : within-category horizontal spacing between groups
%   'Jitter'     (0.30)     : swarm jitter width
%   'Colors'     ([]     )  : Ng×3 RGB; defaults to lines(Ng). Repeats if short.
%   'CatOrder'   ([]     )  : explicit category order (numeric or string/cellstr)
%
% NOTE
% - Handles numeric and string categories. Missing categories per group are allowed.
% - Only the first swarm object per group participates in the legend.

% -------- parse options --------
p = inputParser;
addParameter(p,'columnName',"", @(x)ischar(x) || (isstring(x) && isscalar(x)));
addParameter(p,'ShowBox',false,@islogical);
addParameter(p,'Gap',0.25,@(x)isnumeric(x)&&isscalar(x));
addParameter(p,'Jitter',0.30,@(x)isnumeric(x)&&isscalar(x));
addParameter(p,'Colors',[],@(c)isnumeric(c)&&size(c,2)==3);
addParameter(p,'CatOrder',[],@(x)isnumeric(x)||isstring(x)||iscellstr(x));
addParameter(p,'MaxCatNum',[],@(x)isnumeric(x)&&isscalar(x));
parse(p,varargin{:});
opt = p.Results;

Ng = numel(groupData);
if Ng==0, warning('groupData is empty.'); return; end

catNum = opt.MaxCatNum;

% columnName is required
colName = string(opt.columnName);
if strlength(colName)==0
    error('columnName is required (char or string scalar).');
end

% default colors; repeat if shorter than Ng
if isempty(opt.Colors)
    opt.Colors = lines(Ng);
elseif size(opt.Colors,1) < Ng
    opt.Colors = repmat(opt.Colors, ceil(Ng/size(opt.Colors,1)), 1);
end

% -------- build measurement cell per (group,category) --------
for g = 1:Ng
    % group label for legend
    if ~isfield(groupData(g),'name') || isempty(groupData(g).name)
        groupData(g).name = sprintf('G%d', g);
    end

    % safety: ensure IGNORED exists and is logical
    v = groupData(g).values;
    if ~ismember("IGNORED", v.Properties.VariableNames)
        v.IGNORED = false(height(v),1);
    elseif ~islogical(v.IGNORED)
        v.IGNORED = logical(v.IGNORED);
    end

    % ensure requested column exists
    if ~ismember(colName, string(v.Properties.VariableNames))
        error('Column "%s" does not exist in group %d table.', colName, g);
    end

    % make output cell aligned to listed categories
    c = groupData(g).cats;
    groupData(g).dnaInt = cell(numel(c),1);

    % extract the measurement column once
    colData = v.(colName);

    % fill per category (skip IGNORED rows)
    for i = 1:numel(c)
        % logical mask is clearer than bitand
        mask = (~v.IGNORED) & (v.CYCLE_ID == c(i));

        % remove 0 for log scale plot
        vals = colData(mask);
        ind0 = (vals == 0);
        vals(ind0) = [];
        
        groupData(g).dnaInt{i} = vals;
        %groupData(g).dnaInt{i} = colData(mask);
        
    end

    % write back in case we modified v
    groupData(g).values = v;
end

% -------- union of categories across groups (numeric OR string) --------
allCats = [];
isCatString = false;
for g = 1:Ng
    c = groupData(g).cats;
    if iscell(c) || isstring(c)
        c = cellstr(string(c));      % normalize to cellstr
        isCatString = true;
    end
    if isempty(allCats), allCats = c; else, allCats = union(allCats, c, 'stable'); end
end

% optional explicit order
if ~isempty(opt.CatOrder)
    want = opt.CatOrder;
    if iscell(want) || isstring(want), want = cellstr(string(want)); isCatString = true; end
    present = intersect(want, allCats, 'stable');
    rest    = setdiff(allCats, present, 'stable');
    allCats = [present, rest];
end

% keep a consistent type for category master list
if isCatString && ~iscell(allCats)
    allCats = cellstr(string(allCats));
end

if isempty(catNum)
    Ncat = numel(allCats);
else
    Ncat = catNum;
end
if Ncat==0, warning('No categories to plot.'); return; end

% -------- plotting --------
f = figure; ax = axes('Parent', f); hold(ax,'on');
f.Units = 'centimeters';
f.Position = [2, 2, max(6, Ncat), 5];      % width scales with #categories
centers = linspace(-opt.Gap, opt.Gap, Ng); % per-group horizontal offsets

% store one handle per group to appear in legend
groupLegendHandles = gobjects(Ng,1);

for g = 1:Ng
    cats = groupData(g).cats;
    if iscell(cats) || isstring(cats), cats = cellstr(string(cats)); end

    for j = 1:Ncat
    %for j = 1:numel(cats)
        % locate category index in the master list
        if isCatString
            idx = find(strcmp(allCats, cats{j}), 1);
        else
            idx = find(allCats == cats(j), 1);
        end
        if isempty(idx), continue; end

        x0 = idx + centers(g);
        yv = groupData(g).dnaInt{j};
        if isempty(yv), continue; end

        % first series of each group: visible in legend; others: hidden
        baseArgs = {ax, x0*ones(numel(yv),1), yv, 10, ...
            'filled', 'MarkerFaceColor', opt.Colors(g,:), 'MarkerEdgeColor','none'};
        if ~isgraphics(groupLegendHandles(g))
            s = swarmchart(baseArgs{:}, ...
                'DisplayName', groupData(g).name, ...
                'HandleVisibility','on');
            groupLegendHandles(g) = s;
        else
            s = swarmchart(baseArgs{:}, 'HandleVisibility','off');
        end
        s.XJitter = 'randn'; s.XJitterWidth = opt.Jitter;

        % optional box overlay (hidden from legend)
        if opt.ShowBox
            boxchart(ax, x0*ones(size(yv)), yv, ...
                'BoxFaceColor','k','BoxFaceAlpha',0.35, ...
                'LineWidth',1,'MarkerStyle','none', ...
                'HandleVisibility','off');
        end
    end
end

% axis, bands, labels
xlim(ax, [0.5, Ncat+0.5]); yl = ylim(ax);
for i = 1:Ncat
    if mod(i,2)==0
        h = patch(ax, [i-0.5 i-0.5 i+0.5 i+0.5], [yl(1) yl(2) yl(2) yl(1)], ...
            'k', 'FaceAlpha', 0.08, 'LineStyle', 'none');
        uistack(h, 'bottom');
    end
end
set(ax,'XTick',1:Ncat);
if isCatString, set(ax,'XTickLabel',allCats); end
box(ax,'on'); grid(ax,'off'); ax.Layer = 'top';
xlabel(ax,'Cycles');
ylabel(ax, sprintf('%s (a.u.)', colName), 'Interpreter','none');   % label from columnName

% legend with one entry per group
legend(ax, groupLegendHandles(isgraphics(groupLegendHandles)), 'Location','bestoutside');

end
