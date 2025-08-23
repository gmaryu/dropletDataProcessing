function plotScatterByCategoryGroup(groupData, varargin)
% plotScatterByCategoryGroup
%   Tile scatter plots per category; each tile overlays multiple groups with colors.
%   X/Y/Category/Ignore columns are configurable via Name-Value args.
%
% INPUT
%   groupData : 1×Ng struct array with fields:
%       .values : table containing at least CatVar, XVar, YVar (and optionally IgnoreVar)
%       .cats   : categories present in this group (numeric OR string/cellstr)
%       .name   : (optional) group label; defaults to 'G1','G2',...
%
% NAME-VALUE OPTIONS
%   'XVar'          (required) : char/string scalar, table column used as X
%   'YVar'          (required) : char/string scalar, table column used as Y
%   'CatVar'        ('CYCLE_ID'): char/string scalar, category column in table
%   'IgnoreVar'     ('IGNORED') : char/string scalar; rows with true are skipped
%   'Colors'        ([]        ) : Ng×3 RGB; defaults to lines(Ng). Repeats if short.
%   'CatOrder'      ([]        ) : explicit category order (numeric OR string/cellstr)
%   'LowessSpan'    (0.25     ) : [0,1] span for LOWESS smoothing
%   'MarkerSize'    (8        ) : scatter marker size
%   'LegendLocation'('bestoutside') : legend location string
%   'TileLayout'    ('auto'   ) : 'auto' or [nRows nCols] numeric vector
%
% NOTES
% - Works with numeric or string categories; missing categories per group are fine.
% - Legend shows one entry per group (color key), not per category.
% - Rows with NaN in X or Y are dropped automatically.

% ----- parse options -----
p = inputParser;
addParameter(p,'XVar',"",@(x)ischar(x)||(isstring(x)&&isscalar(x)));
addParameter(p,'YVar',"",@(x)ischar(x)||(isstring(x)&&isscalar(x)));
addParameter(p,'CatVar',"CYCLE_ID",@(x)ischar(x)||(isstring(x)&&isscalar(x)));
addParameter(p,'IgnoreVar',"IGNORED",@(x)ischar(x)||(isstring(x)&&isscalar(x)));
addParameter(p,'Colors',[],@(c)isnumeric(c)&&size(c,2)==3);
addParameter(p,'CatOrder',[],@(x)isnumeric(x)||isstring(x)||iscellstr(x));
addParameter(p,'LowessSpan',0.25,@(x)isnumeric(x)&&isscalar(x)&&x>0&&x<=1);
addParameter(p,'MarkerSize',8,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p,'LegendLocation','bestoutside',@(s)ischar(s)||isstring(s));
addParameter(p,'TileLayout','auto',@(v)(ischar(v)||isstring(v)) || (isnumeric(v)&&numel(v)==2));
parse(p,varargin{:});
opt = p.Results;

Ng = numel(groupData);
if Ng==0, warning('groupData is empty.'); return; end

xName   = string(opt.XVar);
yName   = string(opt.YVar);
catName = string(opt.CatVar);
ignName = string(opt.IgnoreVar);

if strlength(xName)==0 || strlength(yName)==0
    error('XVar and YVar are required (char or string scalar).');
end

% ----- colors (repeat if necessary) -----
if isempty(opt.Colors)
    opt.Colors = lines(Ng);
elseif size(opt.Colors,1) < Ng
    opt.Colors = repmat(opt.Colors, ceil(Ng/size(opt.Colors,1)), 1);
end

% ======================================================================
% Precompute per-group per-category (X,Y) cells after filtering
% ======================================================================
for g = 1:Ng
    % group label
    if ~isfield(groupData(g),'name') || isempty(groupData(g).name)
        groupData(g).name = sprintf('G%d', g);
    end

    T = groupData(g).values;    % table
    % required columns check
    req = unique([catName, xName, yName]);
    missing = req(~ismember(req, string(T.Properties.VariableNames)));
    if ~isempty(missing)
        error('Missing required column(s) in group %d: %s', g, strjoin(cellstr(missing),', '));
    end

    % ensure IgnoreVar exists and is logical (create if missing)
    if ~ismember(ignName, string(T.Properties.VariableNames))
        T.(ignName) = false(height(T),1);
    elseif ~islogical(T.(ignName))
        T.(ignName) = logical(T.(ignName));
    end

    % normalize category list type for later matching
    c = groupData(g).cats;
    isStrCatG = iscell(c) || isstring(c);
    if isStrCatG, c = cellstr(string(c)); end
    groupData(g).cats = c;

    % allocate per-category cells
    groupData(g).XCells = cell(numel(c),1);
    groupData(g).YCells = cell(numel(c),1);

    % extract once
    X = T.(xName); Y = T.(yName); C = T.(catName); IG = T.(ignName);

    % fill per category
    for i = 1:numel(c)
        if isStrCatG
            maskCat = strcmp(string(C), string(c{i}));
        else
            maskCat = (C == c(i));
        end
        mask = (~IG) & maskCat & ~isnan(X) & ~isnan(Y);
        groupData(g).XCells{i} = X(mask);
        groupData(g).YCells{i} = Y(mask);
    end

    % write back table (in case we added IgnoreVar)
    groupData(g).values = T;
end

% ======================================================================
% Build master category list (union across groups); optional ordering
% ======================================================================
allCats = [];
isCatString = false;
for g = 1:Ng
    cg = groupData(g).cats;
    if iscell(cg), isCatString = true; end
    if isempty(allCats), allCats = cg; else, allCats = union(allCats, cg, 'stable'); end
end

if ~isempty(opt.CatOrder)
    want = opt.CatOrder;
    if iscell(want) || isstring(want), want = cellstr(string(want)); isCatString = true; end
    present = intersect(want, allCats, 'stable');
    rest    = setdiff(allCats, present, 'stable');
    allCats = [present, rest];
end

if isCatString && ~iscell(allCats)
    allCats = cellstr(string(allCats));
end
Ncat = numel(allCats);
if Ncat==0, warning('No categories to plot.'); return; end

% ======================================================================
% Layout and plotting (tiledlayout)
% ======================================================================
% decide layout
if (ischar(opt.TileLayout) || isstring(opt.TileLayout)) && strcmpi(string(opt.TileLayout),'auto')
    nRows = ceil(sqrt(Ncat));
    nCols = ceil(Ncat / nRows);
else
    nRows = opt.TileLayout(1);
    nCols = opt.TileLayout(2);
end

f = figure; 
f.Units = 'centimeters';
f.Position = [0, 0, max(12, 4*nCols), max(10, 4*nRows)];  % simple sizing heuristic
t = tiledlayout(f, nRows, nCols, 'TileSpacing','compact', 'Padding','compact');

% one legend handle per group
legendHandles = gobjects(Ng,1);

for j = 1:Ncat
    ax = nexttile(t); hold(ax, 'on');

    for g = 1:Ng
        cg = groupData(g).cats;

        % find index of this master category in group's categories
        if isCatString
            idxG = find(strcmp(cg, allCats{j}), 1);
        else
            idxG = find(cg == allCats(j), 1);
        end
        if isempty(idxG), continue; end

        x = groupData(g).XCells{idxG};
        y = groupData(g).YCells{idxG};
        if isempty(x), continue; end

        % scatter points
        h = scatter(ax, x, y, opt.MarkerSize, opt.Colors(g,:), 'filled', ...
            'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5);

        % LOWESS trend line
        yfit = malowess(x, y, "Span", opt.LowessSpan);
        [xs, I] = sort(x);
        plot(ax, xs, yfit(I), '-', 'LineWidth', 1.5, 'Color', opt.Colors(g,:));
    end

    % tile cosmetics
    if isCatString
        title(ax, sprintf('Category: %s', allCats{j}));
    else
        title(ax, sprintf('Category: %g', allCats(j)));
    end
    xlabel(ax, xName); ylabel(ax, yName);
    grid(ax, 'on'); box(ax, 'on');
end


end
