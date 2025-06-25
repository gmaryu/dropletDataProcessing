function [] = plotDuration(expID, dmin, dmax, YLim, minN, newFig)
visualization.plotFiguresPreamble;

cropWidthCM = 0.45;  % defines axis width

%%
if newFig
    figure;
else
    gcf;
    clf;
end

%data = load(sprintf("%s.mat", expID)).data;
data = expID;
data.info.IGNORED = zeros(size(data.info,1),1);
data.info.SPERM_COUNT(isnan(data.info.SPERM_COUNT)) = 0;

idxCont = data.info(bitand(data.info.IGNORED == 0, bitand(bitand(data.info.SPERM_COUNT == 0, data.info.MEDIAN_DIAMETER > dmin), data.info.MEDIAN_DIAMETER < dmax)), :);
ptmCont = nan * zeros(size(idxCont, 1), 99);
for i = 1:size(idxCont)
    pid = idxCont.POS_ID(i);
    did = idxCont.TRACK_ID(i);
    pcycle = data.cycle(bitand(data.cycle.POS_ID == pid, data.cycle.TRACK_ID == did), :);
    peakTime = [pcycle.START_FRAME; pcycle.END_FRAME(end)]' * data.FrameToMin;
    ptmCont(i, 1:length(peakTime)) = peakTime;
end
ptmCont = ptmCont(:, any(~isnan(ptmCont)));

idx1s = data.info(bitand(data.info.IGNORED == 0, bitand(bitand(data.info.SPERM_COUNT == 1, data.info.MEDIAN_DIAMETER > dmin), data.info.MEDIAN_DIAMETER < dmax)), :);
ptm1s = nan * zeros(size(idx1s, 1), 99);
for i = 1:size(idx1s)
    pid = idx1s.POS_ID(i);
    did = idx1s.TRACK_ID(i);
    pcycle = data.cycle(bitand(data.cycle.POS_ID == pid, data.cycle.TRACK_ID == did), :);
    peakTime = [pcycle.START_FRAME; pcycle.END_FRAME(end)]' * data.FrameToMin;
    ptm1s(i, 1:length(peakTime)) = peakTime;
end
ptm1s = ptm1s(:, any(~isnan(ptm1s)));

dptmCont = diff(ptmCont, 1, 2);
dptm1s = diff(ptm1s, 1, 2);

ax1 = axes();
set(ax1, axSpec{:});
hold on;

x = ptmCont(:, 1:end - 1);
x = x(:);
y = dptmCont(:);
set(ax1, "YLim", YLim);
plot(x, y, "o", "Color", colBlue, "MarkerSize", 3, "LineWidth", LINEWIDTH, "MarkerFaceColor", "w");
hold on;
lowessCont = malowess(x(bitand(~isnan(x), ~isnan(y))), y(bitand(~isnan(x), ~isnan(y))), "Span", 0.25);
xCont = x(bitand(~isnan(x), ~isnan(y)));

x = ptm1s(:, 1:end - 1);
x = x(:);
y = dptm1s(:);
plot(ax1, x, y, "o", "Color", colRed, "MarkerSize", 3, "LineWidth", LINEWIDTH, "MarkerFaceColor", "w");
lowess1s = malowess(x(bitand(~isnan(x), ~isnan(y))), y(bitand(~isnan(x), ~isnan(y))), "Span", 0.25);
x1s = x(bitand(~isnan(x), ~isnan(y)));

[x, i] = sort(xCont);
plot(x, lowessCont(i), "k-", "LineWidth", LINEWIDTH * 5, "Color", [1, 1, 1, 0.75]);
plot(x, lowessCont(i), "k-", "LineWidth", LINEWIDTH, "Color", colBlue);
[x, i] = sort(x1s);
plot(x, lowess1s(i), "k-", "LineWidth", LINEWIDTH * 5, "Color", [1, 1, 1, 0.75]);
plot(x, lowess1s(i), "k-", "LineWidth", LINEWIDTH, "Color", colRed);

set(ax1, "Position", [1.5, 3, 2.5, 2.5]);
xlabel("Cycle start time (min)");
ylabel("Cycle duration (min)");

set(ax1, "XGrid","off","YGrid","on");

%
ax2 = axes();
set(ax2, axSpec{:});
hold on;

nvis = sum(sum(~isnan(dptm1s)) >= minN);
set(ax2, "XLim", [0 + 0.2, nvis + 1 - 0.2], "YLim", YLim);
set(ax2, "XTick", [])
set(ax2, "Position", [5.5, 3, 2 * cropWidthCM * (nvis + 1), 2.5]);

rng(0);
gap = 0.2;
scatterwidth = 0.3;

nvalid0 = sum(~isnan(dptmCont));
nvalid1 = sum(~isnan(dptm1s));

textypos = min([dptmCont(:); dptm1s(:)]) / 2;

tx = text(ax2, 1 - gap, textypos, num2str(nvalid0(1)), "FontSize", 6, "HorizontalAlignment", "center");
txext = tx.Extent;
txextr = txext(1) + txext(3);
tx.String = sprintf("n = %d", nvalid0(1));
tx.HorizontalAlignment = "right";
txext2 = tx.Extent;
txextr2 = txext2(1) + txext2(3);
txpos = tx.Position;
txpos(1) = txpos(1) + (txextr - txextr2);
tx.Position = txpos;
text(ax2, 1 + gap, textypos, num2str(nvalid1(1)), "FontSize", 6, "HorizontalAlignment", "center");

for i = 1:nvis
    x = (rand(length(dptmCont(:, i)), 1) - 0.5) * scatterwidth;
    y1 = dptmCont(:, i);
    plot(ax2, x + i - gap, y1, ".", "Color", colBlue);
    hold on;
    x = (rand(length(dptm1s(:, i)), 1) - 0.5) * scatterwidth;
    y2 = dptm1s(:, i);
    plot(ax2, x + i + gap, y2, ".", "Color", colRed);
    if i > 1
        text(ax2, i - gap, textypos, num2str(nvalid0(i)), "FontSize", 6, "HorizontalAlignment", "center");
        text(ax2, i + gap, textypos, num2str(nvalid1(i)), "FontSize", 6, "HorizontalAlignment", "center");
    end
    %boxchart(ax2, (i - gap) * ones(1, length(y1)), y1, "MarkerStyle", "none", "BoxFaceColor", "w", "LineWidth", LINEWIDTH, "BoxLineColor", "k", "BoxWidth", gap * 1.5, "BoxFaceAlpha", 0.75);
    %boxchart(ax2, (i + gap) * ones(1, length(y2)), y2, "MarkerStyle", "none", "BoxFaceColor", "w", "LineWidth", LINEWIDTH, "BoxLineColor", "k", "BoxWidth", gap * 1.5, "BoxFaceAlpha", 0.75);
    boxchart(ax2, (i - gap) * ones(1, length(y1)), y1, "MarkerStyle", "none", "BoxFaceColor", "w", "LineWidth", LINEWIDTH, "BoxWidth", gap * 1.5, "BoxFaceAlpha", 0.75);
    boxchart(ax2, (i + gap) * ones(1, length(y2)), y2, "MarkerStyle", "none", "BoxFaceColor", "w", "LineWidth", LINEWIDTH, "BoxWidth", gap * 1.5, "BoxFaceAlpha", 0.75);
end

for i = 1:nvis
    y1 = dptmCont(:, i);
    y2 = dptm1s(:, i);

    ymax = max([y1; y2]);
    h = (ax2.YLim(2) - ax2.YLim(1)) * 0.05;
    plot([i - gap, i - gap, i + gap, i + gap], [ymax - h * 0.5, ymax, ymax, ymax - h * 0.5] + h * 2, "k", "LineWidth", LINEWIDTH);
    star = visualization.getStatsStars(y1, y2);
    if strcmp(star, "n.s.")
        text(ax2, i, ymax + 3 * h, star, "FontSize", 6, "HorizontalAlignment", "center");
    else
        text(ax2, i, ymax + 2.5 * h, star, "FontSize", 6, "HorizontalAlignment", "center");
    end
end

for i = 1:nvis
    if mod(i, 2) == 0
        h = patch([i - 0.5, i - 0.5, i + 0.5, i + 0.5], [ax2.YLim(1), ax2.YLim(2), ax2.YLim(2), ax2.YLim(1)], "k", "FaceAlpha", 0.1, "LineStyle", "none");
        uistack(h, "bottom");
    end
end

set(ax2, "XGrid","off","YGrid","on");

%
ax3 = axes();
set(ax3, axSpec{:});
hold on;

set(ax3, "XLim", [0 + 0.2, nvis + 1 - 0.2], "YLim", [0, 1]);
set(ax3, "Position", [5.5, 2.1, 2 * cropWidthCM * (nvis + 1), 1]);
set(ax3, "XTick", 1:nvis);
for i = 1:nvis
    plot(ax3, [i - 1.5 * gap, i + 1.5 * gap], [0.5, 0.5], "k", "LineWidth", LINEWIDTH);
    text(ax3, i - gap, 0.7, "0", "FontSize", 6, "HorizontalAlignment", "center");
    text(ax3, i + gap, 0.7, "1", "FontSize", 6, "HorizontalAlignment", "center");
    text(ax3, i, 0.3, num2str(i), "FontSize", 6, "HorizontalAlignment", "center");
    hold on;
end

axis off;


