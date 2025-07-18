function [] = plotSizeDependence(expID, newFig)
visualization.plotFiguresPreamble;

widthcm = 4;  % defines axis width
heightcm = 2.5;  % defines axis width
horgap = 1.5;
vergap = 1.5;
YLim1 = [0, 400];
YLim2 = [0.3, 0.65];
YLim3 = [0, 0.1];
YLim4 = [0, 1e-7];

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

idx1s = data.info(bitand(data.info.IGNORED == 0, data.info.SPERM_COUNT == 1), :);
cycle = table();
for i = 1:size(idx1s)
    pid = idx1s.POS_ID(i);
    did = idx1s.TRACK_ID(i);
    cycle = [cycle; data.cycle(bitand(data.cycle.POS_ID == pid, data.cycle.TRACK_ID == did), :)];
end

cycle = visualization.calcCycleAdditionalFeatures(cycle, data.FrameToMin, data.PixelToUm);
refmarkers = (log10(visualization.convertDiameterToVolume([80, 120, 160, 200, 240])) - 5.25) * 3;
cycle.MARKERSIZE(cycle.MARKERSIZE < 0.5, :) = 0.5;

maxc = max(unique(cycle.CYCLE_ID));

%
ax1 = axes();
set(ax1, axSpec{:});
hold on;
set(ax1, "Position", [1, 1 + heightcm + vergap, widthcm, heightcm]);

%
ax2 = axes();
set(ax2, axSpec{:});
hold on;
set(ax2, "Position", [1 + widthcm + horgap, 1 + heightcm + vergap, widthcm, heightcm]);

%
ax3 = axes();
set(ax3, axSpec{:});
hold on;
set(ax3, "Position", [1, 1, widthcm, heightcm]);

%
ax4 = axes();
set(ax4, axSpec{:});
hold on;
set(ax4, "Position", [1 + widthcm + horgap, 1, widthcm, heightcm]);

%
axs = {ax1, ax2, ax3, ax4};
YLim = {YLim1, YLim2, YLim3, YLim4};
fields = {"DURATION", "DNA_INC_RATE_COEFF", "NCVR", "DNACR"};
ylabels = {"Cycle duration (min)", "Hoechst acc. rate (a.u.)", "N/C volume ratio", "Hoechst/C ratio (a.u.)"};

for i = 1:length(axs)
    for j = 1:min(7, maxc)
        cyc = cycle(cycle.CYCLE_ID == j, :);
        cyc = sortrows(cyc, "MARKERSIZE", "descend");
        for k = 1:size(cyc, 1)
            plot(axs{i}, cyc(k, :).START_MINUTE, cyc(k, :).(fields{i}), "o", "Color", colcbrewer(j, :), "MarkerSize", cyc(k, :).MARKERSIZE, "LineWidth", LINEWIDTH, "MarkerFaceColor", "w");
        end
    end
    cyc = cycle(cycle.CYCLE_ID > 7, :);
    cyc = sortrows(cyc, "MARKERSIZE", "descend");
    for k = 1:size(cyc, 1)
        plot(axs{i}, cyc(k, :).START_MINUTE, cyc(k, :).(fields{i}), "o", "Color", [0.75, 0.75, 0.75], "MarkerSize", cyc(k, :).MARKERSIZE, "LineWidth", LINEWIDTH, "MarkerFaceColor", "w");
    end
    set(axs{i}, "YLim", YLim{i}, "TickDir", "out");
    
    xlabel(axs{i}, "Cycle start time (min)");
    ylabel(axs{i}, ylabels{i});
end

ax0 = axes();
set(ax0, axSpec{:});
hold on;
set(ax0, "Position", [1 + 2 * (widthcm + horgap), 1, 2, 1]);

for i = 1:length(refmarkers)
    plot(i, 1, "o", "Color", "k", "MarkerSize", refmarkers(i), "LineWidth", LINEWIDTH, "MarkerFaceColor", "w");
end
axis off;
