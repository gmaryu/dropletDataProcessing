function makeCytoplasmReference(X_cytoOnly_raw, cytoTraceInfo, cytoPeakTable, resultDataSavePath)

%% Cytoplasm-only peak-to-peak reference construction
%
% Input:
%   X_cytoOnly_raw
%   cytoPeakTable
%   cytoTraceInfo
%
% Required cytoPeakTable columns:
%   PositionNumber
%   DropletNumber
%   CycleInDroplet
%   Peak1Frame
%   Peak2Frame
%   TroughFrame
%
% Required cytoTraceInfo columns:
%   TraceRow
%   PositionNumber
%   DropletNumber
%
% Output:
%   cytoplasm_only_reference_peak2peak.mat

clc; close all;

%% =========================
%  Step 1: Load data
%  =========================

% load("cytoplasm_only_signals.mat", ...
%     "X_cytoOnly_raw", ...
%     "cytoTraceInfo");
% X_cytoOnly_raw = cyt_Cdk1_signal;
% cytoTraceInfo;
% 
% cytoPeakTable = cyt.cycle;

nDroplets = size(X_cytoOnly_raw, 1);
nTime = size(X_cytoOnly_raw, 2);

fprintf("Loaded %d cytoplasm-only droplets and %d time points.\n", ...
    nDroplets, nTime);

mkdir(resultDataSavePath);

%% =========================
%  Step 2: Parameters
%  =========================

featureSmoothWindow = 1;   % 1 = no smoothing

nPhasePoints = 100;
phase = linspace(0, 1, nPhasePoints);

% Cycle QC
minCycleLength = 5;
maxCycleLength = inf;
minAmplitude = 0.02;

% In peak-to-peak cycle, trough should not be too close to either edge
minTroughPhase = 0.05;
maxTroughPhase = 0.95;

% Minimum cycles to build cycle-specific template
minCyclesPerTemplate = 5;

%% =========================
%  Step 3: Create DropletID and join TraceRow
%  =========================

% cytoTraceInfo.DropletID = strcat( ...
%     "P", string(cytoTraceInfo.PositionNumber), ...
%     "_D", string(cytoTraceInfo.DropletNumber));

cytoPeakTable.DropletID = strcat( ...
    "P", string(cytoPeakTable.POS_ID), ...
    "_D", string(cytoPeakTable.TRACK_ID));

% Check whether traceInfo has duplicate IDs
% Gtrace = groupsummary(cytoTraceInfo, "DropletID", "numel");
% 
% if any(Gtrace.GroupCount > 1)
%     warning("Some cytoplasm-only DropletIDs map to multiple TraceRows.");
%     disp(Gtrace(Gtrace.GroupCount > 1, :));
% end

% Join TraceRow onto peak table
cytoBoundaryTable = outerjoin( ...
    cytoPeakTable, ...
    cytoTraceInfo(:, ["DropletID", "TraceRow"]), ...
    "Keys", "DropletID", ...
    "MergeKeys", true);

% Keep original row order if needed
if ~ismember("OriginalRow", cytoBoundaryTable.Properties.VariableNames)
    cytoBoundaryTable.OriginalRow = (1:height(cytoBoundaryTable))';
end

% Check unmatched rows
missingTrace = isnan(cytoBoundaryTable.TraceRow);

fprintf("Cytoplasm-only cycles without matched trace row: %d / %d\n", ...
    sum(missingTrace), height(cytoBoundaryTable));

if any(missingTrace)
    disp(cytoBoundaryTable(missingTrace, ...
        ["PositionNumber", "DropletNumber", "DropletID"]));
end

%% =========================
%  Step 4: Smooth trace matrix without bridging NaNs
%  =========================

X_cytoOnly_feature = smoothTraceMatrixNoBridge( ...
    X_cytoOnly_raw, featureSmoothWindow);

%% =========================
%  Step 5: Segment peak-to-peak cycles
%  =========================

nCyclesTotal = height(cytoBoundaryTable);

cytoCyclesRaw_peak2peak = cell(nCyclesTotal, 1);
cytoCyclesFeature_peak2peak = cell(nCyclesTotal, 1);

TraceRow = nan(nCyclesTotal, 1);
CycleLength = nan(nCyclesTotal, 1);
Peak1Frame = nan(nCyclesTotal, 1);
Peak2Frame = nan(nCyclesTotal, 1);
TroughFrame = nan(nCyclesTotal, 1);
TroughFrameRelative = nan(nCyclesTotal, 1);
TroughPhase = nan(nCyclesTotal, 1);

Peak1Value = nan(nCyclesTotal, 1);
Peak2Value = nan(nCyclesTotal, 1);
TroughValue = nan(nCyclesTotal, 1);
AmplitudePeak2MinusTrough = nan(nCyclesTotal, 1);
AmplitudeMaxMinusMin = nan(nCyclesTotal, 1);

NaNFractionCycle = nan(nCyclesTotal, 1);

for c = 1:nCyclesTotal

    row = cytoBoundaryTable.TraceRow(c);

    idx1 = cytoBoundaryTable.START_INDEX(c);
    idx2 = cytoBoundaryTable.END_INDEX(c);
    troughIdx = cytoBoundaryTable.TROUGH_INDEX(c);
    % idx1 = cytoBoundaryTable.Peak1Frame(c);
    % idx2 = cytoBoundaryTable.Peak2Frame(c);
    % troughIdx = cytoBoundaryTable.TroughFrame(c);

    if isnan(row) || isnan(idx1) || isnan(idx2) || idx2 <= idx1
        continue;
    end

    row = round(row);
    idx1 = round(idx1);
    idx2 = round(idx2);

    if row < 1 || row > size(X_cytoOnly_raw, 1) || ...
            idx1 < 1 || idx2 > size(X_cytoOnly_raw, 2)
        continue;
    end

    frameIdx = idx1:idx2;

    thisRaw = X_cytoOnly_raw(row, frameIdx);
    thisFeature = X_cytoOnly_feature(row, frameIdx);

    cytoCyclesRaw_peak2peak{c} = thisRaw;
    cytoCyclesFeature_peak2peak{c} = thisFeature;

    TraceRow(c) = row;
    CycleLength(c) = numel(frameIdx);
    Peak1Frame(c) = idx1;
    Peak2Frame(c) = idx2;

    Peak1Value(c) = X_cytoOnly_feature(row, idx1);
    Peak2Value(c) = X_cytoOnly_feature(row, idx2);

    if ~isnan(troughIdx)
        troughIdx = round(troughIdx);
        TroughFrame(c) = troughIdx;

        if troughIdx >= idx1 && troughIdx <= idx2
            troughRel = troughIdx - idx1 + 1;
            TroughFrameRelative(c) = troughRel;
            TroughPhase(c) = troughRel / numel(frameIdx);
            TroughValue(c) = X_cytoOnly_feature(row, troughIdx);
        end
    end

    if ~isnan(TroughValue(c))
        AmplitudePeak2MinusTrough(c) = Peak2Value(c) - TroughValue(c);
    end

    AmplitudeMaxMinusMin(c) = ...
        max(thisFeature, [], "omitnan") - min(thisFeature, [], "omitnan");

    NaNFractionCycle(c) = mean(isnan(thisFeature));
end

%% Create cycle info table

cytoCycleInfo_peak2peak = cytoBoundaryTable;
cytoCycleInfo_peak2peak.CycleInDroplet = cytoBoundaryTable.CYCLE_ID;
cytoCycleInfo_peak2peak.TraceRow = TraceRow;
cytoCycleInfo_peak2peak.CycleLength = CycleLength;
cytoCycleInfo_peak2peak.Peak1Frame = Peak1Frame;
cytoCycleInfo_peak2peak.Peak2Frame = Peak2Frame;
cytoCycleInfo_peak2peak.TroughFrame = TroughFrame;
cytoCycleInfo_peak2peak.TroughFrameRelative = TroughFrameRelative;
cytoCycleInfo_peak2peak.TroughPhase = TroughPhase;

cytoCycleInfo_peak2peak.Peak1Value = Peak1Value;
cytoCycleInfo_peak2peak.Peak2Value = Peak2Value;
cytoCycleInfo_peak2peak.TroughValue = TroughValue;
cytoCycleInfo_peak2peak.AmplitudePeak2MinusTrough = AmplitudePeak2MinusTrough;
cytoCycleInfo_peak2peak.AmplitudeMaxMinusMin = AmplitudeMaxMinusMin;
cytoCycleInfo_peak2peak.NaNFractionCycle = NaNFractionCycle;

fprintf("Segmented %d cytoplasm-only peak-to-peak cycles.\n", ...
    height(cytoCycleInfo_peak2peak));

%% =========================
%  Step 6: Cycle-level QC
%  =========================

isValidCycle = ...
    cytoCycleInfo_peak2peak.CycleLength >= minCycleLength & ...
    cytoCycleInfo_peak2peak.CycleLength <= maxCycleLength & ...
    cytoCycleInfo_peak2peak.AmplitudeMaxMinusMin >= minAmplitude & ...
    cytoCycleInfo_peak2peak.TroughPhase >= minTroughPhase & ...
    cytoCycleInfo_peak2peak.TroughPhase <= maxTroughPhase & ...
    cytoCycleInfo_peak2peak.NaNFractionCycle < 0.2;

cytoCycleInfo_peak2peak.IsValidCycle = isValidCycle;

fprintf("Valid cytoplasm-only peak-to-peak cycles: %d / %d\n", ...
    sum(isValidCycle), height(cytoCycleInfo_peak2peak));

figure("Name", "Cytoplasm-only peak-to-peak cycle QC",'Visible','off');

subplot(2, 2, 1);
histogram(cytoCycleInfo_peak2peak.CycleLength);
xlabel("Cycle length, frames");
ylabel("Count");
title("Peak-to-peak cycle length");
grid on;

subplot(2, 2, 2);
histogram(cytoCycleInfo_peak2peak.TroughPhase);
xlabel("Trough phase");
ylabel("Count");
title("Trough position in peak-to-peak cycle");
grid on;

subplot(2, 2, 3);
histogram(cytoCycleInfo_peak2peak.AmplitudeMaxMinusMin);
xlabel("Amplitude, max - min");
ylabel("Count");
title("Cycle amplitude");
grid on;

subplot(2, 2, 4);
boxchart(cytoCycleInfo_peak2peak.CycleInDroplet, ...
         cytoCycleInfo_peak2peak.CycleLength);
xlabel("Cycle number");
ylabel("Cycle length");
title("Cycle length by cycle number");
grid on;

exportgraphics(gcf, fullfile(resultDataSavePath, 'Cytoplasm_cycle_QC.png'));

%% =========================
%  Step 7: Phase and amplitude normalization
%  =========================

cytoCyclesNorm_peak2peak = nan(nCyclesTotal, nPhasePoints);

for c = 1:nCyclesTotal
    y = cytoCyclesFeature_peak2peak{c};

    if isempty(y)
        continue;
    end

    yInterp = phaseResampleAllowNaN(y, phase);
    yNorm = amplitudeNormalize(yInterp);

    cytoCyclesNorm_peak2peak(c, :) = yNorm;
end

%% Plot normalized cycles by cycle number

maxCycleNumber = max(cytoCycleInfo_peak2peak.CycleInDroplet, [], "omitnan");
nPlotCycles = min(maxCycleNumber, 6);

figure("Name", "Cytoplasm-only normalized peak-to-peak cycles",'Visible','off');

for k = 1:nPlotCycles
    idx = cytoCycleInfo_peak2peak.IsValidCycle & ...
          cytoCycleInfo_peak2peak.CycleInDroplet == k;

    Y = cytoCyclesNorm_peak2peak(idx, :);

    subplot(nPlotCycles, 1, k);
    hold on;

    for i = 1:size(Y, 1)
        plot(phase, Y(i, :), "Color", [0.8 0.8 0.8]);
    end

    if ~isempty(Y)
        plot(phase, mean(Y, 1, "omitnan"), ...
            "k-", "LineWidth", 2);
    end

    ylabel(sprintf("C%d", k));
    ylim([-0.05 1.05]);
    grid on;

    if k == 1
        title("Peak-to-peak normalized cytoplasm-only Cdk1 cycles");
    end

    if k == nPlotCycles
        xlabel("Normalized phase, peak-to-peak");
    end
end
exportgraphics(gcf, fullfile(resultDataSavePath, 'Cytoplasm_normalized_cycle.png'));

%% =========================
%  Step 8: Build cycle-specific templates
%  =========================

cytoTemplateByCycle_peak2peak = nan(maxCycleNumber, nPhasePoints);
cytoTemplateMedianByCycle_peak2peak = nan(maxCycleNumber, nPhasePoints);
cytoTemplateSDByCycle_peak2peak = nan(maxCycleNumber, nPhasePoints);
nCyclesPerTemplate_peak2peak = nan(maxCycleNumber, 1);

for k = 1:maxCycleNumber
    idx = cytoCycleInfo_peak2peak.IsValidCycle & ...
          cytoCycleInfo_peak2peak.CycleInDroplet == k;

    Y = cytoCyclesNorm_peak2peak(idx, :);

    nCyclesPerTemplate_peak2peak(k) = size(Y, 1);

    if size(Y, 1) < minCyclesPerTemplate
        continue;
    end

    cytoTemplateByCycle_peak2peak(k, :) = mean(Y, 1, "omitnan");
    cytoTemplateMedianByCycle_peak2peak(k, :) = median(Y, 1, "omitnan");
    cytoTemplateSDByCycle_peak2peak(k, :) = std(Y, 0, 1, "omitnan");
end

figure("Name", "Cytoplasm-only peak-to-peak templates",'Visible','off');
hold on;

for k = 1:maxCycleNumber
    if all(isnan(cytoTemplateByCycle_peak2peak(k, :)))
        continue;
    end

    plot(phase, cytoTemplateByCycle_peak2peak(k, :), ...
        "LineWidth", 1.5, ...
        "DisplayName", sprintf("Cycle %d, n=%d", ...
        k, nCyclesPerTemplate_peak2peak(k)));
end

xlabel("Normalized phase, peak-to-peak");
ylabel("Normalized Cdk1 activity");
title("Cycle-specific cytoplasm-only peak-to-peak Cdk1 templates");
legend("Location", "best");
grid on;
exportgraphics(gcf, fullfile(resultDataSavePath, 'Cytoplasm_cycle_template.png'));

%% =========================
%  Step 9: Estimate cycle-specific cytoplasmic variability
%  =========================

cytoThresholdCorrByCycle_peak2peak = nan(maxCycleNumber, 1);
cytoThresholdRMSEByCycle_peak2peak = nan(maxCycleNumber, 1);

cytoMedianDistanceCorrByCycle_peak2peak = nan(maxCycleNumber, 1);
cytoMedianDistanceRMSEByCycle_peak2peak = nan(maxCycleNumber, 1);

cytoP95DistanceCorrByCycle_peak2peak = nan(maxCycleNumber, 1);
cytoP95DistanceRMSEByCycle_peak2peak = nan(maxCycleNumber, 1);

cytoVariabilityByCycle_peak2peak = cell(maxCycleNumber, 1);

for k = 1:maxCycleNumber
    idx = cytoCycleInfo_peak2peak.IsValidCycle & ...
          cytoCycleInfo_peak2peak.CycleInDroplet == k;

    cycleIdx = find(idx);
    Y = cytoCyclesNorm_peak2peak(idx, :);

    nY = size(Y, 1);

    if nY < minCyclesPerTemplate
        continue;
    end

    dCorr = nan(nY, 1);
    dRMSE = nan(nY, 1);

    for i = 1:nY
        y = Y(i, :);

        otherIdx = setdiff(1:nY, i);

        if isempty(otherIdx)
            continue;
        end

        templateLOO = mean(Y(otherIdx, :), 1, "omitnan");

        r = corr(y(:), templateLOO(:), "Rows", "complete");
        dCorr(i) = 1 - r;

        dRMSE(i) = sqrt(mean((y - templateLOO).^2, "omitnan"));
    end

    Tvar = table;
    Tvar.CycleNumber = repmat(k, nY, 1);
    Tvar.CycleIndex = cycleIdx(:);
    Tvar.PositionNumber = cytoCycleInfo_peak2peak.POS_ID(cycleIdx);
    Tvar.DropletNumber = cytoCycleInfo_peak2peak.TRACK_ID(cycleIdx);
    Tvar.DropletID = cytoCycleInfo_peak2peak.DropletID(cycleIdx);
    Tvar.DistanceCorr = dCorr;
    Tvar.DistanceRMSE = dRMSE;

    cytoVariabilityByCycle_peak2peak{k} = Tvar;

    cytoMedianDistanceCorrByCycle_peak2peak(k) = median(dCorr, "omitnan");
    cytoMedianDistanceRMSEByCycle_peak2peak(k) = median(dRMSE, "omitnan");

    cytoP95DistanceCorrByCycle_peak2peak(k) = prctile(dCorr, 95);
    cytoP95DistanceRMSEByCycle_peak2peak(k) = prctile(dRMSE, 95);

    cytoThresholdCorrByCycle_peak2peak(k) = ...
        median(dCorr, "omitnan") + 3 * mad(dCorr, 1);

    cytoThresholdRMSEByCycle_peak2peak(k) = ...
        median(dRMSE, "omitnan") + 3 * mad(dRMSE, 1);
end

%% Summary table

cytoReferenceSummary_peak2peak = table;
cytoReferenceSummary_peak2peak.CycleNumber = (1:maxCycleNumber)';
cytoReferenceSummary_peak2peak.Ncycles = nCyclesPerTemplate_peak2peak;
cytoReferenceSummary_peak2peak.MedianDistanceCorr = ...
    cytoMedianDistanceCorrByCycle_peak2peak;
cytoReferenceSummary_peak2peak.P95DistanceCorr = ...
    cytoP95DistanceCorrByCycle_peak2peak;
cytoReferenceSummary_peak2peak.ThresholdCorr_MedianPlus3MAD = ...
    cytoThresholdCorrByCycle_peak2peak;
cytoReferenceSummary_peak2peak.MedianDistanceRMSE = ...
    cytoMedianDistanceRMSEByCycle_peak2peak;
cytoReferenceSummary_peak2peak.P95DistanceRMSE = ...
    cytoP95DistanceRMSEByCycle_peak2peak;
cytoReferenceSummary_peak2peak.ThresholdRMSE_MedianPlus3MAD = ...
    cytoThresholdRMSEByCycle_peak2peak;

disp("Cytoplasm-only peak-to-peak reference summary:");
disp(cytoReferenceSummary_peak2peak);

%% Plot thresholds

figure("Name", "Cytoplasm-only peak-to-peak variability thresholds",'Visible','off');

subplot(2, 1, 1);
hold on;
plot(cytoReferenceSummary_peak2peak.CycleNumber, ...
     cytoReferenceSummary_peak2peak.MedianDistanceCorr, ...
     "-o", "LineWidth", 1.5);
plot(cytoReferenceSummary_peak2peak.CycleNumber, ...
     cytoReferenceSummary_peak2peak.P95DistanceCorr, ...
     "-o", "LineWidth", 1.5);
plot(cytoReferenceSummary_peak2peak.CycleNumber, ...
     cytoReferenceSummary_peak2peak.ThresholdCorr_MedianPlus3MAD, ...
     "-o", "LineWidth", 1.5);
xlabel("Cycle number");
ylabel("1 - Pearson correlation");
title("Peak-to-peak cytoplasmic waveform variability");
legend("Median", "95th percentile", "Median + 3 MAD", "Location", "best");
grid on;

subplot(2, 1, 2);
hold on;
plot(cytoReferenceSummary_peak2peak.CycleNumber, ...
     cytoReferenceSummary_peak2peak.MedianDistanceRMSE, ...
     "-o", "LineWidth", 1.5);
plot(cytoReferenceSummary_peak2peak.CycleNumber, ...
     cytoReferenceSummary_peak2peak.P95DistanceRMSE, ...
     "-o", "LineWidth", 1.5);
plot(cytoReferenceSummary_peak2peak.CycleNumber, ...
     cytoReferenceSummary_peak2peak.ThresholdRMSE_MedianPlus3MAD, ...
     "-o", "LineWidth", 1.5);
xlabel("Cycle number");
ylabel("RMSE");
title("Peak-to-peak cytoplasmic RMSE variability");
legend("Median", "95th percentile", "Median + 3 MAD", "Location", "best");
grid on;

exportgraphics(gcf, fullfile(resultDataSavePath, 'Cytoplasm_cycle_threshold.png'));
%% =========================
%  Step 10: Save output
%  =========================

outputFile = "cytoplasm_only_reference_peak2peak.mat";

save(outputFile, ...
    "X_cytoOnly_raw", ...
    "X_cytoOnly_feature", ...
    "cytoTraceInfo", ...
    "cytoPeakTable", ...
    "cytoBoundaryTable", ...
    "cytoCycleInfo_peak2peak", ...
    "cytoCyclesRaw_peak2peak", ...
    "cytoCyclesFeature_peak2peak", ...
    "cytoCyclesNorm_peak2peak", ...
    "cytoTemplateByCycle_peak2peak", ...
    "cytoTemplateMedianByCycle_peak2peak", ...
    "cytoTemplateSDByCycle_peak2peak", ...
    "nCyclesPerTemplate_peak2peak", ...
    "cytoThresholdCorrByCycle_peak2peak", ...
    "cytoThresholdRMSEByCycle_peak2peak", ...
    "cytoMedianDistanceCorrByCycle_peak2peak", ...
    "cytoMedianDistanceRMSEByCycle_peak2peak", ...
    "cytoP95DistanceCorrByCycle_peak2peak", ...
    "cytoP95DistanceRMSEByCycle_peak2peak", ...
    "cytoVariabilityByCycle_peak2peak", ...
    "cytoReferenceSummary_peak2peak", ...
    "phase", ...
    "nPhasePoints", ...
    "featureSmoothWindow", ...
    "minCycleLength", ...
    "maxCycleLength", ...
    "minAmplitude", ...
    "minTroughPhase", ...
    "maxTroughPhase");

fprintf("Saved peak-to-peak cytoplasm-only reference to %s\n", outputFile);

%% =========================
%  Local functions
%  =========================

function Xsmooth = smoothTraceMatrixNoBridge(Xraw, smoothWindow)

    Xsmooth = nan(size(Xraw));

    for i = 1:size(Xraw, 1)
        y = Xraw(i, :);

        if smoothWindow <= 1
            Xsmooth(i, :) = y;
            continue;
        end

        validIdx = ~isnan(y);

        d = diff([false, validIdx, false]);
        segStart = find(d == 1);
        segEnd = find(d == -1) - 1;

        for s = 1:numel(segStart)
            idx = segStart(s):segEnd(s);

            if numel(idx) < smoothWindow
                Xsmooth(i, idx) = y(idx);
            else
                Xsmooth(i, idx) = smoothdata(y(idx), ...
                    "movmedian", smoothWindow);
            end
        end
    end
end

function yInterp = phaseResampleAllowNaN(y, targetPhase)

    y = y(:)';
    n = numel(y);

    if n < 2 || all(isnan(y))
        yInterp = nan(size(targetPhase));
        return;
    end

    oldPhase = linspace(0, 1, n);
    validIdx = ~isnan(y);

    if sum(validIdx) < 2
        yInterp = nan(size(targetPhase));
        return;
    end

    yInterp = interp1(oldPhase(validIdx), y(validIdx), ...
        targetPhase, "linear", NaN);
end

function yNorm = amplitudeNormalize(y)

    yMin = min(y, [], "omitnan");
    yMax = max(y, [], "omitnan");

    if yMax > yMin
        yNorm = (y - yMin) ./ (yMax - yMin);
    else
        yNorm = nan(size(y));
    end
end
end