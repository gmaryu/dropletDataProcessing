function classifyCdk1dynamics(nucData, dataSet, traceInfo, resultDataPath, resultDataSavePath)
%% Nuclear-reconstituted droplet Cdk1 cycle classification
%
% Classification goal:
%   State 1: cytoplasmic-like Cdk1 activity pattern
%   State 2: high or gradually increasing interphase nuclear Cdk1 activity
%   State 3: low interphase nuclear Cdk1 activity with late G2/M-like activation
%            plus high N/C volume ratio or within-droplet long-cycle context
%
% Key assumptions:
%   - Cytoplasm-only reference was built as peak-to-peak cycles.
%   - Nuclear cycle table defines each cycle as Peak1Frame -> Peak2Frame.
%   - TroughFrame is the minimum between Peak1Frame and Peak2Frame.
%   - Droplets are uniquely identified by Position + Droplet/Track ID.
%   - traceInfo maps each DropletID to the row index of the time-series matrices.
%
% Expected variables in workspace:
%   resultDataPath
%   nucData
%   dataSet
%   traceInfo
%
% Main outputs:
%   nuclear_reconstituted_classification_peak2peak_refactored.mat
%   existing_cycle_table_with_states_peak2peak_refactored.csv

% nucData = nucData_g5;
% dataSet = g5;
% traceInfo = traceInfo_g5;

clc; close all;

mkdir(resultDataSavePath);

P = initParameters();

[cytoRef, phase, nPhasePoints] = loadCytoplasmicReference(resultDataPath);
raw = getNuclearRawMatrices(nucData);
validateRawMatrices(raw);

[nDroplets, nTime] = size(raw.total);
fprintf("Loaded %d nuclear-reconstituted droplets and %d time points.\n", nDroplets, nTime);

[existingCycleTable, cycleBoundaryTable, traceInfo] = prepareCycleBoundaryTable(dataSet, traceInfo, nDroplets);
featuresRaw = smoothAndDeriveSignals(raw, P);

[nucCycleInfo, cycles] = segmentPeakToPeakCycles(cycleBoundaryTable, featuresRaw, nDroplets, nTime);
nucCycleInfo = applyCycleQC(nucCycleInfo, P);
plotCycleQC(nucCycleInfo);

[normCycles, shapeCycles] = normalizeCycleSignals(cycles, phase, nPhasePoints);

[nucCycleInfo, featureContext] = extractCycleFeatures( ...
    nucCycleInfo, normCycles, shapeCycles, featuresRaw, cytoRef, phase, nPhasePoints, P);

nucCycleInfo = addWithinDropletCycleLengthContext(nucCycleInfo, P);
thresholds = computeEffectiveThresholds(nucCycleInfo, P);

nucCycleInfo = classifyStates(nucCycleInfo, thresholds, P);
nucCycleInfo = applyOptionalStateSmoothing(nucCycleInfo, P);

transitionTable = buildTransitionTable(nucCycleInfo);
existingCycleTable_withStates = addResultsBackToOriginalCycleTable(existingCycleTable, nucCycleInfo);

makeDiagnosticPlots(nucCycleInfo, cytoRef, shapeCycles, phase, nPhasePoints, thresholds, P);
saveClassificationOutputs( ...
    nucCycleInfo, existingCycleTable_withStates, transitionTable, cycles, normCycles, shapeCycles, ...
    phase, nPhasePoints, P, thresholds, featureContext);

fprintf("Refactored classification completed.\n");

%% =========================
%  Local functions: setup
%  =========================

function P = initParameters()

    %% Smoothing
    P.detectSmoothWindow = 1;   % retained for compatibility; no independent peak detection is done
    P.featureSmoothWindow = 1;  % 1 = no smoothing for feature extraction

    %% Cycle QC
    P.minCycleLength = 5;
    P.maxCycleLength = inf;
    P.minAmplitude = 0.05;
    P.minTroughPhase = 0.05;
    P.maxTroughPhase = 0.95;
    P.maxNaNFractionTotal = 0.20;

    %% Peak-to-peak trough-based windows
    P.interphaseOffsetAfterTrough = 35;
    P.lateOffsetAfterTrough = 35; % retained for compatibility/documentation
    P.minWindowPoints = 3;

    %% Interphase end detection from NLS-defined nuclear volume
    P.useNucVolForInterphaseEnd = true;
    P.nucVolLossThreshold = 0;
    P.useRelativeNucVolLossThreshold = false;
    P.relativeNucVolLossThreshold = 0.05;
    P.nucVolLossConsecutiveN = 2;
    P.minFramesAfterTroughForLossSearch = 1;
    P.excludeLastNPhasePointsFromInterphase = 3;

    %% State 1: cytoplasmic-like waveform
    P.useCorrForState1 = true;
    P.useRMSEForState1 = true;
    P.state1UseOrLogic = true;     % true: corr OR RMSE, false: corr AND RMSE
    P.state1CorrMultiplier = 1.5;
    P.state1RMSEMultiplier = 1.5;
    P.useMedianCytoTemplateForDistance = false;

    %% State 2: high or gradually increasing interphase nuclear Cdk1
    P.nucRatioHighThreshold = 1.10;
    P.useDataDrivenNucRatioThreshold = false;
    P.useP90ForState2 = true;
    P.useLateInterphaseMeanForState2 = true;
    P.useIncreaseForState2 = true;
    P.useNucMinusCytoForState2 = false;
    P.nucRatioIncreaseThreshold = 0.05;
    P.nucMinusCytoHighThreshold = 0.00;
    P.interphasePercentileForState2 = 90; % historical name kept as P90 output column for compatibility

    %% State 3: late nuclear activation + context gate
    % State 3 is defined as:
    %   late nuclear activation
    %   AND (high N/C volume ratio OR within-droplet long cycle)
    % A global/cytoplasm-reference long-cycle gate is intentionally not used.
    P.lateActivationThreshold = 0.25;
    P.lateTimingOffsetFromTrough = 0.20;  % phase units
    P.useState3ContextGate = true;
    P.useDataDrivenState3ContextThresholds = true;
    P.highNCvolRatioThreshold = 0.05;      % if NaN and data-driven fails, N/C criterion is ignored
    P.state3ContextMADMultiplier = 3;
    P.useInterphaseNCvolRatioForState3 = true;

    %% Within-droplet long-cycle context for State 3
    P.withinDropletFoldChangeThreshold = 1.50;
    P.withinDropletDeltaThreshold = 20;
    P.withinDropletZThreshold = 5.0;

    %% Optional state-sequence post-processing
    P.applyConsecutiveSmoothing = false;
    P.consecutiveN = 2;

    %% Output
    P.outputMatFile = "nuclear_reconstituted_classification_peak2peak_refactored.mat";
    P.outputCsvFile = "existing_cycle_table_with_states_peak2peak_refactored.csv";
end

%%
figure("Name", "State 3 context gate",'Position',[0,0,300,300],'Visible','off');
scatter( ...
nucCycleInfo.InterphaseNCvolRatioMean, ...
nucCycleInfo.CycleLength, ...
35, ...
nucCycleInfo.State_raw, ...
"filled");
hold on;
xscale('log')
xlabel('N/C Ratio');
ylabel('Period length');
xline(0.012,"k--",'in vivo Sync', "LineWidth", 1.5,'Color','blue');
xline(0.025, "k--",'in vivo MBT', "LineWidth", 1.5,'Color','red');
xlim([10^-5,10^0]);
ylim([0,60]);
axis('square')
exportgraphics(gcf, fullfile(resultDataSavePath, 'All_NCVolRatio_vs_Periods.png'));
%%

function [cytoRef, phase, nPhasePoints] = loadCytoplasmicReference(resultDataPath)

    load(fullfile(resultDataPath, "cytoplasm_only_reference_peak2peak.mat"), ...
        "cytoTemplateByCycle_peak2peak", ...
        "cytoTemplateMedianByCycle_peak2peak", ...
        "cytoThresholdCorrByCycle_peak2peak", ...
        "cytoThresholdRMSEByCycle_peak2peak", ...
        "cytoReferenceSummary_peak2peak", ...
        "phase", ...
        "nPhasePoints");

    cytoRef.templateMean = cytoTemplateByCycle_peak2peak;
    cytoRef.templateMedian = cytoTemplateMedianByCycle_peak2peak;
    cytoRef.thresholdCorr = cytoThresholdCorrByCycle_peak2peak;
    cytoRef.thresholdRMSE = cytoThresholdRMSEByCycle_peak2peak;
    cytoRef.summary = cytoReferenceSummary_peak2peak;

    fprintf("Loaded peak-to-peak cytoplasmic reference.\n");
end

function raw = getNuclearRawMatrices(nucData)

    raw.total = nucData.WholeMean_RATIO;
    raw.nuc = nucData.nucMean_RATIO;
    raw.cyto = nucData.cytMean_RATIO;
    raw.nucVol = nucData.nucVol;
    raw.dropletVol = nucData.dropletVol;
    raw.NCvolRatio = nucData.NCVolRatio;
end

function validateRawMatrices(raw)

    baseSize = size(raw.total);
    fields = fieldnames(raw);

    for i = 1:numel(fields)
        f = fields{i};
        assert(isequal(baseSize, size(raw.(f))), ...
            "%s must have the same size as raw.total.", f);
    end
end

function [existingCycleTable, cycleBoundaryTable, traceInfo] = prepareCycleBoundaryTable(dataSet, traceInfo, nDroplets)

    existingCycleTable = dataSet.cycle;
    existingCycleTable.OriginalCycleRow = (1:height(existingCycleTable))';
    existingCycleTable.DropletID = makeDropletID(existingCycleTable.POS_ID, existingCycleTable.TRACK_ID);

    if ~ismember("DropletID", traceInfo.Properties.VariableNames)
        if all(ismember(["PositionNumber", "DropletNumber"], string(traceInfo.Properties.VariableNames)))
            traceInfo.DropletID = makeDropletID(traceInfo.PositionNumber, traceInfo.DropletNumber);
        else
            error("traceInfo must contain DropletID or PositionNumber/DropletNumber.");
        end
    end

    if ~ismember("TraceRow", traceInfo.Properties.VariableNames)
        error("traceInfo must contain TraceRow, the row index in the time-series matrices.");
    end

    checkTraceInfo(traceInfo, nDroplets);

    cycleBoundaryTable = table;
    cycleBoundaryTable.OriginalCycleRow = existingCycleTable.OriginalCycleRow;
    cycleBoundaryTable.Position = existingCycleTable.POS_ID;
    cycleBoundaryTable.Droplet = existingCycleTable.TRACK_ID;
    cycleBoundaryTable.CycleInDroplet = existingCycleTable.CYCLE_ID;
    cycleBoundaryTable.DropletID = existingCycleTable.DropletID;
    cycleBoundaryTable.Peak1Frame = existingCycleTable.START_INDEX;
    cycleBoundaryTable.Peak2Frame = existingCycleTable.END_INDEX;
    cycleBoundaryTable.TroughFrame = existingCycleTable.TROUGH_INDEX;

    optionalColumns = ["Peak1Value", "Peak2Value", "TroughValue"];
    for i = 1:numel(optionalColumns)
        col = optionalColumns(i);
        if ismember(col, string(existingCycleTable.Properties.VariableNames))
            cycleBoundaryTable.(col) = existingCycleTable.(col);
        end
    end

    cycleBoundaryTable = outerjoin( ...
        cycleBoundaryTable, ...
        traceInfo(:, ["DropletID", "TraceRow"]), ...
        "Keys", "DropletID", ...
        "MergeKeys", true);

    cycleBoundaryTable = sortrows(cycleBoundaryTable, "OriginalCycleRow");

    missingTrace = isnan(cycleBoundaryTable.TraceRow);
    fprintf("Cycles without matched TraceRow: %d / %d\n", ...
        sum(missingTrace), height(cycleBoundaryTable));

    if any(missingTrace)
        disp(cycleBoundaryTable(missingTrace, ...
            ["Position", "Droplet", "DropletID", "CycleInDroplet"]));
    end
end

function dropletID = makeDropletID(positionID, trackID)
    dropletID = strcat("P", string(positionID), "_D", string(trackID));
end

function checkTraceInfo(traceInfo, nDroplets)

    Gtrace = groupcounts(traceInfo, "DropletID");
    if any(Gtrace.GroupCount > 1)
        warning("Some DropletIDs map to multiple TraceRows. Inspect traceInfo before trusting classification.");
        disp(Gtrace(Gtrace.GroupCount > 1, :));
    end

    badRows = traceInfo.TraceRow < 1 | traceInfo.TraceRow > nDroplets | isnan(traceInfo.TraceRow);
    if any(badRows)
        warning("Some traceInfo.TraceRow values are outside the time-series matrix row range.");
        disp(traceInfo(badRows, :));
    end
end

%% =========================
%  Local functions: preprocessing
%  =========================

function featuresRaw = smoothAndDeriveSignals(raw, P)

    featuresRaw.totalDetect = smoothTraceMatrixNoBridge(raw.total, P.detectSmoothWindow);

    featuresRaw.total = smoothTraceMatrixNoBridge(raw.total, P.featureSmoothWindow);
    featuresRaw.nuc = smoothTraceMatrixNoBridge(raw.nuc, P.featureSmoothWindow);
    featuresRaw.cyto = smoothTraceMatrixNoBridge(raw.cyto, P.featureSmoothWindow);
    featuresRaw.nucVol = smoothTraceMatrixNoBridge(raw.nucVol, P.featureSmoothWindow);
    featuresRaw.dropletVol = smoothTraceMatrixNoBridge(raw.dropletVol, P.featureSmoothWindow);
    featuresRaw.NCvolRatio = smoothTraceMatrixNoBridge(raw.NCvolRatio, P.featureSmoothWindow);

    epsilon = 1e-6;
    featuresRaw.nucCytoRatio = featuresRaw.nuc ./ (featuresRaw.cyto + epsilon);
    featuresRaw.nucMinusCyto = featuresRaw.nuc - featuresRaw.cyto;
end

function [nucCycleInfo, cycles] = segmentPeakToPeakCycles(cycleBoundaryTable, X, nDroplets, nTime)

    nCyclesTotal = height(cycleBoundaryTable);
    cycles = initializeCycleCells(nCyclesTotal);
    nucCycleInfo = cycleBoundaryTable;

    meta = initializeSegmentationMeta(nCyclesTotal);

    for c = 1:nCyclesTotal
        row = nucCycleInfo.TraceRow(c);
        idx1 = nucCycleInfo.Peak1Frame(c);
        idx2 = nucCycleInfo.Peak2Frame(c);
        troughIdx = nucCycleInfo.TroughFrame(c);

        if ~isUsableBoundary(row, idx1, idx2, nDroplets, nTime)
            continue;
        end

        row = round(row);
        idx1 = round(idx1);
        idx2 = round(idx2);
        frameIdx = idx1:idx2;

        cycles.total{c} = X.total(row, frameIdx);
        cycles.nuc{c} = X.nuc(row, frameIdx);
        cycles.cyto{c} = X.cyto(row, frameIdx);
        cycles.nucVol{c} = X.nucVol(row, frameIdx);
        cycles.dropletVol{c} = X.dropletVol(row, frameIdx);
        cycles.NCvolRatio{c} = X.NCvolRatio(row, frameIdx);
        cycles.nucCytoRatio{c} = X.nucCytoRatio(row, frameIdx);
        cycles.nucMinusCyto{c} = X.nucMinusCyto(row, frameIdx);

        meta = updateSegmentationMeta(meta, c, row, idx1, idx2, troughIdx, X.total, frameIdx);
        meta.NaNFractionNucVol(c) = mean(isnan(cycles.nucVol{c}));
        meta.NaNFractionNCvolRatio(c) = mean(isnan(cycles.NCvolRatio{c}));
    end

    nucCycleInfo = addStructFieldsToTable(nucCycleInfo, meta);
end

function cycles = initializeCycleCells(n)
    names = ["total", "nuc", "cyto", "nucVol", "dropletVol", "NCvolRatio", "nucCytoRatio", "nucMinusCyto"];
    for i = 1:numel(names)
        cycles.(names(i)) = cell(n, 1);
    end
end

function meta = initializeSegmentationMeta(n)
    meta.CycleLength = nan(n, 1);
    meta.Peak1Frame = nan(n, 1);
    meta.Peak2Frame = nan(n, 1);
    meta.TroughFrame = nan(n, 1);
    meta.TroughFrameRelative = nan(n, 1);
    meta.TroughPhase = nan(n, 1);
    meta.TotalPeak1Value = nan(n, 1);
    meta.TotalPeak2Value = nan(n, 1);
    meta.TotalTroughValue = nan(n, 1);
    meta.TotalAmplitudeFromTroughToPeak2 = nan(n, 1);
    meta.TotalAmplitudeMaxMinusMin = nan(n, 1);
    meta.NaNFractionTotal = nan(n, 1);
    meta.NaNFractionNucVol = nan(n, 1);
    meta.NaNFractionNCvolRatio = nan(n, 1);
end

function tf = isUsableBoundary(row, idx1, idx2, nDroplets, nTime)
    tf = ~(isnan(row) || isnan(idx1) || isnan(idx2) || idx2 <= idx1);
    if ~tf
        return;
    end
    row = round(row);
    idx1 = round(idx1);
    idx2 = round(idx2);
    tf = row >= 1 && row <= nDroplets && idx1 >= 1 && idx2 <= nTime;
end

function meta = updateSegmentationMeta(meta, c, row, idx1, idx2, troughIdx, Xtotal, frameIdx)

    thisTotal = Xtotal(row, frameIdx);

    meta.CycleLength(c) = numel(frameIdx);
    meta.Peak1Frame(c) = idx1;
    meta.Peak2Frame(c) = idx2;
    meta.TotalPeak1Value(c) = Xtotal(row, idx1);
    meta.TotalPeak2Value(c) = Xtotal(row, idx2);

    if ~isnan(troughIdx)
        troughIdx = round(troughIdx);
        meta.TroughFrame(c) = troughIdx;

        if troughIdx >= idx1 && troughIdx <= idx2
            troughRel = troughIdx - idx1 + 1;
            meta.TroughFrameRelative(c) = troughRel;
            meta.TroughPhase(c) = troughRel / numel(frameIdx);
            meta.TotalTroughValue(c) = Xtotal(row, troughIdx);
        end
    end

    if ~isnan(meta.TotalTroughValue(c))
        meta.TotalAmplitudeFromTroughToPeak2(c) = meta.TotalPeak2Value(c) - meta.TotalTroughValue(c);
    end

    meta.TotalAmplitudeMaxMinusMin(c) = max(thisTotal, [], "omitnan") - min(thisTotal, [], "omitnan");
    meta.NaNFractionTotal(c) = mean(isnan(thisTotal));
end

function nucCycleInfo = applyCycleQC(nucCycleInfo, P)

    isValidCycle = ...
        nucCycleInfo.CycleLength >= P.minCycleLength & ...
        nucCycleInfo.CycleLength <= P.maxCycleLength & ...
        nucCycleInfo.TotalAmplitudeMaxMinusMin >= P.minAmplitude & ...
        nucCycleInfo.TroughPhase >= P.minTroughPhase & ...
        nucCycleInfo.TroughPhase <= P.maxTroughPhase & ...
        nucCycleInfo.NaNFractionTotal < P.maxNaNFractionTotal;

    nucCycleInfo.IsValidCycle = isValidCycle;

    fprintf("Valid cycles: %d / %d\n", sum(isValidCycle), height(nucCycleInfo));
end

function [normCycles, shapeCycles] = normalizeCycleSignals(cycles, phase, nPhasePoints)

    nCyclesTotal = numel(cycles.total);
    names = fieldnames(cycles);

    for i = 1:numel(names)
        name = names{i};
        normCycles.(name) = nan(nCyclesTotal, nPhasePoints);
        for c = 1:nCyclesTotal
            normCycles.(name)(c, :) = phaseResampleAllowNaN(cycles.(name){c}, phase);
        end
    end

    shapeCycles.total = nan(nCyclesTotal, nPhasePoints);
    shapeCycles.nuc = nan(nCyclesTotal, nPhasePoints);
    shapeCycles.cyto = nan(nCyclesTotal, nPhasePoints);

    for c = 1:nCyclesTotal
        shapeCycles.total(c, :) = amplitudeNormalize(normCycles.total(c, :));
        shapeCycles.nuc(c, :) = amplitudeNormalize(normCycles.nuc(c, :));
        shapeCycles.cyto(c, :) = amplitudeNormalize(normCycles.cyto(c, :));
    end
end

%% =========================
%  Local functions: feature extraction
%  =========================

function [nucCycleInfo, ctx] = extractCycleFeatures(nucCycleInfo, normCycles, shapeCycles, X, cytoRef, phase, nPhasePoints, P)

    n = height(nucCycleInfo);
    F = initializeFeatureStruct(n);

    for c = 1:n
        F = computeCytoplasmicDistance(F, c, nucCycleInfo.CycleInDroplet(c), shapeCycles.total(c, :), cytoRef, P);

        [windows, F] = defineCycleWindows(F, c, nucCycleInfo, X.nucVol, nPhasePoints, P);
        if windows.skip
            continue;
        end

        F = computeCdk1WindowFeatures(F, c, normCycles, shapeCycles, phase, windows, P);
        F = computeVolumeWindowFeatures(F, c, normCycles, windows);
    end

    nucCycleInfo = addStructFieldsToTable(nucCycleInfo, F);

    ctx.interphaseWindow = [];
    ctx.lateWindow = [];
end

function F = initializeFeatureStruct(n)

    names = [ ...
        "CytoDistanceCorr", "CytoDistanceRMSE", "CytoThresholdCorr", "CytoThresholdRMSE", ...
        "TroughIdxNorm", "InterphaseStartIdx", "InterphaseEndIdx", "LateStartIdx", "LateEndIdx", ...
        "InterphaseStartFrame", "InterphaseEndFrame", ...
        "InterphaseNucMean", "InterphaseCytoMean", "InterphaseNucCytoRatioMean", ...
        "InterphaseNucCytoRatioMax", "InterphaseNucCytoRatioP90", ...
        "EarlyInterphaseNucCytoRatioMean", "LateInterphaseNucCytoRatioMean", ...
        "InterphaseNucCytoRatioIncrease", "InterphaseNucCytoRatioSlope", ...
        "InterphaseNucMinusCytoMean", "LateNucMean", "LateCytoMean", ...
        "LateActivationScore_nucShape", "TimeToHalfActivation_nuc", "PeakPhase_nuc", "MaxSlopeLate_nucShape", ...
        "InterphaseNucVolMean", "InterphaseDropletVolMean", "InterphaseNCvolRatioMean", ...
        "CycleMeanNucVol", "CycleMeanNCvolRatio", ...
        "NucVolValidFraction", "NCvolRatioValidFraction", ...
        "InterphaseNucVolValidFraction", "InterphaseNCvolRatioValidFraction"];

    for i = 1:numel(names)
        F.(names(i)) = nan(n, 1);
    end

    F.InterphaseEndReason = strings(n, 1);
end

function F = computeCytoplasmicDistance(F, c, k, y, cytoRef, P)

    if isnan(k) || k < 1 || k > size(cytoRef.templateMean, 1)
        return;
    end

    if P.useMedianCytoTemplateForDistance
        template = cytoRef.templateMedian(k, :);
    else
        template = cytoRef.templateMean(k, :);
    end

    if all(isnan(template))
        return;
    end

    r = corr(y(:), template(:), "Rows", "complete");
    F.CytoDistanceCorr(c) = 1 - r; % close to 0, if the dynamics looks like cytoplasmic activation
    F.CytoDistanceRMSE(c) = sqrt(mean((y - template).^2, "omitnan"));
    F.CytoThresholdCorr(c) = cytoRef.thresholdCorr(k);
    F.CytoThresholdRMSE(c) = cytoRef.thresholdRMSE(k);
end

function [windows, F] = defineCycleWindows(F, c, nucCycleInfo, XnucVol, nPhasePoints, P)

    windows = struct("skip", true);

    idx1 = nucCycleInfo.Peak1Frame(c);
    idx2 = nucCycleInfo.Peak2Frame(c);
    troughFrame = nucCycleInfo.TroughFrame(c);
    traceRow = nucCycleInfo.TraceRow(c);

    if isnan(idx1) || isnan(idx2) || isnan(troughFrame) || isnan(traceRow)
        return;
    end

    idx1 = round(idx1);
    idx2 = min(size(XnucVol, 2), round(idx2));
    troughFrame = max(idx1, min(idx2, round(troughFrame)));
    traceRow = round(traceRow);

    troughIdxNorm = frameToPhaseIndex(troughFrame, idx1, idx2, nPhasePoints);
    if isnan(troughIdxNorm)
        return;
    end

    interphaseEndFrame = detectInterphaseEndFrameFromNucVol(traceRow, idx1, idx2, troughFrame, troughIdxNorm, XnucVol, nPhasePoints, P);
    interphaseStart = frameToPhaseIndex(troughFrame, idx1, idx2, nPhasePoints);
    interphaseEnd = frameToPhaseIndex(interphaseEndFrame, idx1, idx2, nPhasePoints);

    interphaseStart = max(1, min(nPhasePoints, interphaseStart));
    interphaseEnd = max(interphaseStart, min(nPhasePoints, interphaseEnd));
    interphaseEnd = min(interphaseEnd, nPhasePoints - P.excludeLastNPhasePointsFromInterphase);
    interphaseEnd = max(interphaseStart, interphaseEnd);

    if interphaseEnd - interphaseStart + 1 < P.minWindowPoints
        interphaseStart = max(1, interphaseEnd - P.minWindowPoints + 1);
    end

    lateStart = min(interphaseEnd + 1, nPhasePoints);
    lateEnd = nPhasePoints;
    if lateEnd - lateStart + 1 < P.minWindowPoints
        lateStart = max(1, nPhasePoints - P.minWindowPoints + 1);
    end

    windows.skip = false;
    windows.troughIdxNorm = troughIdxNorm;
    windows.interphase = interphaseStart:interphaseEnd;
    windows.late = lateStart:lateEnd;

    F.TroughIdxNorm(c) = troughIdxNorm;
    F.InterphaseStartIdx(c) = interphaseStart;
    F.InterphaseEndIdx(c) = interphaseEnd;
    F.LateStartIdx(c) = lateStart;
    F.LateEndIdx(c) = lateEnd;
    F.InterphaseStartFrame(c) = troughFrame;
    F.InterphaseEndFrame(c) = interphaseEndFrame;

    if P.useNucVolForInterphaseEnd
        F.InterphaseEndReason(c) = getInterphaseEndReason(traceRow, idx1, idx2, troughFrame, XnucVol, P);
    else
        F.InterphaseEndReason(c) = "fixedOffset";
    end
end

function endFrame = detectInterphaseEndFrameFromNucVol(traceRow, idx1, idx2, troughFrame, troughIdxNorm, XnucVol, nPhasePoints, P)

    fallbackEndIdx = min(troughIdxNorm + P.interphaseOffsetAfterTrough, ...
        nPhasePoints - P.excludeLastNPhasePointsFromInterphase);
    fallbackEndIdx = max(troughIdxNorm, fallbackEndIdx);
    fallbackEndFrame = phaseIndexToFrame(fallbackEndIdx, idx1, idx2, nPhasePoints);

    if ~P.useNucVolForInterphaseEnd
        endFrame = fallbackEndFrame;
        return;
    end

    nucVolCycleRaw = XnucVol(traceRow, idx1:idx2);
    frameIdxCycle = idx1:idx2;
    searchStartFrame = max(idx1, min(idx2, troughFrame + P.minFramesAfterTroughForLossSearch));
    searchMask = frameIdxCycle >= searchStartFrame;

    nucVolSearch = nucVolCycleRaw(searchMask);
    frameSearch = frameIdxCycle(searchMask);

    if P.useRelativeNucVolLossThreshold
        nucVolCycleMax = max(nucVolCycleRaw, [], "omitnan");
        nucVolLossThresholdThisCycle = P.relativeNucVolLossThreshold * nucVolCycleMax;
    else
        nucVolLossThresholdThisCycle = P.nucVolLossThreshold;
    end

    lostNucleus = isnan(nucVolSearch) | nucVolSearch <= nucVolLossThresholdThisCycle;
    lossLocalIdx = firstConsecutiveTrue(lostNucleus, P.nucVolLossConsecutiveN);

    if ~isnan(lossLocalIdx)
        firstLossFrame = frameSearch(lossLocalIdx);
        endFrame = max(troughFrame, min(idx2, firstLossFrame - 1));
    else
        endFrame = fallbackEndFrame;
    end
end

function reason = getInterphaseEndReason(traceRow, idx1, idx2, troughFrame, XnucVol, P)

    nucVolCycleRaw = XnucVol(traceRow, idx1:idx2);
    frameIdxCycle = idx1:idx2;
    searchStartFrame = max(idx1, min(idx2, troughFrame + P.minFramesAfterTroughForLossSearch));
    searchMask = frameIdxCycle >= searchStartFrame;
    nucVolSearch = nucVolCycleRaw(searchMask);

    if P.useRelativeNucVolLossThreshold
        nucVolCycleMax = max(nucVolCycleRaw, [], "omitnan");
        threshold = P.relativeNucVolLossThreshold * nucVolCycleMax;
    else
        threshold = P.nucVolLossThreshold;
    end

    lostNucleus = isnan(nucVolSearch) | nucVolSearch <= threshold;
    if ~isnan(firstConsecutiveTrue(lostNucleus, P.nucVolLossConsecutiveN))
        reason = "nucVolLoss";
    else
        reason = "fallbackFixedOffset";
    end
end

function F = computeCdk1WindowFeatures(F, c, normCycles, shapeCycles, phase, windows, P)

    iw = windows.interphase;
    lw = windows.late;

    F.InterphaseNucMean(c) = mean(normCycles.nuc(c, iw), "omitnan");
    F.InterphaseCytoMean(c) = mean(normCycles.cyto(c, iw), "omitnan");
    F.InterphaseNucCytoRatioMean(c) = mean(normCycles.nucCytoRatio(c, iw), "omitnan");
    F.InterphaseNucMinusCytoMean(c) = mean(normCycles.nucMinusCyto(c, iw), "omitnan");

    ratioTrace = normCycles.nucCytoRatio(c, iw);
    F.InterphaseNucCytoRatioMax(c) = max(ratioTrace, [], "omitnan");
    F.InterphaseNucCytoRatioP90(c) = prctileOmitNaN(ratioTrace, P.interphasePercentileForState2);

    [earlyW, lateInterphaseW] = splitWindowInHalf(iw);
    F.EarlyInterphaseNucCytoRatioMean(c) = mean(normCycles.nucCytoRatio(c, earlyW), "omitnan");
    F.LateInterphaseNucCytoRatioMean(c) = mean(normCycles.nucCytoRatio(c, lateInterphaseW), "omitnan");
    F.InterphaseNucCytoRatioIncrease(c) = ...
        F.LateInterphaseNucCytoRatioMean(c) - F.EarlyInterphaseNucCytoRatioMean(c);

    F.InterphaseNucCytoRatioSlope(c) = computeLinearSlope(phase(iw)', ratioTrace(:));

    F.LateNucMean(c) = mean(normCycles.nuc(c, lw), "omitnan");
    F.LateCytoMean(c) = mean(normCycles.cyto(c, lw), "omitnan");
    F.LateActivationScore_nucShape(c) = ...
        mean(shapeCycles.nuc(c, lw), "omitnan") - mean(shapeCycles.nuc(c, iw), "omitnan");

    yNucShape = shapeCycles.nuc(c, :);
    [F.TimeToHalfActivation_nuc(c), F.PeakPhase_nuc(c), F.MaxSlopeLate_nucShape(c)] = ...
        computeLateActivationTiming(yNucShape, windows.troughIdxNorm, lw, numel(yNucShape));
end

function [earlyW, lateW] = splitWindowInHalf(w)
    nW = numel(w);
    splitPoint = max(1, floor(nW / 2));
    earlyW = w(1:splitPoint);
    lateW = w(min(splitPoint + 1, nW):end);
end

function slope = computeLinearSlope(x, y)
    valid = ~isnan(x) & ~isnan(y);
    if sum(valid) >= 3
        p = polyfit(x(valid), y(valid), 1);
        slope = p(1);
    else
        slope = NaN;
    end
end

function [timeToHalf, peakPhase, maxSlopeLate] = computeLateActivationTiming(y, searchStart, lateWindow, nPhasePoints)

    timeToHalf = NaN;
    peakPhase = NaN;
    maxSlopeLate = NaN;

    if all(isnan(y))
        return;
    end

    [~, peakIdx] = max(y);
    peakPhase = peakIdx / nPhasePoints;

    halfMax = min(y, [], "omitnan") + 0.5 * (max(y, [], "omitnan") - min(y, [], "omitnan"));
    idxHalf = find(y(searchStart:end) >= halfMax, 1, "first");

    if ~isempty(idxHalf)
        timeToHalf = (searchStart + idxHalf - 1) / nPhasePoints;
    end

    if numel(lateWindow) >= 2
        dyLate = diff(y(lateWindow));
        maxSlopeLate = max(dyLate, [], "omitnan");
    end
end

function F = computeVolumeWindowFeatures(F, c, normCycles, windows)

    iw = windows.interphase;

    F.InterphaseNucVolMean(c) = mean(normCycles.nucVol(c, iw), "omitnan");
    F.InterphaseDropletVolMean(c) = mean(normCycles.dropletVol(c, iw), "omitnan");
    F.InterphaseNCvolRatioMean(c) = mean(normCycles.NCvolRatio(c, iw), "omitnan");

    F.CycleMeanNucVol(c) = mean(normCycles.nucVol(c, :), "omitnan");
    F.CycleMeanNCvolRatio(c) = mean(normCycles.NCvolRatio(c, :), "omitnan");

    F.NucVolValidFraction(c) = mean(~isnan(normCycles.nucVol(c, :)));
    F.NCvolRatioValidFraction(c) = mean(~isnan(normCycles.NCvolRatio(c, :)));
    F.InterphaseNucVolValidFraction(c) = mean(~isnan(normCycles.nucVol(c, iw)));
    F.InterphaseNCvolRatioValidFraction(c) = mean(~isnan(normCycles.NCvolRatio(c, iw)));
end

%% =========================
%  Local functions: classification
%  =========================

function nucCycleInfo = addWithinDropletCycleLengthContext(nucCycleInfo, P)

    n = height(nucCycleInfo);
    CycleLengthFoldChangePrev = nan(n, 1);
    CycleLengthDeltaPrev = nan(n, 1);
    CycleLengthZWithinDroplet = nan(n, 1);
    IsWithinDropletLongCycle = false(n, 1);

    uniqueDropletIDs = unique(string(nucCycleInfo.DropletID));

    for i = 1:numel(uniqueDropletIDs)
        thisID = uniqueDropletIDs(i);
        idx = find(string(nucCycleInfo.DropletID) == thisID & nucCycleInfo.IsValidCycle);

        if numel(idx) < 2
            continue;
        end

        [~, order] = sort(nucCycleInfo.CycleInDroplet(idx));
        idx = idx(order);

        L = nucCycleInfo.CycleLength(idx);
        prevL = [NaN; L(1:end-1)];

        CycleLengthFoldChangePrev(idx) = L ./ prevL;
        CycleLengthDeltaPrev(idx) = L - prevL;

        medL = median(L, "omitnan");
        madL = mad(L, 1);

        if madL > 0
            CycleLengthZWithinDroplet(idx) = (L - medL) ./ (1.4826 * madL);
        end

        IsWithinDropletLongCycle(idx) = ...
            CycleLengthFoldChangePrev(idx) >= P.withinDropletFoldChangeThreshold | ...
            CycleLengthDeltaPrev(idx) >= P.withinDropletDeltaThreshold | ...
            CycleLengthZWithinDroplet(idx) >= P.withinDropletZThreshold;
    end

    nucCycleInfo.CycleLengthFoldChangePrev = CycleLengthFoldChangePrev;
    nucCycleInfo.CycleLengthDeltaPrev = CycleLengthDeltaPrev;
    nucCycleInfo.CycleLengthZWithinDroplet = CycleLengthZWithinDroplet;
    nucCycleInfo.IsWithinDropletLongCycle = IsWithinDropletLongCycle;
end

function thresholds = computeEffectiveThresholds(nucCycleInfo, P)

    idxState1Candidate = nucCycleInfo.IsValidCycle & ...
        (nucCycleInfo.CytoDistanceCorr <= P.state1CorrMultiplier * nucCycleInfo.CytoThresholdCorr | ...
         nucCycleInfo.CytoDistanceRMSE <= P.state1RMSEMultiplier * nucCycleInfo.CytoThresholdRMSE);

    thresholds.effectiveNucRatioHighThreshold = P.nucRatioHighThreshold;
    if P.useDataDrivenNucRatioThreshold
        baseRatio = nucCycleInfo.InterphaseNucCytoRatioP90(idxState1Candidate);
        if sum(~isnan(baseRatio)) >= 5
            thresholds.effectiveNucRatioHighThreshold = median(baseRatio, "omitnan") + 3 * mad(baseRatio, 1);
        else
            warning("Not enough State-1-like cycles for data-driven ratio threshold. Using manual threshold.");
        end
    end

    thresholds.highNCvolRatioThreshold = P.highNCvolRatioThreshold;
    if P.useDataDrivenState3ContextThresholds
        if P.useInterphaseNCvolRatioForState3
            baseNC = nucCycleInfo.InterphaseNCvolRatioMean(idxState1Candidate);
        else
            baseNC = nucCycleInfo.CycleMeanNCvolRatio(idxState1Candidate);
        end

        if sum(~isnan(baseNC)) >= 5
            thresholds.highNCvolRatioThreshold = median(baseNC, "omitnan") + ...
                P.state3ContextMADMultiplier * mad(baseNC, 1);
        else
            warning("Not enough State-1-like cycles for data-driven N/C threshold. Keeping manual value.");
        end
    end

    fprintf("State 2 nuclear/cytoplasmic ratio threshold = %.4f\n", ...
        thresholds.effectiveNucRatioHighThreshold);
    fprintf("State 3 context thresholds:\n");
    fprintf("  highNCvolRatioThreshold = %.4f\n", thresholds.highNCvolRatioThreshold);
    fprintf("  withinDropletFoldChangeThreshold = %.2f\n", P.withinDropletFoldChangeThreshold);
    fprintf("  withinDropletDeltaThreshold = %.2f frames\n", P.withinDropletDeltaThreshold);
    fprintf("  withinDropletZThreshold = %.2f\n", P.withinDropletZThreshold);
end

function nucCycleInfo = classifyStates(nucCycleInfo, thresholds, P)

    n = height(nucCycleInfo);

    State_raw = nan(n, 1);
    IsCytoLike = false(n, 1);
    IsHighInterphaseNuc = false(n, 1);
    IsLateG2Mlike = false(n, 1);
    IsHighNCvolRatio = false(n, 1);
    IsWithinDropletLongCycleForState3 = false(n, 1);
    HasState3Context = false(n, 1);

    for c = 1:n
        if ~nucCycleInfo.IsValidCycle(c)
            continue;
        end

        isCytoLike = classifyState1(nucCycleInfo, c, P);
        isHighInterphaseNuc = classifyState2(nucCycleInfo, c, thresholds, P);
        [isLateG2Mlike, isHighNCvolRatio, isWithinLong, hasContext] = classifyState3(nucCycleInfo, c, thresholds, P);

        IsCytoLike(c) = isCytoLike;
        IsHighInterphaseNuc(c) = isHighInterphaseNuc;
        IsLateG2Mlike(c) = isLateG2Mlike;
        IsHighNCvolRatio(c) = isHighNCvolRatio;
        IsWithinDropletLongCycleForState3(c) = isWithinLong;
        HasState3Context(c) = hasContext;

        if isLateG2Mlike
            State_raw(c) = 3;
        elseif isHighInterphaseNuc
            State_raw(c) = 2;
        elseif isCytoLike
            State_raw(c) = 1;
        else
            State_raw(c) = NaN;
        end
    end

    nucCycleInfo.State_raw = State_raw;
    nucCycleInfo.IsCytoLike = IsCytoLike;
    nucCycleInfo.IsHighInterphaseNuc = IsHighInterphaseNuc;
    nucCycleInfo.IsLateG2Mlike = IsLateG2Mlike;
    nucCycleInfo.IsHighNCvolRatio = IsHighNCvolRatio;
    nucCycleInfo.IsWithinDropletLongCycleForState3 = IsWithinDropletLongCycleForState3;
    nucCycleInfo.HasState3Context = HasState3Context;
end

function isCytoLike = classifyState1(T, c, P)

    isCytoLike_corr = false;
    if P.useCorrForState1
        isCytoLike_corr = T.CytoDistanceCorr(c) <= P.state1CorrMultiplier * T.CytoThresholdCorr(c);
    end

    isCytoLike_rmse = false;
    if P.useRMSEForState1
        isCytoLike_rmse = T.CytoDistanceRMSE(c) <= P.state1RMSEMultiplier * T.CytoThresholdRMSE(c);
    end

    if P.state1UseOrLogic
        isCytoLike = isCytoLike_corr | isCytoLike_rmse;
    else
        isCytoLike = isCytoLike_corr & isCytoLike_rmse;
    end
end

function isHighInterphaseNuc = classifyState2(T, c, thresholds, P)

    thr = thresholds.effectiveNucRatioHighThreshold;

    isHighMeanRatio = T.InterphaseNucCytoRatioMean(c) >= thr;

    isHighP90Ratio = false;
    if P.useP90ForState2
        isHighP90Ratio = T.InterphaseNucCytoRatioP90(c) >= thr;
    end

    isHighLateInterphaseMean = false;
    if P.useLateInterphaseMeanForState2
        isHighLateInterphaseMean = T.LateInterphaseNucCytoRatioMean(c) >= thr;
    end

    isGradualNuclearIncrease = false;
    if P.useIncreaseForState2
        isGradualNuclearIncrease = T.InterphaseNucCytoRatioIncrease(c) >= P.nucRatioIncreaseThreshold;
    end

    isHighMinusCyto = false;
    if P.useNucMinusCytoForState2
        isHighMinusCyto = T.InterphaseNucMinusCytoMean(c) >= P.nucMinusCytoHighThreshold;
    end

    isHighInterphaseNuc = ...
        isHighP90Ratio | ...
        isGradualNuclearIncrease | ...
        isHighMinusCyto;
    % isHighInterphaseNuc = ...
    %     isHighMeanRatio | ...
    %     isHighP90Ratio | ...
    %     isHighLateInterphaseMean | ...
    %     isGradualNuclearIncrease | ...
    %     isHighMinusCyto;
end

function [isLateG2Mlike, isHighNCvolRatio, isWithinLong, hasContext] = classifyState3(T, c, thresholds, P)

    thr = thresholds.effectiveNucRatioHighThreshold;

    isLateTiming = T.TimeToHalfActivation_nuc(c) >= T.TroughPhase(c) + P.lateTimingOffsetFromTrough;

    isLowInterphaseForState3 = ...
        T.InterphaseNucCytoRatioP90(c) < thr & ...
        T.LateInterphaseNucCytoRatioMean(c) < thr;

    if P.useInterphaseNCvolRatioForState3
        ncForState3 = T.InterphaseNCvolRatioMean(c);
    else
        ncForState3 = T.CycleMeanNCvolRatio(c);
    end

    isHighNCvolRatio = false;
    if ~isnan(thresholds.highNCvolRatioThreshold)
        isHighNCvolRatio = ncForState3 >= thresholds.highNCvolRatioThreshold;
    end

    isWithinLong = T.IsWithinDropletLongCycle(c);

    if P.useState3ContextGate
        hasContext = isHighNCvolRatio & isWithinLong;
    else
        hasContext = true;
    end

    isLateG2Mlike = ...
        hasContext & ...
        isLowInterphaseForState3 & ...
        T.LateActivationScore_nucShape(c) >= P.lateActivationThreshold & ...
        isLateTiming;
end

function nucCycleInfo = applyOptionalStateSmoothing(nucCycleInfo, P)

    State_smooth = nucCycleInfo.State_raw;

    if P.applyConsecutiveSmoothing
        uniqueDropletIDs = unique(string(nucCycleInfo.DropletID));

        for i = 1:numel(uniqueDropletIDs)
            thisID = uniqueDropletIDs(i);
            idx = find(string(nucCycleInfo.DropletID) == thisID);
            if isempty(idx)
                continue;
            end

            [~, order] = sort(nucCycleInfo.CycleInDroplet(idx));
            idx = idx(order);

            states = nucCycleInfo.State_raw(idx);
            statesSmooth = smoothStateSequenceLeftToRight(states, P.consecutiveN);
            State_smooth(idx) = statesSmooth;
        end
    end

    nucCycleInfo.State_smooth = State_smooth;
end

function transitionTable = buildTransitionTable(nucCycleInfo)

    transitionRows = cell(0, 1);
    uniqueDropletIDs = unique(string(nucCycleInfo.DropletID));

    for i = 1:numel(uniqueDropletIDs)
        thisID = uniqueDropletIDs(i);
        idx = find(string(nucCycleInfo.DropletID) == thisID & nucCycleInfo.IsValidCycle);
        if isempty(idx)
            continue;
        end

        [~, order] = sort(nucCycleInfo.CycleInDroplet(idx));
        idx = idx(order);

        cycles = nucCycleInfo.CycleInDroplet(idx);
        states = nucCycleInfo.State_smooth(idx);

        newRow = table;
        newRow.DropletID = thisID;
        newRow.Position = nucCycleInfo.Position(idx(1));
        newRow.Droplet = nucCycleInfo.Droplet(idx(1));
        newRow.TraceRow = nucCycleInfo.TraceRow(idx(1));
        newRow.FirstState1Cycle = firstCycleOfState(cycles, states, 1);
        newRow.FirstState2Cycle = firstCycleOfState(cycles, states, 2);
        newRow.FirstState3Cycle = firstCycleOfState(cycles, states, 3);
        newRow.NumValidCycles = numel(cycles);
        newRow.MeanNCvolRatio = mean(nucCycleInfo.CycleMeanNCvolRatio(idx), "omitnan");
        newRow.MaxNCvolRatio = max(nucCycleInfo.CycleMeanNCvolRatio(idx), [], "omitnan");

        transitionRows{end+1, 1} = newRow;
    end

    if isempty(transitionRows)
        transitionTable = table;
    else
        transitionTable = vertcat(transitionRows{:});
    end

    fprintf("Transition table preview:\n");
    disp(transitionTable(1:min(10, height(transitionTable)), :));
end

%% =========================
%  Local functions: reporting / output
%  =========================

function existingCycleTable_withStates = addResultsBackToOriginalCycleTable(existingCycleTable, nucCycleInfo)

    columnsToReturn = getColumnsToReturn();
    resultTable = nucCycleInfo(:, ["OriginalCycleRow", columnsToReturn]);

    existingCycleTable_withStates = outerjoin( ...
        existingCycleTable, ...
        resultTable, ...
        "Keys", "OriginalCycleRow", ...
        "MergeKeys", true, ...
        "Type", "left");

    existingCycleTable_withStates = sortrows(existingCycleTable_withStates, "OriginalCycleRow");
end

function columnsToReturn = getColumnsToReturn()

    columnsToReturn = [ ...
        "TraceRow", "IsValidCycle", "TroughPhase", "TroughIdxNorm", ...
        "InterphaseStartIdx", "InterphaseEndIdx", "InterphaseStartFrame", "InterphaseEndFrame", ...
        "InterphaseEndReason", "LateStartIdx", "LateEndIdx", ...
        "CytoDistanceCorr", "CytoDistanceRMSE", "CytoThresholdCorr", "CytoThresholdRMSE", ...
        "InterphaseNucCytoRatioMean", "InterphaseNucCytoRatioMax", "InterphaseNucCytoRatioP90", ...
        "EarlyInterphaseNucCytoRatioMean", "LateInterphaseNucCytoRatioMean", ...
        "InterphaseNucCytoRatioIncrease", "InterphaseNucCytoRatioSlope", "InterphaseNucMinusCytoMean", ...
        "LateActivationScore_nucShape", "TimeToHalfActivation_nuc", "PeakPhase_nuc", "MaxSlopeLate_nucShape", ...
        "InterphaseNucVolMean", "InterphaseNCvolRatioMean", "CycleMeanNucVol", "CycleMeanNCvolRatio", ...
        "NucVolValidFraction", "NCvolRatioValidFraction", ...
        "InterphaseNucVolValidFraction", "InterphaseNCvolRatioValidFraction", ...
        "IsCytoLike", "IsHighInterphaseNuc", "IsLateG2Mlike", "IsHighNCvolRatio", ...
        "CycleLengthFoldChangePrev", "CycleLengthDeltaPrev", "CycleLengthZWithinDroplet", ...
        "IsWithinDropletLongCycle", "IsWithinDropletLongCycleForState3", "HasState3Context", ...
        "State_raw", "State_smooth"];
end

function makeDiagnosticPlots(nucCycleInfo, cytoRef, shapeCycles, phase, nPhasePoints, thresholds, P)

    plotStateFeatureSpace(nucCycleInfo, thresholds, P);
    plotCytoplasmicDistance(nucCycleInfo, cytoRef);
    plotStateDistribution(nucCycleInfo);
    plotState3Context(nucCycleInfo, thresholds, P);
    plotRepresentativeWaveforms(nucCycleInfo, cytoRef, shapeCycles, phase, nPhasePoints);
end

function plotCycleQC(nucCycleInfo)

    figure("Name", "Peak-to-peak cycle QC",'Visible','off');

    subplot(2, 2, 1);
    histogram(nucCycleInfo.CycleLength);
    xlabel("Cycle length, frames"); ylabel("Count");
    title("Peak-to-peak cycle length"); grid on;

    subplot(2, 2, 2);
    histogram(nucCycleInfo.TroughPhase);
    xlabel("Trough phase"); ylabel("Count");
    title("Trough position within peak-to-peak cycle"); grid on;

    subplot(2, 2, 3);
    histogram(nucCycleInfo.TotalAmplitudeMaxMinusMin);
    xlabel("Amplitude, max - min"); ylabel("Count");
    title("Cycle amplitude"); grid on;

    subplot(2, 2, 4);
    boxchart(nucCycleInfo.CycleInDroplet, nucCycleInfo.CycleLength);
    xlabel("Cycle number"); ylabel("Cycle length");
    title("Cycle length by cycle number"); grid on;

    exportgraphics(gcf, fullfile(resultDataSavePath, 'Cdk1_Classification_1.png'));
end

function plotStateFeatureSpace(T, thresholds, P)

    figure("Name", "Feature space: State 2 vs State 3",'Visible','off');
    scatter(T.InterphaseNucCytoRatioMean, T.LateActivationScore_nucShape, 30, T.CycleInDroplet, "filled");
    hold on;
    xline(thresholds.effectiveNucRatioHighThreshold, "k--", "LineWidth", 1.5);
    yline(P.lateActivationThreshold, "k--", "LineWidth", 1.5);
    xlabel("Interphase nuclear/cytoplasmic Cdk1 ratio");
    ylabel("Late activation score, nuclear Cdk1 shape");
    title("State 2 vs State 3 feature space");
    cb = colorbar; ylabel(cb, "Cycle number"); grid on;
    
    exportgraphics(gcf, fullfile(resultDataSavePath, 'Cdk1_Classification_2.png'));

    figure("Name", "Interphase enrichment metrics",'Visible','off');
    scatter(T.LateInterphaseNucCytoRatioMean, T.InterphaseNucCytoRatioIncrease, 30, T.CycleInDroplet, "filled");
    hold on;
    xline(thresholds.effectiveNucRatioHighThreshold, "k--", "LineWidth", 1.5);
    yline(P.nucRatioIncreaseThreshold, "k--", "LineWidth", 1.5);
    xlabel("Late-interphase nuclear/cytoplasmic Cdk1 ratio");
    ylabel("Late - early interphase ratio increase");
    title("Gradual interphase nuclear Cdk1 enrichment");
    cb = colorbar; ylabel(cb, "Cycle number"); grid on;
    exportgraphics(gcf, fullfile(resultDataSavePath, 'Cdk1_Classification_3.png'));


end


function plotCytoplasmicDistance(T, cytoRef)

    figure("Name", "Cytoplasmic-template distance vs nuclear enrichment",'Visible','off');
    scatter(T.CycleInDroplet, T.CytoDistanceCorr, 30, T.InterphaseNucCytoRatioMean, "filled");
    hold on;
    if ismember("ThresholdCorr_MedianPlus3MAD", string(cytoRef.summary.Properties.VariableNames))
        plot(cytoRef.summary.CycleNumber, cytoRef.summary.ThresholdCorr_MedianPlus3MAD, "k-o", "LineWidth", 1.5);
    end
    xlabel("Cycle number");
    ylabel("Distance from cytoplasmic template, 1 - Pearson r");
    title("Nuclear cycles vs peak-to-peak cytoplasmic reference");
    cb = colorbar; ylabel(cb, "Interphase nuc/cyto Cdk1 ratio"); grid on;
    exportgraphics(gcf, fullfile(resultDataSavePath, 'Cdk1_Classification_4.png'));
end

function plotStateDistribution(T)

    figure("Name", "State distribution by cycle number",'Visible','off');
    validIdx = T.IsValidCycle & ~isnan(T.State_smooth);
    boxchart(T.State_smooth(validIdx), T.CycleInDroplet(validIdx));
    xlabel("State"); ylabel("Cycle number");
    title("Cycle numbers assigned to each state"); grid on;
    exportgraphics(gcf, fullfile(resultDataSavePath, 'Cdk1_Classification_5.png'));
end

function plotState3Context(T, thresholds, P)

    figure("Name", "State 3 context gate",'Visible','off');
    scatter(T.InterphaseNCvolRatioMean, T.CycleLengthFoldChangePrev, 35, T.State_raw, "filled");
    hold on;
    if ~isnan(thresholds.highNCvolRatioThreshold)
        xline(thresholds.highNCvolRatioThreshold, "k--", "LineWidth", 1.5);
    end
    yline(P.withinDropletFoldChangeThreshold, "k--", "LineWidth", 1.5);
    xlabel("Interphase N/C volume ratio");
    ylabel("Cycle length fold-change vs previous cycle");
    title("State 3 context: high N/C OR within-droplet long cycle");
    cb = colorbar; ylabel(cb, "Raw state"); grid on;

    figure("Name", "Within-droplet cycle-length context");
    scatter(T.CycleLengthDeltaPrev, T.CycleLengthZWithinDroplet, 35, T.State_raw, "filled");
    hold on;
    xline(P.withinDropletDeltaThreshold, "k--", "LineWidth", 1.5);
    yline(P.withinDropletZThreshold, "k--", "LineWidth", 1.5);
    xlabel("Cycle length increase vs previous cycle, frames");
    ylabel("Within-droplet cycle-length robust z-score");
    title("Within-droplet long-cycle State 3 context");
    cb = colorbar; ylabel(cb, "Raw state"); grid on;
    exportgraphics(gcf, fullfile(resultDataSavePath, 'Cdk1_Classification_6.png'));
end

function plotRepresentativeWaveforms(T, cytoRef, shapeCycles, phase, nPhasePoints)

    figure("Name", "Representative normalized waveforms",'Visible','off');
    validCycleIdx = find(T.IsValidCycle);
    nShowCycles = min(12, numel(validCycleIdx));

    for ii = 1:nShowCycles
        c = validCycleIdx(ii);
        k = T.CycleInDroplet(c);

        subplot(nShowCycles, 1, ii);
        hold on;

        if k <= size(cytoRef.templateMean, 1)
            plot(phase, cytoRef.templateMean(k, :), "k-", "LineWidth", 2);
        end

        plot(phase, shapeCycles.total(c, :), "LineWidth", 1.4);
        plot(phase, shapeCycles.nuc(c, :), "LineWidth", 1.2);

        xline(T.TroughPhase(c), "--");
        if ~isnan(T.InterphaseStartIdx(c))
            xline(T.InterphaseStartIdx(c) / nPhasePoints, ":");
            xline(T.InterphaseEndIdx(c) / nPhasePoints, ":");
            xline(T.LateStartIdx(c) / nPhasePoints, "-.");
        end

        ylim([-0.05 1.05]);
        grid on;

        if ii == 1
            title("Cytoplasmic template vs nuclear-reconstituted waveforms");
            legend("Cyto template", "Total Cdk1", "Nuclear Cdk1", "Trough", "Location", "best");
        end
        if ii == nShowCycles
            xlabel("Normalized phase, peak-to-peak");
        end
    end
    exportgraphics(gcf, fullfile(resultDataSavePath, 'Cdk1_Classification_7.png'));
end

function saveClassificationOutputs(nucCycleInfo, existingCycleTable_withStates, transitionTable, cycles, normCycles, shapeCycles, phase, nPhasePoints, P, thresholds, featureContext)

    save(P.outputMatFile, ...
        "nucCycleInfo", ...
        "existingCycleTable_withStates", ...
        "transitionTable", ...
        "cycles", ...
        "normCycles", ...
        "shapeCycles", ...
        "phase", ...
        "nPhasePoints", ...
        "P", ...
        "thresholds", ...
        "featureContext");

    writetable(existingCycleTable_withStates, P.outputCsvFile);

    fprintf("Saved classification to %s\n", P.outputMatFile);
    fprintf("Saved cycle table with states to %s\n", P.outputCsvFile);
end

%% =========================
%  Local functions: generic utilities
%  =========================

function T = addStructFieldsToTable(T, S)
    names = fieldnames(S);
    for i = 1:numel(names)
        T.(names{i}) = S.(names{i});
    end
end

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
                Xsmooth(i, idx) = smoothdata(y(idx), "movmedian", smoothWindow);
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

    yInterp = interp1(oldPhase(validIdx), y(validIdx), targetPhase, "linear", NaN);
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

function q = prctileOmitNaN(y, p)

    y = y(:);
    y = y(~isnan(y));

    if isempty(y)
        q = NaN;
    else
        q = prctile(y, p);
    end
end

function statesOut = smoothStateSequenceLeftToRight(statesIn, consecutiveN)

    statesOut = statesIn(:);
    n = numel(statesOut);

    if n < consecutiveN
        return;
    end

    for i = 2:(n - 1)
        if ~isnan(statesOut(i-1)) && ~isnan(statesOut(i+1)) && ...
                statesOut(i-1) == statesOut(i+1) && ...
                statesOut(i) ~= statesOut(i-1)
            statesOut(i) = statesOut(i-1);
        end
    end

    maxStateSoFar = NaN;
    for i = 1:n
        if isnan(statesOut(i))
            continue;
        end

        if isnan(maxStateSoFar)
            maxStateSoFar = statesOut(i);
        elseif statesOut(i) < maxStateSoFar
            statesOut(i) = maxStateSoFar;
        else
            maxStateSoFar = statesOut(i);
        end
    end
end

function firstCyc = firstCycleOfState(cycles, states, stateValue)

    idx = find(states == stateValue, 1, "first");

    if isempty(idx)
        firstCyc = NaN;
    else
        firstCyc = cycles(idx);
    end
end

function idxNorm = frameToPhaseIndex(frame, peak1Frame, peak2Frame, nPhasePoints)

    if isnan(frame) || peak2Frame <= peak1Frame
        idxNorm = NaN;
        return;
    end

    phaseValue = (frame - peak1Frame) / (peak2Frame - peak1Frame);
    idxNorm = round(phaseValue * (nPhasePoints - 1)) + 1;
    idxNorm = max(1, min(nPhasePoints, idxNorm));
end

function frame = phaseIndexToFrame(idxNorm, peak1Frame, peak2Frame, nPhasePoints)

    if isnan(idxNorm) || nPhasePoints <= 1 || peak2Frame <= peak1Frame
        frame = NaN;
        return;
    end

    phaseValue = (idxNorm - 1) / (nPhasePoints - 1);
    frame = round(peak1Frame + phaseValue * (peak2Frame - peak1Frame));
    frame = max(peak1Frame, min(peak2Frame, frame));
end

function firstIdx = firstConsecutiveTrue(logicalVec, consecutiveN)

    logicalVec = logicalVec(:)';
    firstIdx = NaN;

    if isempty(logicalVec) || consecutiveN < 1 || numel(logicalVec) < consecutiveN
        return;
    end

    for i = 1:(numel(logicalVec) - consecutiveN + 1)
        if all(logicalVec(i:(i + consecutiveN - 1)))
            firstIdx = i;
            return;
        end
    end
end
end