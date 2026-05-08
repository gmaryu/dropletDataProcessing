function dataSet = calcCycleAdditionalFeatures(dataSet)
% calcCycleAdditionalFeatures
% Adds cycle-level nuclear-growth and DNA-replication features to dataSet.cycle.
%
% Main changes from the older version:
%   1. Calculations are organized by droplet/track, not by cycle row.
%      This avoids repeatedly extracting the same time-series table.
%   2. Nuclear growth is quantified from NUC_VOLUMEUM3 using shape-aware metrics:
%      overall slope, early slope, late slope, plateau index, late gain fraction,
%      time-to-max fraction, piecewise improvement, and growth type.
%   3. DNA replication is quantified using Q90-to-previous-Q90 fold change,
%      plus robust within-cycle start/end medians and slopes.
%   4. The legacy typo SUM_NUCLEUS_HORCHST_INT_MOD is aliased to the corrected
%      name SUM_NUCLEUS_HOECHST_INT_MOD when needed.

    dnarenormfactor = 1e7; % scaling factor for DNA summed intensity-derived quantities

    % --- Extract variables ---
    tp = dataSet.cycle;
    tm = dataSet.timeSeries;
    FrameToMin = dataSet.FrameToMin;
    PixelToUm = dataSet.PixelToUm;

    % --- Add IGNORE columns if necessary ---
    if ~ismember("IGNORED", dataSet.info.Properties.VariableNames)
        disp("IGNORED is added to dataSet.info");
        dataSet.info.IGNORED = zeros(height(dataSet.info), 1);
    end

    if ~ismember("IGNORED", tp.Properties.VariableNames)
        disp("IGNORED is added to dataSet.cycle");
        tp.IGNORED = zeros(height(tp), 1);
    end

    if ~ismember("IGNORED", tm.Properties.VariableNames)
        disp("IGNORED is added to dataSet.timeSeries");
        tm.IGNORED = zeros(height(tm), 1);
    end

    % --- Fix / alias Hoechst typo in timeSeries ---
    % The older pipeline generated SUM_NUCLEUS_HORCHST_INT_MOD.
    % Use the corrected name internally while preserving the old column.
    if ismember('SUM_NUCLEUS_HORCHST_INT_MOD', tm.Properties.VariableNames) && ...
       ~ismember('SUM_NUCLEUS_HOECHST_INT_MOD', tm.Properties.VariableNames)
        disp('Aliased legacy SUM_NUCLEUS_HORCHST_INT_MOD to SUM_NUCLEUS_HOECHST_INT_MOD');
        tm.SUM_NUCLEUS_HOECHST_INT_MOD = tm.SUM_NUCLEUS_HORCHST_INT_MOD;
    end

    nCycles = height(tp);

    % ---------------------------------------------------------------------
    % Allocate new / updated cycle-level variables
    % ---------------------------------------------------------------------
    nucVolOverallSlope      = nan(nCycles, 1); % um^3 / min
    nucVolEarlySlope        = nan(nCycles, 1); % um^3 / min
    nucVolLateSlope         = nan(nCycles, 1); % um^3 / min
    nucVolPlateauIndex      = nan(nCycles, 1);
    nucVolGrowthFraction    = nan(nCycles, 1);
    nucVolPiecewiseImprove  = nan(nCycles, 1);
    nucVolPiecewiseBreakIdx = nan(nCycles, 1);
    nucVolPiecewiseSlope1   = nan(nCycles, 1); % um^3 / min
    nucVolPiecewiseSlope2   = nan(nCycles, 1); % um^3 / min
    nucVolTotalGrowth       = nan(nCycles, 1); % um^3
    nucVolDynamicRange      = nan(nCycles, 1); % um^3
    nucVolEarlyGainFraction = nan(nCycles, 1);
    nucVolLateGainFraction  = nan(nCycles, 1);
    nucVolTimeToMaxFraction = nan(nCycles, 1);
    nucVolEarlySlopeNorm    = nan(nCycles, 1);
    nucVolLateSlopeNorm     = nan(nCycles, 1);
    nucVolStartValue        = nan(nCycles, 1); % um^3
    nucVolEndValue          = nan(nCycles, 1); % um^3
    nucVolMaxValue          = nan(nCycles, 1); % um^3
    nucVolDurationMin       = nan(nCycles, 1);
    nucVolGrowthType        = strings(nCycles, 1);
    nucVolGrowthType(:)     = missing;

    nucNpixelMODQ90         = nan(nCycles, 1);
    nucDNAintMODQ90         = nan(nCycles, 1);
    nucDNAintDNAQ90         = nan(nCycles, 1);
    nucSurfQ90              = nan(nCycles, 1);
    nucVolQ90               = nan(nCycles, 1);

    dnaIncRateCoeff         = nan(nCycles, 1); % scaled a.u. / min
    dnaDeltaCycle           = nan(nCycles, 1); % scaled a.u.
    dnaFoldChangeCycle      = nan(nCycles, 1); % current-cycle Q90 / previous-cycle Q90
    dnaFoldChangeWithin     = nan(nCycles, 1); % end median / start median within the same cycle
    dnaPrevQ90              = nan(nCycles, 1); % previous-cycle Q90 used as denominator
    dnaCurrentQ90           = nan(nCycles, 1); % current-cycle Q90 used as numerator
    dnaRateCycle            = nan(nCycles, 1); % scaled a.u. / min
    dnaStartMedian          = nan(nCycles, 1); % raw summed intensity
    dnaEndMedian            = nan(nCycles, 1); % raw summed intensity
    dnaDurationMin          = nan(nCycles, 1);

    % ---------------------------------------------------------------------
    % Only do nuclear/DNA feature calculation when relevant time-series
    % columns exist. The function still computes general cycle features later.
    % ---------------------------------------------------------------------
    hasNucVolume = ismember('NUC_VOLUMEUM3', tm.Properties.VariableNames);
    hasNucFeature = ismember('NUC_NPIXELS_Q90', tp.Properties.VariableNames) && ...
                    any(isfinite(tp.NUC_NPIXELS_Q90));

    if hasNucFeature || hasNucVolume

        % Group cycle rows by POS_ID and TRACK_ID. This is the main speed-up.
        validRows = find(~logical(tp.IGNORED));

        if ~isempty(validRows)
            groupKey = string(tp.POS_ID(validRows)) + "_" + string(tp.TRACK_ID(validRows));
            uniqueGroups = unique(groupKey, 'stable');

            for g = 1:numel(uniqueGroups)
                rowsInGroup = validRows(groupKey == uniqueGroups(g));
                if isempty(rowsInGroup)
                    continue
                end

                % Ensure cycle-to-cycle calculations, especially DNA_FC_CYCLE,
                % use chronological order within each droplet/track.
                [~, rowOrder] = sort(tp.START_FRAME(rowsInGroup));
                rowsInGroup = rowsInGroup(rowOrder);

                posid = tp.POS_ID(rowsInGroup(1));
                trackid = tp.TRACK_ID(rowsInGroup(1));

                tmp_tm = tm(tm.POS_ID == posid & tm.TRACK_ID == trackid, :);
                if isempty(tmp_tm)
                    continue
                end
                tmp_tm = sortrows(tmp_tm, 'FRAME');

                nLocalCycles = numel(rowsInGroup);
                cycleStartIdx = nan(nLocalCycles, 1);
                cycleEndIdx   = nan(nLocalCycles, 1);

                % Convert cycle frame boundaries into local time-series indices.
                for j = 1:nLocalCycles
                    k = rowsInGroup(j);
                    [startFrame, endFrame] = getCycleFrameWindow(tp, k);

                    startidx = find(tmp_tm.FRAME == startFrame, 1, 'first');
                    endidx   = find(tmp_tm.FRAME == endFrame,   1, 'first');

                    if isempty(startidx) || isempty(endidx) || endidx < startidx
                        continue
                    end

                    cycleStartIdx(j) = startidx;
                    cycleEndIdx(j) = endidx;
                end

                % Nuclear volume growth-shape quantification.
                if hasNucVolume
                    if ismember('NPIXEL_NUC_MOD', tmp_tm.Properties.VariableNames)
                        nucleusDetected = tmp_tm.NPIXEL_NUC_MOD > 0 & isfinite(tmp_tm.NUC_VOLUMEUM3);
                    elseif ismember('NPIXEL_NUC', tmp_tm.Properties.VariableNames)
                        nucleusDetected = tmp_tm.NPIXEL_NUC > 0 & isfinite(tmp_tm.NUC_VOLUMEUM3);
                    else
                        nucleusDetected = isfinite(tmp_tm.NUC_VOLUMEUM3);
                    end

                    growthResult = quantifyNuclearGrowthShape( ...
                        tmp_tm.NUC_VOLUMEUM3, ...
                        nucleusDetected, ...
                        cycleStartIdx, ...
                        cycleEndIdx, ...
                        FrameToMin);

                    for j = 1:nLocalCycles
                        k = rowsInGroup(j);
                        if j > height(growthResult)
                            continue
                        end

                        nucVolOverallSlope(k)      = growthResult.overallSlope(j);
                        nucVolEarlySlope(k)        = growthResult.earlySlope(j);
                        nucVolLateSlope(k)         = growthResult.lateSlope(j);
                        nucVolPlateauIndex(k)      = growthResult.plateauIndex(j);
                        nucVolGrowthFraction(k)    = growthResult.growthFraction(j);
                        nucVolPiecewiseImprove(k)  = growthResult.piecewiseImprovement(j);
                        nucVolPiecewiseBreakIdx(k) = growthResult.piecewiseBreakIdx(j);
                        nucVolPiecewiseSlope1(k)   = growthResult.piecewiseEarlySlope(j);
                        nucVolPiecewiseSlope2(k)   = growthResult.piecewiseLateSlope(j);
                        nucVolTotalGrowth(k)       = growthResult.totalGrowth(j);
                        nucVolDynamicRange(k)      = growthResult.dynamicRange(j);
                        nucVolEarlyGainFraction(k) = growthResult.earlyGainFraction(j);
                        nucVolLateGainFraction(k)  = growthResult.lateGainFraction(j);
                        nucVolTimeToMaxFraction(k) = growthResult.timeToMaxFraction(j);
                        nucVolEarlySlopeNorm(k)    = growthResult.earlySlopeNorm(j);
                        nucVolLateSlopeNorm(k)     = growthResult.lateSlopeNorm(j);
                        nucVolStartValue(k)        = growthResult.startValue(j);
                        nucVolEndValue(k)          = growthResult.endValue(j);
                        nucVolMaxValue(k)          = growthResult.maxValue(j);
                        nucVolDurationMin(k)       = growthResult.durationMin(j);
                        nucVolGrowthType(k)        = growthResult.growthType(j);
                    end
                end

                % Per-cycle Q90 and DNA replication features.
                for j = 1:nLocalCycles
                    k = rowsInGroup(j);
                    if ~isfinite(cycleStartIdx(j)) || ~isfinite(cycleEndIdx(j))
                        continue
                    end

                    cycleData = tmp_tm(cycleStartIdx(j):cycleEndIdx(j), :);

                    nucNpixelMODQ90(k) = safeQuantileFromTable(cycleData, 'NPIXEL_NUC_MOD', 0.9);
                    nucDNAintMODQ90(k) = safeQuantileFromTable(cycleData, 'SUM_NUCLEUS_HOECHST_INT_MOD', 0.9);
                    nucDNAintDNAQ90(k) = safeQuantileFromTable(cycleData, 'SUM_SPERM_HOECHST_INT', 0.9);
                    nucSurfQ90(k)      = safeQuantileFromTable(cycleData, 'NUC_SURF_AREA', 0.9);
                    nucVolQ90(k)       = safeQuantileFromTable(cycleData, 'NUC_VOLUMEUM3', 0.9);

                    % DNA replication metrics. Prefer the corrected MOD signal.
                    dnaSignalName = chooseFirstExistingVariable(cycleData, { ...
                        'SUM_NUCLEUS_HOECHST_INT_MOD', ...
                        'SUM_NUCLEUS_HOECHST_INT', ...
                        'SUM_SPERM_HOECHST_INT'});

                    if strlength(dnaSignalName) > 0
                        dnaStats = quantifyDNAReplicationCycle( ...
                            cycleData.FRAME, ...
                            cycleData.(dnaSignalName), ...
                            FrameToMin);

                        dnaIncRateCoeff(k)    = dnaStats.slope / dnarenormfactor;
                        dnaDeltaCycle(k)      = dnaStats.delta / dnarenormfactor;
                        dnaFoldChangeWithin(k) = dnaStats.foldChange;
                        dnaRateCycle(k)       = dnaStats.rate / dnarenormfactor;
                        dnaStartMedian(k)     = dnaStats.startValue;
                        dnaEndMedian(k)       = dnaStats.endValue;
                        dnaDurationMin(k)     = dnaStats.durationMin;
                    end
                end

                % DNA fold change based on cycle Q90 values.
                % This preserves the older logic: current cycle Q90 / previous cycle Q90
                % within the same droplet/track. This is more robust than using the
                % first frame of the current cycle, which can be artificially low due
                % to chromosome condensation, segmentation, or focus effects.
                dnaQ90ForFC = chooseDNAQ90ForFoldChange( ...
                    nucDNAintMODQ90(rowsInGroup), ...
                    nucDNAintDNAQ90(rowsInGroup));

                for j = 1:nLocalCycles
                    k = rowsInGroup(j);
                    dnaCurrentQ90(k) = dnaQ90ForFC(j);
                    if j > 1
                        prevQ90 = dnaQ90ForFC(j-1);
                        dnaPrevQ90(k) = prevQ90;
                        if isfinite(dnaQ90ForFC(j)) && isfinite(prevQ90) && prevQ90 > 0
                            dnaFoldChangeCycle(k) = dnaQ90ForFC(j) / prevQ90;
                        end
                    end
                end
            end
        end

        % -----------------------------------------------------------------
        % Write nuclear and DNA features back to the cycle table.
        % Existing columns are intentionally updated to keep features
        % consistent with the current implementation.
        % -----------------------------------------------------------------
        tp.NUC_INC_RATE_COEFF = nucVolOverallSlope;
        tp.NUC_VOL_OVERALL_SLOPE = nucVolOverallSlope;
        tp.NUC_VOL_EARLY_SLOPE = nucVolEarlySlope;
        tp.NUC_VOL_LATE_SLOPE = nucVolLateSlope;
        tp.NUC_VOL_PLATEAU_INDEX = nucVolPlateauIndex;
        tp.NUC_VOL_GROWTH_FRACTION = nucVolGrowthFraction;
        tp.NUC_VOL_PIECEWISE_IMPROVEMENT = nucVolPiecewiseImprove;
        tp.NUC_VOL_PIECEWISE_BREAK_IDX = nucVolPiecewiseBreakIdx;
        tp.NUC_VOL_PIECEWISE_EARLY_SLOPE = nucVolPiecewiseSlope1;
        tp.NUC_VOL_PIECEWISE_LATE_SLOPE = nucVolPiecewiseSlope2;
        tp.NUC_VOL_TOTAL_GROWTH = nucVolTotalGrowth;
        tp.NUC_VOL_DYNAMIC_RANGE = nucVolDynamicRange;
        tp.NUC_VOL_EARLY_GAIN_FRACTION = nucVolEarlyGainFraction;
        tp.NUC_VOL_LATE_GAIN_FRACTION = nucVolLateGainFraction;
        tp.NUC_VOL_TIME_TO_MAX_FRACTION = nucVolTimeToMaxFraction;
        tp.NUC_VOL_EARLY_SLOPE_NORM = nucVolEarlySlopeNorm;
        tp.NUC_VOL_LATE_SLOPE_NORM = nucVolLateSlopeNorm;
        tp.NUC_VOL_START_VALUE = nucVolStartValue;
        tp.NUC_VOL_END_VALUE = nucVolEndValue;
        tp.NUC_VOL_MAX_VALUE = nucVolMaxValue;
        tp.NUC_VOL_DURATION_MIN = nucVolDurationMin;
        tp.NUC_VOL_GROWTH_TYPE = nucVolGrowthType;

        tp.NUC_NPIXELS_MOD_Q90 = nucNpixelMODQ90;
        tp.DNA_SUM_INT_MOD_Q90 = nucDNAintMODQ90;
        tp.DNA_SUM_INT_DNA_Q90 = nucDNAintDNAQ90;
        tp.NSURF_Q90 = nucSurfQ90;
        tp.NUC_VOLUMEUM3_Q90 = nucVolQ90;

        tp.DNA_INC_RATE_COEFF = dnaIncRateCoeff;
        tp.DNA_DELTA_CYCLE = dnaDeltaCycle;
        tp.DNA_FC_CYCLE = dnaFoldChangeCycle;
        tp.DNA_FC_WITHIN_CYCLE = dnaFoldChangeWithin;
        tp.DNA_PREV_Q90_FOR_FC = dnaPrevQ90;
        tp.DNA_CURRENT_Q90_FOR_FC = dnaCurrentQ90;
        tp.DNA_RATE_CYCLE = dnaRateCycle;
        tp.DNA_START_MEDIAN = dnaStartMedian;
        tp.DNA_END_MEDIAN = dnaEndMedian;
        tp.DNA_DURATION_MIN = dnaDurationMin;
    end

    % ---------------------------------------------------------------------
    % General cycle features
    % ---------------------------------------------------------------------
    tp.START_MINUTE = tp.START_FRAME * FrameToMin;
    tp.DURATION = (tp.END_FRAME - tp.START_FRAME) * FrameToMin;

    if ismember('NUC_NPIXELS_Q90', tp.Properties.VariableNames) && any(isfinite(tp.NUC_NPIXELS_Q90))
        dropletVolume = visualization.convertAreaPixelsToVolume(tp.AREA_NPIXELS_MEDIAN, PixelToUm);

        tp.NCVR_ORI = power(tp.NUC_NPIXELS_Q90 ./ tp.AREA_NPIXELS_MEDIAN, 3 / 2);

        if ismember('NUC_NPIXELS_MOD_Q90', tp.Properties.VariableNames)
            tp.NCVR = power(tp.NUC_NPIXELS_MOD_Q90 ./ tp.AREA_NPIXELS_MEDIAN, 3 / 2);
        end

        if ismember('DNA_SUM_INT_Q90', tp.Properties.VariableNames)
            tp.DNACR_ORI = tp.DNA_SUM_INT_Q90 ./ dropletVolume / dnarenormfactor;
        end
        if ismember('DNA_SUM_INT_MOD_Q90', tp.Properties.VariableNames)
            tp.DNACR_NUC = tp.DNA_SUM_INT_MOD_Q90 ./ dropletVolume / dnarenormfactor;
        end
        if ismember('DNA_SUM_INT_DNA_Q90', tp.Properties.VariableNames)
            tp.DNACR_DNA = tp.DNA_SUM_INT_DNA_Q90 ./ dropletVolume / dnarenormfactor;
        end

        % Keep the existing foldChangeDNA calculation if available.
        tp.FC_DNA = visualization.foldChangeDNA(tm, tp);
    end

    tp.FC_Period = visualization.foldChangePeriod(tp);
    tp.VOLUMEUM3 = visualization.convertAreaPixelsToVolume(tp.AREA_NPIXELS_MEDIAN, PixelToUm);
    tp.MARKERSIZE = (log10(tp.VOLUMEUM3) - 5.25) * 5; % log um3 volume roughly within 5-7
    tp.MARKERSIZE(tp.MARKERSIZE < 0.5, :) = 0.5;

    dataSet.cycle = tp;
    dataSet.timeSeries = tm;
end

% ========================================================================
% Local helper functions
% ========================================================================

function [startFrame, endFrame] = getCycleFrameWindow(tp, k)
% Use interphase window when available; otherwise use the full Cdk1-defined cycle.

    if ismember('INTERPHASE_START_FRAME', tp.Properties.VariableNames) && ...
       isfinite(tp.INTERPHASE_START_FRAME(k))
        startFrame = tp.INTERPHASE_START_FRAME(k);
    else
        startFrame = tp.START_FRAME(k);
    end

    if ismember('INTERPHASE_END_FRAME', tp.Properties.VariableNames) && ...
       isfinite(tp.INTERPHASE_END_FRAME(k))
        endFrame = tp.INTERPHASE_END_FRAME(k);
    else
        endFrame = tp.END_FRAME(k);
    end
end

function result = quantifyNuclearGrowthShape(y, nucleusDetected, cycleStartIdx, cycleEndIdx, frameToMin)
% quantifyNuclearGrowthShape
% Quantifies nuclear-volume growth shape for multiple cycles at once.
%
% The function returns one row per input cycle. Short or invalid cycles are
% retained as rows with NaNs and growthType = "tooShort" or "invalid".

    y = y(:);
    nucleusDetected = logical(nucleusDetected(:));
    nCycles = numel(cycleStartIdx);

    % Parameters. These can be tuned using control data.
    smoothWindow = 3;
    earlyFraction = 0.30;
    lateFraction = 0.50;
    minOverallPoints = 2;
    minRobustPoints = 5;
    minClassifyPoints = 8;
    minSegmentLength = 5;
    minPiecewiseImprovement = 0.25;

    ySmoothAll = smoothdata(y, 'movmedian', smoothWindow, 'omitnan');

    result = table();
    result.cycle = (1:nCycles)';
    result.startIdx = nan(nCycles, 1);
    result.endIdx = nan(nCycles, 1);
    result.durationMin = nan(nCycles, 1);
    result.overallSlope = nan(nCycles, 1);
    result.earlySlope = nan(nCycles, 1);
    result.lateSlope = nan(nCycles, 1);
    result.plateauIndex = nan(nCycles, 1);
    result.growthFraction = nan(nCycles, 1);
    result.piecewiseImprovement = nan(nCycles, 1);
    result.piecewiseBreakIdx = nan(nCycles, 1);
    result.piecewiseEarlySlope = nan(nCycles, 1);
    result.piecewiseLateSlope = nan(nCycles, 1);
    result.totalGrowth = nan(nCycles, 1);
    result.dynamicRange = nan(nCycles, 1);
    result.earlyGainFraction = nan(nCycles, 1);
    result.lateGainFraction = nan(nCycles, 1);
    result.timeToMaxFraction = nan(nCycles, 1);
    result.earlySlopeNorm = nan(nCycles, 1);
    result.lateSlopeNorm = nan(nCycles, 1);
    result.startValue = nan(nCycles, 1);
    result.endValue = nan(nCycles, 1);
    result.maxValue = nan(nCycles, 1);
    result.growthType = strings(nCycles, 1);
    result.growthType(:) = "invalid";

    for c = 1:nCycles
        if ~isfinite(cycleStartIdx(c)) || ~isfinite(cycleEndIdx(c)) || cycleEndIdx(c) < cycleStartIdx(c)
            continue
        end

        idxCycle = cycleStartIdx(c):cycleEndIdx(c);
        idxCycle = idxCycle(idxCycle >= 1 & idxCycle <= numel(y));
        idxNuc = idxCycle(nucleusDetected(idxCycle) & isfinite(y(idxCycle)));
        idxNuc = keepLongestContinuousSegment(idxNuc);

        if numel(idxNuc) < minOverallPoints
            result.growthType(c) = "tooShort";
            continue
        end

        yRaw = y(idxNuc);
        ySmooth = ySmoothAll(idxNuc);
        tau = (idxNuc(:) - idxNuc(1)) * frameToMin;
        durationMin = tau(end) - tau(1);

        result.startIdx(c) = idxNuc(1);
        result.endIdx(c) = idxNuc(end);
        result.durationMin(c) = durationMin;

        result.overallSlope(c) = safeSlopeFit(tau, yRaw, minRobustPoints);

        nPts = numel(idxNuc);
        if nPts < 4 || durationMin <= 0
            result.growthType(c) = "tooShort";
            continue
        end

        nEarly = max(2, round(earlyFraction * nPts));
        nLate  = max(2, round(lateFraction  * nPts));
        nEarly = min(nEarly, floor(nPts / 2));
        nLate  = min(nLate,  floor(nPts / 2));

        earlyIdxLocal = 1:nEarly;
        lateIdxLocal = (nPts - nLate + 1):nPts;

        earlySlope = safeSlopeFit(tau(earlyIdxLocal), yRaw(earlyIdxLocal), minRobustPoints);
        lateSlope  = safeSlopeFit(tau(lateIdxLocal),  yRaw(lateIdxLocal),  minRobustPoints);

        result.earlySlope(c) = earlySlope;
        result.lateSlope(c) = lateSlope;

        [pieceImprove, pieceBreakIdx, pieceSlope1, pieceSlope2] = ...
            bestTwoPhaseFit(tau, yRaw, idxNuc, minSegmentLength, minRobustPoints);

        result.piecewiseImprovement(c) = pieceImprove;
        result.piecewiseBreakIdx(c) = pieceBreakIdx;
        result.piecewiseEarlySlope(c) = pieceSlope1;
        result.piecewiseLateSlope(c) = pieceSlope2;

        % Shape metrics based mainly on smoothed trajectory, while slopes use raw data.
        startValue = median(ySmooth(1:min(3, nPts)), 'omitnan');
        endValue = median(ySmooth(max(1, nPts-2):nPts), 'omitnan');
        maxValue = max(ySmooth, [], 'omitnan');
        minValue = min(ySmooth, [], 'omitnan');

        totalGrowth = maxValue - startValue;
        dynamicRange = maxValue - minValue;

        result.startValue(c) = startValue;
        result.endValue(c) = endValue;
        result.maxValue(c) = maxValue;
        result.totalGrowth(c) = totalGrowth;
        result.dynamicRange(c) = dynamicRange;

        if ~(isfinite(totalGrowth) && totalGrowth > 0 && isfinite(dynamicRange) && dynamicRange > 0)
            result.growthType(c) = "noGrowth";
            continue
        end

        earlyEndLocal = max(2, round(0.30 * nPts));
        midLocal = max(2, round(0.50 * nPts));
        earlyEndLocal = min(earlyEndLocal, nPts);
        midLocal = min(midLocal, nPts);

        earlyEndValue = median(ySmooth(max(1, earlyEndLocal-1):min(nPts, earlyEndLocal+1)), 'omitnan');
        midValue = median(ySmooth(max(1, midLocal-1):min(nPts, midLocal+1)), 'omitnan');

        earlyGain = earlyEndValue - startValue;
        lateGain = endValue - midValue;

        earlyGainFraction = earlyGain / totalGrowth;
        lateGainFraction = lateGain / totalGrowth;

        [~, maxLocalIdx] = max(ySmooth);
        timeToMaxFraction = maxLocalIdx / nPts;

        overallGrowthRate = totalGrowth / durationMin;
        if overallGrowthRate > 0
            earlySlopeNorm = earlySlope / overallGrowthRate;
            lateSlopeNorm = lateSlope / overallGrowthRate;
        else
            earlySlopeNorm = NaN;
            lateSlopeNorm = NaN;
        end

        earlySlopePos = max(earlySlope, eps);
        lateSlopePos = max(lateSlope, 0);
        plateauIndex = 1 - lateSlopePos / earlySlopePos;
        plateauIndex = max(0, min(1, plateauIndex));

        dy = diff(ySmooth) / frameToMin;
        derivativeThreshold = 0.05 * dynamicRange / max(durationMin, eps);
        growthFraction = mean(dy > derivativeThreshold, 'omitnan');

        result.plateauIndex(c) = plateauIndex;
        result.growthFraction(c) = growthFraction;
        result.earlyGainFraction(c) = earlyGainFraction;
        result.lateGainFraction(c) = lateGainFraction;
        result.timeToMaxFraction(c) = timeToMaxFraction;
        result.earlySlopeNorm(c) = earlySlopeNorm;
        result.lateSlopeNorm(c) = lateSlopeNorm;

        if nPts < minClassifyPoints
            result.growthType(c) = "tooShort";
            continue
        end

        % Conservative classification:
        % plateau requires weak late growth AND little gain after the mid-point
        % AND early arrival at maximum. This avoids calling ordinary late
        % saturation a true plateau.
        isPlateau = ...
            isfinite(lateSlopeNorm) && isfinite(lateGainFraction) && isfinite(timeToMaxFraction) && ...
            lateSlopeNorm < 0.15 && ...
            lateGainFraction < 0.20 && ...
            timeToMaxFraction < 0.65;

        isTwoPhase = ...
            isfinite(pieceImprove) && isfinite(earlySlopeNorm) && isfinite(lateSlopeNorm) && ...
            pieceImprove > minPiecewiseImprovement && ...
            earlySlopeNorm > 1.5 && ...
            lateSlopeNorm > 0.15 && ...
            lateGainFraction >= 0.20;

        isLateSaturation = ...
            isfinite(lateSlopeNorm) && isfinite(lateGainFraction) && ...
            lateSlopeNorm < 0.15 && ...
            lateGainFraction >= 0.20;

        if isPlateau
            result.growthType(c) = "plateau";
        elseif isTwoPhase
            result.growthType(c) = "twoPhase";
        elseif isLateSaturation
            result.growthType(c) = "lateSaturation";
        else
            result.growthType(c) = "continuous";
        end
    end
end

function stats = quantifyDNAReplicationCycle(frame, dnaSignal, frameToMin)
% quantifyDNAReplicationCycle
% Quantifies DNA replication within one cycle using robust edge medians and slope.

    frame = frame(:);
    dnaSignal = dnaSignal(:);
    valid = isfinite(frame) & isfinite(dnaSignal);
    frame = frame(valid);
    dnaSignal = dnaSignal(valid);

    stats.slope = NaN;
    stats.startValue = NaN;
    stats.endValue = NaN;
    stats.delta = NaN;
    stats.foldChange = NaN;
    stats.rate = NaN;
    stats.durationMin = NaN;

    if numel(dnaSignal) < 2
        return
    end

    t = (frame - frame(1)) * frameToMin;
    stats.durationMin = t(end) - t(1);
    stats.slope = safeSlopeFit(t, dnaSignal, 5);

    n = numel(dnaSignal);
    nEdge = max(2, round(0.20 * n));
    nEdge = min(nEdge, n);

    stats.startValue = median(dnaSignal(1:nEdge), 'omitnan');
    stats.endValue = median(dnaSignal((n-nEdge+1):n), 'omitnan');
    stats.delta = stats.endValue - stats.startValue;

    if isfinite(stats.startValue) && stats.startValue > 0
        stats.foldChange = stats.endValue / stats.startValue;
    end

    if isfinite(stats.durationMin) && stats.durationMin > 0
        stats.rate = stats.delta / stats.durationMin;
    end
end

function slope = safeSlopeFit(t, y, minRobustPoints)
% safeSlopeFit
% Uses robustfit when possible, falls back to polyfit for short segments,
% and returns NaN when the slope is not estimable.

    if nargin < 3
        minRobustPoints = 5;
    end

    t = t(:);
    y = y(:);
    valid = isfinite(t) & isfinite(y);
    t = t(valid);
    y = y(valid);

    if numel(y) < 2 || numel(unique(t)) < 2
        slope = NaN;
        return
    end

    if numel(y) >= minRobustPoints && exist('robustfit', 'file') == 2
        try
            b = robustfit(t, y);
            slope = b(2);
            return
        catch
            % Fall through to polyfit.
        end
    end

    try
        p = polyfit(t, y, 1);
        slope = p(1);
    catch
        slope = NaN;
    end
end

function [bestImprovement, bestBreakIdx, bestSlope1, bestSlope2] = bestTwoPhaseFit(t, y, originalIdx, minSegmentLength, minRobustPoints)
% bestTwoPhaseFit
% Finds the best two-segment linear fit and reports RSS improvement over a
% single linear fit.

    t = t(:);
    y = y(:);
    n = numel(y);

    bestImprovement = NaN;
    bestBreakIdx = NaN;
    bestSlope1 = NaN;
    bestSlope2 = NaN;

    if n < 2 * minSegmentLength
        return
    end

    singleSlope = safeSlopeFit(t, y, minRobustPoints);
    if ~isfinite(singleSlope)
        return
    end

    pSingle = polyfit(t, y, 1);
    yPredSingle = polyval(pSingle, t);
    RSSsingle = sum((y - yPredSingle).^2, 'omitnan');

    if ~(isfinite(RSSsingle) && RSSsingle > 0)
        return
    end

    bestRSS = Inf;

    for b = minSegmentLength:(n - minSegmentLength)
        t1 = t(1:b);
        y1 = y(1:b);
        t2 = t((b+1):end);
        y2 = y((b+1):end);

        slope1 = safeSlopeFit(t1, y1, minRobustPoints);
        slope2 = safeSlopeFit(t2, y2, minRobustPoints);
        if ~isfinite(slope1) || ~isfinite(slope2)
            continue
        end

        p1 = polyfit(t1, y1, 1);
        p2 = polyfit(t2, y2, 1);
        yPred1 = polyval(p1, t1);
        yPred2 = polyval(p2, t2);
        RSS = sum((y1 - yPred1).^2, 'omitnan') + sum((y2 - yPred2).^2, 'omitnan');

        if RSS < bestRSS
            bestRSS = RSS;
            bestBreakIdx = originalIdx(b);
            bestSlope1 = slope1;
            bestSlope2 = slope2;
        end
    end

    if isfinite(bestRSS)
        bestImprovement = (RSSsingle - bestRSS) / RSSsingle;
    end
end

function idxOut = keepLongestContinuousSegment(idxIn)
% keepLongestContinuousSegment
% Keeps the longest consecutive index segment.

    idxIn = idxIn(:);
    if isempty(idxIn)
        idxOut = idxIn;
        return
    end

    breaks = find(diff(idxIn) > 1);
    segStarts = [1; breaks + 1];
    segEnds = [breaks; numel(idxIn)];
    segLengths = segEnds - segStarts + 1;
    [~, bestSeg] = max(segLengths);
    idxOut = idxIn(segStarts(bestSeg):segEnds(bestSeg));
end

function q = safeQuantileFromTable(T, varName, prob)
% safeQuantileFromTable
% Returns quantile while handling missing variables and all-NaN vectors.

    q = NaN;
    if ~ismember(varName, T.Properties.VariableNames)
        return
    end

    x = T.(varName);
    x = x(isfinite(x));
    if isempty(x)
        return
    end

    q = quantile(x, prob);
end


function dnaQ90 = chooseDNAQ90ForFoldChange(dnaModQ90, dnaDnaQ90)
% chooseDNAQ90ForFoldChange
% Chooses the Q90 DNA signal used for cycle-to-cycle fold change.
%
% Preferred signal:
%   1. SUM_NUCLEUS_HOECHST_INT_MOD-derived Q90
%   2. SUM_SPERM_HOECHST_INT-derived Q90
%
% The returned vector has the same length as the input cycle group.

    dnaQ90 = dnaModQ90(:);

    if nargin >= 2 && ~isempty(dnaDnaQ90)
        dnaDnaQ90 = dnaDnaQ90(:);
        replaceIdx = ~isfinite(dnaQ90) & isfinite(dnaDnaQ90);
        dnaQ90(replaceIdx) = dnaDnaQ90(replaceIdx);
    end
end

function varName = chooseFirstExistingVariable(T, candidates)
% chooseFirstExistingVariable
% Returns the first variable name in candidates that exists in table T.

    varName = "";
    for i = 1:numel(candidates)
        if ismember(candidates{i}, T.Properties.VariableNames)
            varName = string(candidates{i});
            return
        end
    end
end
